     '********       PROGRAM QUAD        **********
     '*     2-D STRESS ANALYSIS USING 4-NODE      *
     '*  QUADRILATERAL ELEMENTS WITH TEMPERATURE  *
     '*    T.R.Chandrupatla and A.D.Belegundu     *
     '*********************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "FEA - 4 NODED QUADRILATERAL"; SPACE$(24);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     PRINT "  1) Plane Stress Analysis"
     PRINT "  2) Plane Strain Analysis"
     INPUT "     Choose 1 or 2 "; LC
     IF LC < 1 OR LC > 2 THEN LC = 1   '--- default is Plane Stress
     PRINT : INPUT "Data File Name <path fn.ext>  ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 3   'Material properties E, Nu, Alpha
     INPUT "File Name for Output ", FILE2$
     PRINT : PRINT "      PLOT CHOICE"
     PRINT "  1) No Plot Data"
     PRINT "  2) Create Data File for in-plane Shear Stress"
     PRINT "  3) Create Data File for Von Mises Stress"
     INPUT "     Choose 1 or 2 or 3"; IPL
     IF IPL < 1 OR IPL > 3 THEN IPL = 1   '--- default is no data
     IF IPL > 1 THEN INPUT "File Name for Plot Data "; FILE3$
     '----- Total dof is  NQ
     NQ = NDN * NN
     DIM X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     DIM TH(NE), DT(NE), NU(ND), U(ND), F(NQ)
     DIM D(3, 3), B(3, 8), DB(3, 8), SE(8, 8), Q(8), STR(3)
     DIM TL(8), XNI(4, 2), A(3, 4), G(4, 8), MPC(NMPC, 2), BT(NMPC, 3)
     GOSUB GETDATA
'----- Bandwidth NBW from Connectivity NOC()
     NBW = 0
     FOR I = 1 TO NE
        NMIN = NOC(I, 1): NMAX = NOC(I, 1)
        FOR J = 2 TO NEN
           IF NMIN > NOC(I, J) THEN NMIN = NOC(I, J)
           IF NMAX < NOC(I, J) THEN NMAX = NOC(I, J)
        NEXT J
        NTMP = NDN * (NMAX - NMIN + 1)
        IF NBW < NTMP THEN NBW = NTMP
     NEXT I
     FOR I = 1 TO NMPC
        NABS = ABS(MPC(I, 1) - MPC(I, 2)) + 1
        IF NBW < NABS THEN NBW = NABS
     NEXT I
     PRINT "The Bandwidth is"; NBW
     DIM S(NQ, NBW)
'---------- Stiffness Matrix ----------
     '----- Corner Nodes and Integration Points
     GOSUB INTEG
     FOR N = 1 TO NE
        PRINT "... Forming Stiffness Matrix of Element "; N
        GOSUB DMATX: '<--- D-Matrix
        GOSUB ELSTIF:   '<---Get Element Stiffness
        PRINT ".... Placing in Banded Locations"
        FOR II = 1 TO NEN
           NRT = NDN * (NOC(N, II) - 1)
           FOR IT = 1 TO NDN
              NR = NRT + IT
              I = NDN * (II - 1) + IT
              FOR JJ = 1 TO NEN
                 NCT = NDN * (NOC(N, JJ) - 1)
                 FOR JT = 1 TO NDN
                    J = NDN * (JJ - 1) + JT
                    NC = NCT + JT - NR + 1
                    IF NC > 0 THEN
                       S(NR, NC) = S(NR, NC) + SE(I, J)
                    END IF
                 NEXT JT
              NEXT JJ
              F(NR) = F(NR) + TL(I)
           NEXT IT
        NEXT II
     NEXT N
'----- Decide Penalty Parameter CNST -----
     CNST = 0
     FOR I = 1 TO NQ
        IF CNST < S(I, 1) THEN CNST = S(I, 1)
     NEXT I
     CNST = CNST * 10000
'----- Modify for Boundary Conditions -----
        '--- Displacement BC ---
     FOR I = 1 TO ND
        N = NU(I)
        S(N, 1) = S(N, 1) + CNST
        F(N) = F(N) + CNST * U(I)
     NEXT I
        '--- Multi-point Constraints ---
     FOR I = 1 TO NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1)
        S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2)
        IR = I1: IF IR > I2 THEN IR = I2
        IC = ABS(I2 - I1) + 1
        S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
     NEXT I
'----- Equation Solving
     PRINT ".... Solving Equations"
     GOSUB BANSOL
     OPEN FILE2$ FOR OUTPUT AS #2
     PRINT #2, "Output for Input Data in File --- "; FILE1$
     PRINT #2, TITLE$
     PRINT TITLE$
     IF LC = 1 THEN PRINT #2, "Plane Stress Analysis"
     IF LC = 2 THEN PRINT #2, "Plane Strain Analysis"
     PRINT #2, "NODE#   X-Displ     Y-Displ"
     PRINT "NODE#   X-Displ     Y-Displ"
     FOR I = 1 TO NN
        PRINT USING " ###"; I;
        PRINT USING "  ##.###^^^^"; F(2 * I - 1); F(2 * I)
        PRINT #2, USING " ###"; I;
        PRINT #2, USING "  ##.###^^^^"; F(2 * I - 1); F(2 * I)
     NEXT I
'----- Reaction Calculation
     PRINT "DOF#        Reaction"
     PRINT #2, "DOF#        Reaction"
     FFF1$ = "  ###     ##.####^^^^"
     FOR I = 1 TO ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
        PRINT USING FFF1$; N; R
        PRINT #2, USING FFF1$; N; R
     NEXT I
     IF IPL > 1 THEN
        OPEN FILE3$ FOR OUTPUT AS #3
        IF IPL = 2 THEN PRINT #3, "Max. in-plane Shear Stress ";
        IF IPL = 3 THEN PRINT #3, "Von Mises Stress ";
        PRINT #3, "(Element) for Data in File "; FILE1$
     END IF
'-----  Stress Calculations -----
     '--- Stresses at Integration Points
     PRINT #2, "ELEM#    Von Mises Stresses at 4 Integration points"
     FOR N = 1 TO NE
        PRINT #2, USING "###  "; N;
        FOR IP = 1 TO 4
           XI = XNI(IP, 1): ETA = XNI(IP, 2)
           GOSUB DMATX: '--- D Matrix
           GOSUB DBMAT: '--- Get DB Matrix
           '--- Stress Evaluation
           FOR I = 1 TO NEN
              IN = NDN * (NOC(N, I) - 1)
              II = NDN * (I - 1)
              FOR J = 1 TO NDN
                 Q(II + J) = F(IN + J)
              NEXT J
           NEXT I
           C1 = AL * DT(N): IF LC = 2 THEN C1 = C1 * (1 + PNU)
           FOR I = 1 TO 3
              C = 0
              FOR K = 1 TO 8
                 C = C + DB(I, K) * Q(K)
              NEXT K
              STR(I) = C - C1 * (D(I, 1) + D(I, 2))
           NEXT I
        '--- Von Mises Stress at Integration Point
           C = 0: IF LC = 2 THEN C = PNU * (STR(1) + STR(2))
           C1 = (STR(1) - STR(2)) ^ 2 + (STR(2) - C) ^ 2 + (C - STR(1)) ^ 2
           SV = SQR(.5 * C1 + 3 * STR(3) ^ 2)
           PRINT #2, USING " ##.####^^^^"; SV;
        '--- Maximum Shear Stress R
           R = SQR(.25 * (STR(1) - STR(2)) ^ 2 + (STR(3)) ^ 2)
           IF IPL = 2 THEN PRINT #3, R;
           IF IPL = 3 THEN PRINT #3, SV;
        NEXT IP
        PRINT #2,
        IF IPL > 1 THEN PRINT #3,
     NEXT N
        CLOSE #2: PRINT
     PRINT "-----    All Calculations are done    -----"
     PRINT "The Results are available in the text file "; FILE2$
     PRINT "View using a text processor"
     IF IPL > 1 THEN
        CLOSE #3: PRINT "Element Stress Data in file "; FILE3$
        PRINT "Run BESTFITQ and then CONTOURA or CONTOURB to plot stresses"
     END IF
     END
GETDATA:
'-----------  READ DATA  ------------
     '----- Coordinates
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N
        FOR J = 1 TO NDIM
           INPUT #1, X(N, J)
        NEXT J
     NEXT I
     '----- Connectivity, Material, Thickness, Temp-change
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), TH(N), DT(N)
     NEXT I
     '----- Displacement BC
     LINE INPUT #1, D$
     FOR I = 1 TO ND: INPUT #1, NU(I), U(I): NEXT I
     '----- Component Loads
     LINE INPUT #1, D$
     FOR I = 1 TO NL: INPUT #1, N, F(N): NEXT I
     '----- Material Properties
     LINE INPUT #1, D$
     FOR I = 1 TO NM
        INPUT #1, N
        FOR J = 1 TO NPR
           INPUT #1, PM(N, J)
        NEXT J
     NEXT I
     IF NMPC > 0 THEN
        '-----  Multi-point Constraints
        LINE INPUT #1, D$
        FOR I = 1 TO NMPC
           INPUT #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        NEXT I
     END IF
     CLOSE #1
     RETURN
INTEG:
'------- Integration Points XNI() --------
     C = .57735026919#
     XNI(1, 1) = -C: XNI(1, 2) = -C
     XNI(2, 1) = C: XNI(2, 2) = -C
     XNI(3, 1) = C: XNI(3, 2) = C
     XNI(4, 1) = -C: XNI(4, 2) = C
     RETURN
DMATX:
'-----  D() Matrix and Element Nodal Coordinates -----
     '--- Material Properties
     MATN = MAT(N)
     E = PM(MATN, 1): PNU = PM(MATN, 2)
     AL = PM(MATN, 3)
     '--- D() Matrix
     IF LC = 1 THEN
        '--- Plane Stress
        C1 = E / (1 - PNU ^ 2): C2 = C1 * PNU
     ELSE
        '--- Plane Strain
        C = E / ((1 + PNU) * (1 - 2 * PNU))
        C1 = C * (1 - PNU): C2 = C * PNU
     END IF
     C3 = .5 * E / (1 + PNU)
     D(1, 1) = C1: D(1, 2) = C2: D(1, 3) = 0
     D(2, 1) = C2: D(2, 2) = C1: D(2, 3) = 0
     D(3, 1) = 0: D(3, 2) = 0: D(3, 3) = C3
     RETURN
ELSTIF:
'--------  Element Stiffness and Temperature Load  -----
     FOR I = 1 TO 8: FOR J = 1 TO 8: SE(I, J) = 0: NEXT J: TL(I) = 0: NEXT I
     DTE = DT(N)
     '--- Weight Factor is ONE
     '--- Loop on Integration Points
     FOR IP = 1 TO 4
        '---  Get DB Matrix at Integration Point IP
        XI = XNI(IP, 1): ETA = XNI(IP, 2)
        GOSUB DBMAT
        '--- Element Stiffness Matrix  SE
        FOR I = 1 TO 8
           FOR J = 1 TO 8
              C = 0
              FOR K = 1 TO 3
                 C = C + B(K, I) * DB(K, J) * DJ * THICK
              NEXT K
              SE(I, J) = SE(I, J) + C
           NEXT J
        NEXT I
        '--- Determine Temperature Load TL
        C = AL * DTE: IF LC = 2 THEN C = (1 + PNU) * C
        FOR I = 1 TO 8
           TL(I) = TL(I) + THICK * DJ * C * (DB(1, I) + DB(2, I))
        NEXT I
     NEXT IP
     RETURN
DBMAT:
'-------  DB()  MATRIX  ------
     '--- Nodal Coordinates
     THICK = TH(N)
     N1 = NOC(N, 1): N2 = NOC(N, 2)
     N3 = NOC(N, 3): N4 = NOC(N, 4)
     X1 = X(N1, 1): Y1 = X(N1, 2)
     X2 = X(N2, 1): Y2 = X(N2, 2)
     X3 = X(N3, 1): Y3 = X(N3, 2)
     X4 = X(N4, 1): Y4 = X(N4, 2)
     '--- Formation of Jacobian  TJ
     TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4
     TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4
     TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4
     TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4
     '--- Determinant of the JACOBIAN
     DJ = TJ11 * TJ22 - TJ12 * TJ21
     '--- A(3,4) Matrix relates Strains to
     '--- Local Derivatives of u
     A(1, 1) = TJ22 / DJ: A(2, 1) = 0: A(3, 1) = -TJ21 / DJ
     A(1, 2) = -TJ12 / DJ: A(2, 2) = 0: A(3, 2) = TJ11 / DJ
     A(1, 3) = 0: A(2, 3) = -TJ21 / DJ: A(3, 3) = TJ22 / DJ
     A(1, 4) = 0: A(2, 4) = TJ11 / DJ: A(3, 4) = -TJ12 / DJ
     '--- G(4,8) Matrix relates Local Derivatives of u
     '--- to Local Nodal Displacements q(8)
     FOR I = 1 TO 4: FOR J = 1 TO 8
     G(I, J) = 0: NEXT J: NEXT I
     G(1, 1) = -(1 - ETA) / 4: G(2, 1) = -(1 - XI) / 4
     G(3, 2) = -(1 - ETA) / 4: G(4, 2) = -(1 - XI) / 4
     G(1, 3) = (1 - ETA) / 4: G(2, 3) = -(1 + XI) / 4
     G(3, 4) = (1 - ETA) / 4: G(4, 4) = -(1 + XI) / 4
     G(1, 5) = (1 + ETA) / 4: G(2, 5) = (1 + XI) / 4
     G(3, 6) = (1 + ETA) / 4: G(4, 6) = (1 + XI) / 4
     G(1, 7) = -(1 + ETA) / 4: G(2, 7) = (1 - XI) / 4
     G(3, 8) = -(1 + ETA) / 4: G(4, 8) = (1 - XI) / 4
     '--- B(3,8) Matrix Relates Strains to q
     FOR I = 1 TO 3
        FOR J = 1 TO 8
           C = 0
           FOR K = 1 TO 4
              C = C + A(I, K) * G(K, J)
           NEXT K
           B(I, J) = C
        NEXT J
     NEXT I
     '--- DB(3,8) Matrix relates Stresses to q(8)
     FOR I = 1 TO 3
        FOR J = 1 TO 8
           C = 0
           FOR K = 1 TO 3
              C = C + D(I, K) * B(K, J)
           NEXT K:
           DB(I, J) = C
        NEXT J
     NEXT I
     RETURN
BANSOL:
'-----  Band Solver  -----
     N1 = NQ - 1
     PRINT "--- Forward Elimination in Progress ---"
     FOR K = 1 TO N1
        NK = NQ - K + 1
        IF NK > NBW THEN NK = NBW
        FOR I = 2 TO NK
           C1 = S(K, I) / S(K, 1)
           I1 = K + I - 1
           FOR J = I TO NK
              J1 = J - I + 1
              S(I1, J1) = S(I1, J1) - C1 * S(K, J)
           NEXT J
           F(I1) = F(I1) - C1 * F(K)
        NEXT I
     NEXT K
     F(NQ) = F(NQ) / S(NQ, 1)
     PRINT "--- Back Substitution in Progress ---"
     FOR KK = 1 TO N1
        K = NQ - KK
        C1 = 1 / S(K, 1)
        F(K) = C1 * F(K)
        NK = NQ - K + 1
        IF NK > NBW THEN NK = NBW
        FOR J = 2 TO NK
           F(K) = F(K) - C1 * S(K, J) * F(K + J - 1)
        NEXT J
     NEXT KK
     RETURN

