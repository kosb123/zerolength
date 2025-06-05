     '***************************************
     '*           PROGRAM  CST              *
     '*   2-D  CONSTANT STRAIN TRIANGLE     *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "2-D CONSTANT STRAIN TRIANGLE"; SPACE$(23);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     PRINT "  1) Plane Stress Analysis"
     PRINT "  2) Plane Strain Analysis"
     INPUT "     Choose 1 or 2 "; LC
     IF LC < 1 OR LC > 2 THEN LC = 1   '--- default is Plane Stress
     PRINT : INPUT "File Name for Input Data ", FILE1$
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
     DIM TH(NE), DT(NE), NU(ND), U(ND), F(NQ), MPC(NMPC, 2), BT(NMPC, 3)
     DIM D(3, 3), B(3, 6), DB(3, 6), SE(6, 6), Q(6), STR(3), TL(6)
     GOSUB GETDATA
'----- Bandwidth NBW from Connectivity NOC() and MPC
     NBW = 0
     FOR I = 1 TO NE
        NMIN = NOC(I, 1): NMAX = NOC(I, 1)
        FOR J = 2 TO 3
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
'----- Global Stiffness Matrix
     FOR N = 1 TO NE
        PRINT "Forming Stiffness Matrix of Element "; N
        GOSUB DBMAT
     '--- Element Stiffness
        FOR I = 1 TO 6
           FOR J = 1 TO 6
              C = 0
              FOR K = 1 TO 3
                 C = C + .5 * ABS(DJ) * B(K, I) * DB(K, J) * TH(N)
              NEXT K
              SE(I, J) = C
           NEXT J
        NEXT I
     '--- Temperature Load Vector
        C = AL * DT(N): IF LC = 2 THEN C = C * (1 + PNU)
        FOR I = 1 TO 6
           TL(I) = .5 * C * TH(N) * ABS(DJ) * (DB(1, I) + DB(2, I))
        NEXT I
        PRINT ".... Placing in Global Locations"
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
'-----  Stress Calculations
     PRINT #2, "ELEM#    SX         SY         TXY";
     PRINT #2, "        S1         S2     ANGLE SX->S1"
     FOR N = 1 TO NE
        GOSUB DBMAT
        GOSUB STRESS
     '--- Principal Stress Calculations
        IF STR(3) = 0 THEN
           S1 = STR(1): S2 = STR(2): ANG = 0
           IF S2 > S1 THEN
              S1 = STR(2): S2 = STR(1): ANG = 90
           END IF
        ELSE
           C = .5 * (STR(1) + STR(2))
           R = SQR(.25 * (STR(1) - STR(2)) ^ 2 + (STR(3)) ^ 2)
           S1 = C + R: S2 = C - R
           IF C > STR(1) THEN
              ANG = 57.2957795# * ATN(STR(3) / (S1 - STR(1)))
              IF STR(3) > 0 THEN ANG = 90 - ANG
              IF STR(3) < 0 THEN ANG = -90 - ANG
           ELSE
              ANG = 57.29577951# * ATN(STR(3) / (STR(1) - S2))
           END IF
        END IF
        PRINT #2, USING " ###"; N;
        PRINT #2, USING " ##.###^^^^"; STR(1); STR(2); STR(3);
        PRINT #2, USING " ##.###^^^^"; S1; S2; ANG
        IF IPL = 2 THEN PRINT #3, .5 * (S1 - S2)
        IF IPL = 3 THEN
           S3 = 0: IF LC = 2 THEN S3 = PNU * (S1 + S2)
           C = (S1 - S2) ^ 2 + (S2 - S3) ^ 2 + (S3 - S1) ^ 2
           PRINT #3, SQR(.5 * C)
        END IF
     NEXT N
     CLOSE #2
     PRINT "Complete results are in file "; FILE2$
     IF IPL > 1 THEN
        CLOSE #3: PRINT "Element Stress Data in file "; FILE3$
        PRINT "Run BESTFIT and then CONTOURA or CONTOURB to plot stresses"
     END IF
     END
GETDATA:
     '===============  READ DATA  ====================
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
DBMAT:
     '----- D(), B() and DB() matrices
     '--- First the D-Matrix
     M = MAT(N): E = PM(M, 1): PNU = PM(M, 2): AL = PM(M, 3)
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
     '--- Strain-Displacement Matrix B()
     I1 = NOC(N, 1): I2 = NOC(N, 2): I3 = NOC(N, 3)
     X1 = X(I1, 1): Y1 = X(I1, 2)
     X2 = X(I2, 1): Y2 = X(I2, 2)
     X3 = X(I3, 1): Y3 = X(I3, 2)
     X21 = X2 - X1: X32 = X3 - X2: X13 = X1 - X3
     Y12 = Y1 - Y2: Y23 = Y2 - Y3: Y31 = Y3 - Y1
     DJ = X13 * Y23 - X32 * Y31   'DJ is determinant of Jacobian
     '--- Definition of B() Matrix
     B(1, 1) = Y23 / DJ: B(2, 1) = 0: B(3, 1) = X32 / DJ
     B(1, 2) = 0: B(2, 2) = X32 / DJ: B(3, 2) = Y23 / DJ
     B(1, 3) = Y31 / DJ: B(2, 3) = 0: B(3, 3) = X13 / DJ
     B(1, 4) = 0: B(2, 4) = X13 / DJ: B(3, 4) = Y31 / DJ
     B(1, 5) = Y12 / DJ: B(2, 5) = 0: B(3, 5) = X21 / DJ
     B(1, 6) = 0: B(2, 6) = X21 / DJ: B(3, 6) = Y12 / DJ
     '--- DB Matrix DB = D*B
     FOR I = 1 TO 3
        FOR J = 1 TO 6
           C = 0
           FOR K = 1 TO 3
              C = C + D(I, K) * B(K, J)
           NEXT K
           DB(I, J) = C
        NEXT J
     NEXT I
     RETURN
STRESS:
     '----- Stress Evaluation
     Q(1) = F(2 * I1 - 1): Q(2) = F(2 * I1)
     Q(3) = F(2 * I2 - 1): Q(4) = F(2 * I2)
     Q(5) = F(2 * I3 - 1): Q(6) = F(2 * I3)
     C1 = AL * DT(N): IF LC = 2 THEN C1 = C1 * (1 + PNU)
     FOR I = 1 TO 3
        C = 0
        FOR K = 1 TO 6
           C = C + DB(I, K) * Q(K)
        NEXT K
        STR(I) = C - C1 * (D(I, 1) + D(I, 2))
     NEXT I
     RETURN
BANSOL:
     '----- Band Solver -----
     N1 = NQ - 1
     '--- Forward Elimination
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
     '--- Back-substitution
     F(NQ) = F(NQ) / S(NQ, 1)
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

