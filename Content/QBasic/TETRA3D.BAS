     '*******       PROGRAM TETRA3D       *******
     '*    Three Dimensional Stress Analysis    *
     '*          Tetrahedral Elements           *
     '*   T.R.Chandrupatla and A.D.Belegundu    *
     '*******************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "3-D TETRAHEDRAL ELEMENT"; SPACE$(28);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25
     PRINT : INPUT "File Name for Input Data ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 3   'Material properties E, Nu, Alpha
     INPUT "File Name for Output ", FILE2$
     '----- Total dof is  NQ
     NQ = NDN * NN
     DIM X(NN, NDIM), F(NQ), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     DIM NU(ND), U(ND), MPC(NMPC, 2), BT(NMPC, 3), QT(12)
     DIM D(6, 6), B(6, 12), DB(6, 12), SE(12, 12), STR(6), DT(NE)
     GOSUB GETDATA
'----- Bandwidth NBW from Connectivity NOC() and MPC
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
'----- Global Stiffness Matrix
     FOR N = 1 TO NE
        PRINT "Forming Stiffness Matrix of Element "; N
        GOSUB DBMAT
     '--- Element Stiffness
     FOR I = 1 TO 12
        FOR J = 1 TO 12
        SE(I, J) = 0
        FOR K = 1 TO 6
           SE(I, J) = SE(I, J) + B(K, I) * DB(K, J) * ABS(DJ) / 6
        NEXT K: NEXT J: NEXT I
     '--- Temperature Load Vector QT()
        C = AL * DT(N)
        FOR I = 1 TO 12
           DSUM = DB(1, I) + DB(2, I) + DB(3, I)
           QT(I) = C * ABS(DJ) * DSUM / 6
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
              F(NR) = F(NR) + QT(I)
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
     PRINT TITLE$
     PRINT #2, TITLE$
     PRINT " NODE#  X-Displ    Y-Displ    Z-Displ"
     PRINT #2, " NODE#  X-Displ    Y-Displ    Z-Displ"
     FOR I = 1 TO NN
     PRINT USING "#### "; I;
     PRINT USING " #.####^^^^"; F(3 * I - 2); F(3 * I - 1); F(3 * I)
     PRINT #2, USING " ####"; I;
     PRINT #2, USING " #.####^^^^"; F(3 * I - 2); F(3 * I - 1); F(3 * I)
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
'--- Stress Calculations
     PI = 3.141593
     FOR N = 1 TO NE
        GOSUB DBMAT
        GOSUB STRESS
        '--- Principal Stress Calculations
        AI1 = STR(1) + STR(2) + STR(3)
        AI21 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1)
        AI22 = STR(4) * STR(4) + STR(5) * STR(5) + STR(6) * STR(6)
        AI2 = AI21 - AI22
        AI31 = STR(1) * STR(2) * STR(3) + 2 * STR(4) * STR(5) * STR(6)
        AI32 = STR(1) * STR(4) ^ 2 + STR(2) * STR(5) ^ 2 + STR(3) * STR(6) ^ 2
        AI3 = AI31 - AI32
        C1 = AI2 - AI1 ^ 2 / 3
        C2 = -2 * (AI1 / 3) ^ 3 + AI1 * AI2 / 3 - AI3
        C3 = 2 * SQR(-C1 / 3)
        TH = -3 * C2 / (C1 * C3): TH2 = ABS(1 - TH * TH)
        IF TH = 0 THEN TH = PI / 2
        IF TH > 0 THEN TH = ATN(SQR(TH2) / TH)
        IF TH < 0 THEN TH = PI - ATN(SQR(TH2) / TH)
        TH = TH / 3
        '--- Principal Stresses
        P1 = AI1 / 3 + C3 * COS(TH)
        P2 = AI1 / 3 + C3 * COS(TH + 2 * PI / 3)
        P3 = AI1 / 3 + C3 * COS(TH + 4 * PI / 3)
        PRINT #2,
        PRINT #2, "STRESSES IN ELEMENT NO.  "; N
        PRINT #2, "  Normal Stresses SX,SY,SZ"
        PRINT #2, USING "  #.####^^^^"; STR(1); STR(2); STR(3)
        PRINT #2, "  Shear Stresses TYZ,TXZ,TXY"
        PRINT #2, USING "  #.####^^^^"; STR(4); STR(5); STR(6)
        PRINT #2,
        PRINT #2, "  Principal Stresses"
        PRINT #2, USING "  #.####^^^^"; P1; P2; P3
     NEXT N
     CLOSE #2
     PRINT "Results are in the file  "; FILE2$
     END
DBMAT:
     '----- D(), B() and DB() matrices
     '--- First the D-Matrix
     M = MAT(N): E = PM(M, 1): PNU = PM(M, 2): AL = PM(M, 3)
     C4 = E / ((1 + PNU) * (1 - 2 * PNU))
     C1 = C4 * (1 - PNU): C2 = C4 * PNU: C3 = .5 * E / (1 + PNU)
     FOR I = 1 TO 6: FOR J = 1 TO 6: D(I, J) = 0: NEXT J: NEXT I
     D(1, 1) = C1: D(1, 2) = C2: D(1, 3) = C2
     D(2, 1) = C2: D(2, 2) = C1: D(2, 3) = C2
     D(3, 1) = C2: D(3, 2) = C2: D(3, 3) = C1
     D(4, 4) = C3: D(5, 5) = C3: D(6, 6) = C3
     '--- Strain-Displacement Matrix B()
     I1 = NOC(N, 1): I2 = NOC(N, 2): I3 = NOC(N, 3): I4 = NOC(N, 4)
     X14 = X(I1, 1) - X(I4, 1): X24 = X(I2, 1) - X(I4, 1): X34 = X(I3, 1) - X(I4, 1)
     Y14 = X(I1, 2) - X(I4, 2): Y24 = X(I2, 2) - X(I4, 2): Y34 = X(I3, 2) - X(I4, 2)
     Z14 = X(I1, 3) - X(I4, 3): Z24 = X(I2, 3) - X(I4, 3): Z34 = X(I3, 3) - X(I4, 3)
     DJ1 = X14 * (Y24 * Z34 - Z24 * Y34)
     DJ2 = Y14 * (Z24 * X34 - X24 * Z34)
     DJ3 = Z14 * (X24 * Y34 - Y24 * X34)
     DJ = DJ1 + DJ2 + DJ3
     A11 = (Y24 * Z34 - Z24 * Y34) / DJ
     A21 = (Z24 * X34 - X24 * Z34) / DJ
     A31 = (X24 * Y34 - Y24 * X34) / DJ
     A12 = (Y34 * Z14 - Z34 * Y14) / DJ
     A22 = (Z34 * X14 - X34 * Z14) / DJ
     A32 = (X34 * Y14 - Y34 * X14) / DJ
     A13 = (Y14 * Z24 - Z14 * Y24) / DJ
     A23 = (Z14 * X24 - X14 * Z24) / DJ
     A33 = (X14 * Y24 - Y14 * X24) / DJ
     '---  B Matrix
     FOR I = 1 TO 6: FOR J = 1 TO 12: B(I, J) = 0: NEXT J: NEXT I
     B(1, 1) = A11: B(1, 4) = A12: B(1, 7) = A13: B(1, 10) = -A11 - A12 - A13
     B(2, 2) = A21: B(2, 5) = A22: B(2, 8) = A23: B(2, 11) = -A21 - A22 - A23
     B(3, 3) = A31: B(3, 6) = A32: B(3, 9) = A33: B(3, 12) = -A31 - A32 - A33
     B(4, 2) = A31: B(4, 3) = A21: B(4, 5) = A32: B(4, 6) = A22: B(4, 8) = A33
                B(4, 9) = A23: B(4, 11) = B(3, 12): B(4, 12) = B(2, 11)
     B(5, 1) = A31: B(5, 3) = A11: B(5, 4) = A32: B(5, 6) = A12: B(5, 7) = A33
                B(5, 9) = A13: B(5, 10) = B(3, 12): B(5, 12) = B(1, 10)
     B(6, 1) = A21: B(6, 2) = A11: B(6, 4) = A22: B(6, 5) = A12: B(6, 7) = A23
                B(6, 8) = A13: B(6, 10) = B(2, 11): B(6, 11) = B(1, 10)
     '--- DB Matrix DB = D*B
     FOR I = 1 TO 6
     FOR J = 1 TO 12
     DB(I, J) = 0
     FOR K = 1 TO 6
     DB(I, J) = DB(I, J) + D(I, K) * B(K, J)
     NEXT K: NEXT J: NEXT I
     RETURN
STRESS:
     'Stress Evaluation (Element Nodal Displacements stored in QT() )
     FOR I = 1 TO 4
        IN = 3 * (NOC(N, I) - 1): II = 3 * (I - 1)
        FOR J = 1 TO 3
           QT(II + J) = F(IN + J): NEXT J: NEXT I
           C1 = AL * DT(N)
           FOR I = 1 TO 6
              STR(I) = 0
              FOR K = 1 TO 12
                 STR(I) = STR(I) + DB(I, K) * QT(K)
              NEXT K
              STR(I) = STR(I) - C1 * (D(I, 1) + D(I, 2) + D(I, 3))
           NEXT I
     RETURN
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
     '----- Connectivity, Material, Temp-change
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), DT(N)
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

