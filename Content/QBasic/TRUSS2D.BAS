     '***************************************
     '*          PROGRAM TRUSS2D            *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '*            (C) 2001                 *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "2-D TRUSS ANALYSIS"; SPACE$(33);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Input Data File Name <DOS file name> ", FILE1$
     INPUT "Output Data File Name ", FILE2$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 2   'Material properties E, Alpha
     NQ = 2 * NN
     DIM X(NN, 2), NOC(NE, NEN), F(NQ), AREA(NE), MAT(NE)
     DIM PM(NM, 2), NU(ND), U(ND), DT(NE), SE(4, 4), TL(4)
     DIM MPC(NMPC, NPR), BT(NMPC, 3)
     '----- Coordinates -----
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N, X(N, 1), X(N, 2)
     NEXT I
     '----- Connectivity -----
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N, NOC(N, 1), NOC(N, 2)
        INPUT #1, MAT(N), AREA(N)
        INPUT #1, DT(N)
     NEXT I
     '----- Specified Displacements -----
     LINE INPUT #1, D$
     FOR I = 1 TO ND
        INPUT #1, NU(I), U(I)
     NEXT I
     '----- Component Loads -----
     LINE INPUT #1, D$
     FOR I = 1 TO NL
        INPUT #1, N, F(N)
     NEXT I
     '----- Material Properties -----
     LINE INPUT #1, D$
     FOR I = 1 TO NM
        INPUT #1, N
        FOR J = 1 TO NPR
           INPUT #1, PM(N, J)
        NEXT J
     NEXT I
     '----- Multi-point Constraints B1*Qi+B2*Qj=B0
     IF NMPC > 0 THEN
        LINE INPUT #1, A$
        FOR I = 1 TO NMPC
           INPUT #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        NEXT I
     END IF
     CLOSE #1
     '----- Bandwidth Evaluation -----
     NBW = 0
     FOR N = 1 TO NE
        NABS = NDN * (ABS(NOC(N, 1) - NOC(N, 2)) + 1)
        IF NBW < NABS THEN NBW = NABS
     NEXT N
     FOR I = 1 TO NMPC
        NABS = ABS(MPC(I, 1) - MPC(I, 2)) + 1
        IF NBW < NABS THEN NBW = NABS
     NEXT I
     PRINT "Bandwidth = "; NBW
     DIM S(NQ, NBW)
     '-----  Stiffness Matrix -----
     FOR N = 1 TO NE
        I1 = NOC(N, 1): I2 = NOC(N, 2)
        I3 = MAT(N)
        X21 = X(I2, 1) - X(I1, 1)
        Y21 = X(I2, 2) - X(I1, 2)
        EL = SQR(X21 * X21 + Y21 * Y21)
        EAL = PM(I3, 1) * AREA(N) / EL
        CS = X21 / EL: SN = Y21 / EL
     '----------- Element Stiffness Matrix SE() -----------
        SE(1, 1) = CS * CS * EAL
        SE(1, 2) = CS * SN * EAL: SE(2, 1) = SE(1, 2)
        SE(1, 3) = -CS * CS * EAL: SE(3, 1) = SE(1, 3)
        SE(1, 4) = -CS * SN * EAL: SE(4, 1) = SE(1, 4)
        SE(2, 2) = SN * SN * EAL
        SE(2, 3) = -CS * SN * EAL: SE(3, 2) = SE(2, 3)
        SE(2, 4) = -SN * SN * EAL: SE(4, 2) = SE(2, 4)
        SE(3, 3) = CS * CS * EAL
        SE(3, 4) = CS * SN * EAL: SE(4, 3) = SE(3, 4)
        SE(4, 4) = SN * SN * EAL
     '-------------- Temperature Load TL() ---------------
        EE0 = PM(I3, 2) * DT(N) * PM(I3, 1) * AREA(N)
        TL(1) = -EE0 * CS: TL(2) = -EE0 * SN
        TL(3) = EE0 * CS: TL(4) = EE0 * SN
        PRINT "..... Adding Stiffness to Global Locations"
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
     '----- Equation Solving using Band Solver -----
     GOSUB BANSOL
     OPEN FILE2$ FOR OUTPUT AS #2
     PRINT #2, "Output for Input Data from File "; FILE1$
     PRINT TITLE$
     PRINT #2, TITLE$
     PRINT "Node#    X-Displ       Y-Displ"
     PRINT #2, "Node#    X-Displ       Y-Displ"
     FOR I = 1 TO NN
        PRINT USING " ###"; I;
        PRINT #2, USING " ###"; I;
        PRINT USING "   ##.####^^^^"; F(2 * I - 1); F(2 * I)
        PRINT #2, USING "   ##.####^^^^"; F(2 * I - 1); F(2 * I)
     NEXT I
     '----- Stress Calculation -----
     PRINT "Elem#    Stress"
     PRINT #2, "Elem#    Stress"
     FOR I = 1 TO NE
        I1 = NOC(I, 1): I2 = NOC(I, 2)
        I3 = MAT(I)
        X21 = X(I2, 1) - X(I1, 1): Y21 = X(I2, 2) - X(I1, 2)
        EL = SQR(X21 * X21 + Y21 * Y21)
        CS = X21 / EL: SN = Y21 / EL
        J2 = 2 * I1: J1 = J2 - 1: K2 = 2 * I2: K1 = K2 - 1
        DLT = (F(K1) - F(J1)) * CS + (F(K2) - F(J2)) * SN
        STRESS = PM(I3, 1) * (DLT / EL - PM(I3, 2) * DT(I))
        PRINT USING " ###"; I;
        PRINT #2, USING " ###"; I;
        PRINT USING "   ##.####^^^^"; STRESS
        PRINT #2, USING "   ##.###^^^^"; STRESS
     NEXT I
     '----- Reaction Calculation -----
     PRINT "DOF#      Reaction"
     PRINT #2, "DOF#      Reaction"
     FFF1$ = " ###   ##.####^^^^"
     FOR I = 1 TO ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
        PRINT USING FFF1$; N; R
        PRINT #2, USING FFF1$; N; R
     NEXT I
     CLOSE #2
     PRINT "Results are now available in the text file "; FILE2$
     END
BANSOL:
     N1 = NQ - 1
     '----- Forward Elimination -----
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
     '----- Back Substitution -----
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

