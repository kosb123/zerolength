     '***************************************
     '*           PROGRAM FEM1D             *
     '*    WITH MULTI-POINT CONSTRAINTS     *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "1-D BAR PROBLEMS"; SPACE$(35);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Input Data File Name <DOS file name> ", FILE1$
     INPUT "Output Data File Name ", FILE2$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$
     LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 2   'Material Properties E, Alpha
     DIM X(NN), NOC(NE, NEN), F(NN), AREA(NE), MAT(NE), DT(NE)
     DIM PM(NM, NPR), NU(ND), U(ND), MPC(NMPC, 2), BT(NMPC, 3)
     '----- Coordinates -----
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N, X(N)
     NEXT I
     '----- Connectivity -----
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N, NOC(N, 1), NOC(N, 2)
        INPUT #1, MAT(N), AREA(N), DT(N)
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
        LINE INPUT #1, D$
        FOR I = 1 TO NMPC
           INPUT #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        NEXT I
     ENDIF
     CLOSE #1
     '----- Bandwidth Evaluation -----
     NBW = 0
     FOR N = 1 TO NE
        NABS = ABS(NOC(N, 1) - NOC(N, 2)) + 1
        IF NBW < NABS THEN NBW = NABS
     NEXT N
     FOR I = 1 TO NMPC
        NABS = ABS(MPC(I, 1) - MPC(I, 2)) + 1
        IF NBW < NABS THEN NBW = NABS
     NEXT I
     DIM S(NN, NBW)
     '----- Stiffness Matrix -----
     FOR N = 1 TO NE
        N1 = NOC(N, 1): N2 = NOC(N, 2): N3 = MAT(N)
        X21 = X(N2) - X(N1): EL = ABS(X21)
        EAL = PM(N3, 1) * AREA(N) / EL
        IF NPR > 1 THEN C = PM(N3, 2)
        TL = PM(N3, 1) * C * DT(N) * AREA(N) * EL / X21
        '----- Temperature Loads -----
        F(N1) = F(N1) - TL
        F(N2) = F(N2) + TL
        '----- Element Stiffness in Global Locations -----
        S(N1, 1) = S(N1, 1) + EAL
        S(N2, 1) = S(N2, 1) + EAL
        IR = N1: IF IR > N2 THEN IR = N2
        IC = ABS(N2 - N1) + 1
        S(IR, IC) = S(IR, IC) - EAL
     NEXT N
     '----- Decide Penalty Parameter CNST -----
     CNST = 0
     FOR I = 1 TO NN
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
     PRINT TITLE$
     PRINT #2, TITLE$
     PRINT "NODE NO.", "DISPLACEMENT"
     PRINT #2, "NODE NO.", "DISPLACEMENT"
     FOR I = 1 TO NN
        PRINT I, F(I)
        PRINT #2, I, F(I)
     NEXT I
     '----- Stress Calculation -----
     PRINT "ELEM NO.", "STRESS"
     PRINT #2, "ELEM NO.", "STRESS"
     FOR N = 1 TO NE
        N1 = NOC(N, 1): N2 = NOC(N, 2): N3 = MAT(N)
        EPS = (F(N2) - F(N1)) / (X(N2) - X(N1))
        IF NPR > 1 THEN C = PM(N3, 2)
        STRESS = PM(N3, 1) * (EPS - C * DT(N))
        PRINT N, STRESS
        PRINT #2, N, STRESS
     NEXT N
     '----- Reaction Calculation -----
     PRINT "NODE NO.", "REACTION"
     PRINT #2, "NODE NO.", "REACTION"
     FOR I = 1 TO ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
        PRINT N, R
        PRINT #2, N, R
     NEXT I
     CLOSE #2
     PRINT "RESULTS ARE IN FILE "; FILE2$
     END
BANSOL:
     N = NN
     '----- Forward Elimination -----
     FOR K = 1 TO N - 1
        NBK = N - K + 1
        IF N - K + 1 > NBW THEN NBK = NBW
        FOR I = K + 1 TO NBK + K - 1
           I1 = I - K + 1
           C = S(K, I1) / S(K, 1)
           FOR J = I TO NBK + K - 1
              J1 = J - I + 1
              J2 = J - K + 1
              S(I, J1) = S(I, J1) - C * S(K, J2)
           NEXT J
           F(I) = F(I) - C * F(K)
        NEXT I
     NEXT K
     '----- Back Substitution -----
     F(N) = F(N) / S(N, 1)
     FOR II = 1 TO N - 1
        I = N - II
        NBI = N - I + 1
        IF N - I + 1 > NBW THEN NBI = NBW
        SUM = 0!
        FOR J = 2 TO NBI
           SUM = SUM + S(I, J) * F(I + J - 1)
        NEXT J
        F(I) = (F(I) - SUM) / S(I, 1)
     NEXT II
     RETURN

