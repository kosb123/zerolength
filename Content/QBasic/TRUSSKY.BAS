     '****************************************
     '*          PROGRAM TRUSSKY             *
     '*  T.R.Chandrupatla and A.D.Belegundu  *
     '****************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "2-D TRUSS ANALYSIS (SKYLINE)"; SPACE$(23);
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
     NPR = 2   'Material properties E, Alpha
     NQ = 2 * NN
     DIM X(NN, 2), NOC(NE, NEN), F(NQ), AREA(NE), MAT(NE)
     DIM PM(NM, NPR), NU(ND), U(ND), DT(NE), SE(4, 4), TL(4)
     DIM MPC(NMPC, 2), BT(NMPC, 3), ID(NQ)
     '----- Coordinates -----
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N, X(N, 1), X(N, 2)
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
        LINE INPUT #1, A$
        FOR I = 1 TO NMPC
           INPUT #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        NEXT I
     END IF
     CLOSE #1
     '----- Preparation for Skyline Storage -----
     FOR I = 1 TO NE
        II = NOC(I, 1): IF II > NOC(I, 2) THEN II = NOC(I, 2)
        II1 = NDN * (II - 1) + 1
        FOR J = 1 TO NEN
           NIJ = NDN * (NOC(I, J) - 1)
           FOR JJ = 1 TO NDN
              NDG = NIJ + JJ: NHT = NDG - II1 + 1
              IF NHT > ID(NDG) THEN ID(NDG) = NHT
           NEXT JJ
        NEXT J
     NEXT I
     '----- Skyline Height adjustment for MPC
     FOR I = 1 TO NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        NDG = I1: IF NDG < I2 THEN NDG = I2
        NHT = ABS(I2 - I1) + 1
        IF NHT > ID(NDG) THEN ID(NDG) = NHT
     NEXT I
     FOR I = 2 TO NQ
        ID(I) = ID(I) + ID(I - 1)
     NEXT I
     NSUM = ID(NQ)
     DIM A(NSUM)
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
           NCT = NDN * (NOC(N, II) - 1)
           FOR IT = 1 TO NDN
              NC = NCT + IT: IID = ID(NC)
              I = NDN * (II - 1) + IT
              FOR JJ = 1 TO NEN
                 NRT = NDN * (NOC(N, JJ) - 1)
                 FOR JT = 1 TO 2
                    J = NDN * (JJ - 1) + JT
                    NR = NRT + JT
                    IF NR <= NC THEN
                       NLOC = IID - (NC - NR)
                       A(NLOC) = A(NLOC) + SE(I, J)
                    END IF
                 NEXT JT
              NEXT JJ
              F(NC) = F(NC) + TL(I)
           NEXT IT
        NEXT II
     NEXT N
     '----- Decide Penalty Parameter CNST -----
     CNST = 0
     FOR I = 1 TO NQ
        II = ID(I)
        IF CNST < A(II) THEN CNST = A(II)
     NEXT I
     CNST = CNST * 10000
     '----- Modify for Boundary Conditions -----
        '--- Displacement BC ---
     FOR I = 1 TO ND
        N = NU(I): II = ID(N)
        A(II) = A(II) + CNST
        F(N) = F(N) + CNST * U(I)
     NEXT I
        '--- Multi-point Constraints ---
     FOR I = 1 TO NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        A(ID(I1)) = A(ID(I1)) + CNST * BT(I, 1) * BT(I, 1)
        A(ID(I2)) = A(ID(I2)) + CNST * BT(I, 2) * BT(I, 2)
        II = I1: IF II < I2 THEN II = I2
        IL = ID(II) - ABS(I2 - I1)
        A(IL) = A(IL) + CNST * BT(I, 1) * BT(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
     NEXT I
     '----- Equation Solving using Skyline Solver -----
     GOSUB SKYLINE
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
SKYLINE:
     '-----  Sky Line Solver -----
     '--- Forward Elimination ---
     FOR J = 2 TO NQ
        NJ = ID(J) - ID(J - 1)
        IF NJ <> 1 THEN
           K1 = 0: NJ = J - NJ + 1
           FOR K = NJ TO J - 1
              K1 = K1 + 1: KJ = ID(J - 1) + K1: KK = ID(K)
              C = A(KJ) / A(KK)
              FOR I = K + 1 TO J
                 NI = ID(I) - ID(I - 1)
                 IF (I - K + 1) <= NI THEN
                    IJ = ID(J) - J + I
                    KI = ID(I) - I + K
                    A(IJ) = A(IJ) - C * A(KI)
                 END IF
              NEXT I
           NEXT K
        END IF
     NEXT J
     FOR K = 1 TO NQ - 1: KK = ID(K)
        C = F(K) / A(KK)
        FOR I = K + 1 TO NQ
           NI = ID(I) - ID(I - 1)
           IF (I - K + 1) <= NI THEN
              KI = ID(I) - I + K
              F(I) = F(I) - C * A(KI)
           END IF
        NEXT I
     NEXT K
     '--- Back-substitution ---
     NS = ID(NQ): F(NQ) = F(NQ) / A(NS)
     FOR I1 = 1 TO NQ - 1
        I = NQ - I1
        II = ID(I)
        C = 1 / A(II): F(I) = C * F(I)
        FOR J = I + 1 TO NQ
           J1 = J - I + 1: NJ = ID(J) - ID(J - 1)
           IF J1 <= NJ THEN
              IJ = ID(J) - J + I
              F(I) = F(I) - C * A(IJ) * F(J)
           END IF
        NEXT J
     NEXT I1
     RETURN

