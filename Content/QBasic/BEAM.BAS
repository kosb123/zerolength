     '**************************************"
     '*             PROGRAM BEAM           *"
     '*        Beam Bending Analysis       *"
     '* T.R.Chandrupatla and A.D.Belegundu *"
     '**************************************"
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "BEAM BENDING ANALYSIS"; SPACE$(30);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "File Name for Input Data ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: LINE INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     INPUT "File Name for Output ", FILE2$
     '----- Total dof is  NQ
     NPR = 1   'Material property E
     NQ = 2 * NN
     DIM X(NN), NOC(NE, NEN), MAT(NE), PM(NM, NPR), SMI(NE)
     DIM NU(ND), U(ND), F(NQ), SE(4, 4), MPC(NMPC, 2), BT(NMPC, 3)
     GOSUB GETDATA
'----- Bandwidth NBW from Connectivity NOC()
     NBW = 0
     FOR I = 1 TO NE
        NTMP = NDN * (ABS(NOC(I, 1) - NOC(I, 2)) + 1)
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
         N1 = NOC(N, 1)
         N2 = NOC(N, 2)
         M = MAT(N)
         EL = ABS(X(N1) - X(N2))
         EIL = PM(M, 1) * SMI(N) / EL ^ 3
         SE(1, 1) = 12 * EIL
         SE(1, 2) = EIL * 6 * EL
         SE(1, 3) = -12 * EIL
         SE(1, 4) = EIL * 6 * EL
            SE(2, 1) = SE(1, 2)
            SE(2, 2) = EIL * 4 * EL * EL
            SE(2, 3) = -EIL * 6 * EL
            SE(2, 4) = EIL * 2 * EL * EL
         SE(3, 1) = SE(1, 3)
         SE(3, 2) = SE(2, 3)
         SE(3, 3) = EIL * 12
         SE(3, 4) = -EIL * 6 * EL
            SE(4, 1) = SE(1, 4)
            SE(4, 2) = SE(2, 4)
            SE(4, 3) = SE(3, 4)
            SE(4, 4) = EIL * 4 * EL * EL
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
     PRINT #2, "Output for Input Data in File --- ", FILE1$
     PRINT #2, TITLE$
     PRINT TITLE$
     PRINT "Node#   Displ.       Rotation"
     PRINT #2, "Node#   Displ.       Rotation"
     FOR I = 1 TO NN
     PRINT USING "  ###"; I;
     PRINT #2, USING "  ###"; I;
     PRINT USING "  ##.####^^^^"; F(2 * I - 1); F(2 * I)
     PRINT #2, USING "  ##.####^^^^"; F(2 * I - 1); F(2 * I)
     NEXT I
'----- Reaction Calculation
     PRINT "DOF#   Reaction"
     PRINT #2, "DOF#   Reaction"
     FOR I = 1 TO ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
        PRINT USING " ###"; N; : PRINT USING "  ##.####^^^^  "; R
        PRINT #2, USING " ###"; N; : PRINT #2, USING "  ##.####^^^^  "; R
     NEXT I
     END
GETDATA:
     '===============  READ DATA  ====================
     '----- Coordinates
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N
           INPUT #1, X(N)
     NEXT I
     '----- Connectivity, Material, Mom_Inertia, Dummy
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), SMI(N)
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

