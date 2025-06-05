     '********    PROGRAM FRAME2D      ********
     '*      2-D   FRAME ANALYSIS BY FEM      *
     '*  T.R.Chandrupatla and A.D.Belegundu   *
     '*****************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "2-D FRAME ANALYSIS"; SPACE$(33);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "File Name for Input Data ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 1   'Material property E
     INPUT "File Name for Output ", FILE2$
     '----- Total dof is  NQ
     NQ = NDN * NN
     DIM X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR), ARIN(NE, 2)
     DIM NU(ND), U(ND), F(NQ), SEP(6, 6), SE(6, 6), MPC(NMPC, 2), UDL(NE)
     DIM BT(NMPC, 3), DCOS(3, 3), ALAMBDA(6, 6), ED(6), EDP(6), EF(6)
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
        ISTF = 2
        GOSUB ELSTIF
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
'----- Loads due to uniformly distributed load on element
     FOR N = 1 TO NE
     IF ABS(UDL(N)) > 0 THEN
       ISTF = 1
       GOSUB ELSTIF
       I1 = NOC(N, 1): I2 = NOC(N, 2)
       ED(1) = 0: ED(4) = 0
       ED(2) = UDL(N) * EL / 2: ED(5) = ED(2)
       ED(3) = UDL(N) * EL ^ 2 / 12: ED(6) = -ED(3)
       FOR I = 1 TO 6
         EDP(I) = 0
         FOR K = 1 TO 6
           EDP(I) = EDP(I) + ALAMBDA(K, I) * ED(K)
         NEXT K
       NEXT I
       FOR I = 1 TO 3
         F(3 * I1 - 3 + I) = F(3 * I1 - 3 + I) + EDP(I)
         F(3 * I2 - 3 + I) = F(3 * I2 - 3 + I) + EDP(I + 3)
       NEXT I
     END IF
     NEXT N
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
     PRINT TITLE$
     PRINT #2, TITLE$
     PRINT "NODE#   X-Displ.    Y-Displ.    Z-Rot."
     PRINT #2, "NODE#   X-Displ.    Y-Displ.    Z-Rot."
     FOR I = 1 TO NN
       PRINT USING " ###"; I;
       PRINT #2, USING " ###"; I;
       I1 = 3 * I - 2: I2 = I1 + 1: I3 = I1 + 2
       PRINT USING "  ##.###^^^^"; F(I1); F(I2); F(I3)
       PRINT #2, USING "  ##.###^^^^"; F(I1); F(I2); F(I3)
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
'---- Member End-Actions
     PRINT #2, " Member End-Forces "
     FOR N = 1 TO NE
       ISTF = 1
       GOSUB ELSTIF
       I1 = NOC(N, 1): I2 = NOC(N, 2)
       FOR I = 1 TO 3
         ED(I) = F(3 * I1 - 3 + I): ED(I + 3) = F(3 * I2 - 3 + I)
       NEXT I
       FOR I = 1 TO 6
         EDP(I) = 0
         FOR K = 1 TO 6
           EDP(I) = EDP(I) + ALAMBDA(I, K) * ED(K)
         NEXT K
       NEXT I
' END FORCES DUE TO DISTRIBUTED LOADS
       IF ABS(UDL(N)) > 0 THEN
         ED(1) = 0: ED(4) = 0
         ED(2) = -UDL(N) * EL / 2: ED(5) = ED(2)
         ED(3) = -UDL(N) * EL ^ 2 / 12: ED(6) = -ED(3)
       ELSE
         FOR K = 1 TO 6: ED(K) = 0: NEXT K
       END IF
       FOR I = 1 TO 6
         EF(I) = ED(I)
         FOR K = 1 TO 6
           EF(I) = EF(I) + SEP(I, K) * EDP(K)
         NEXT K
       NEXT I
       PRINT #2, " Member #"; N
       FOR I = 1 TO 2
         II = (I - 1) * 3
         PRINT #2, USING "  ##.###^^^^"; EF(II + 1); EF(II + 2); EF(II + 3)
       NEXT I
     NEXT N
     END
GETDATA:
     '===============  READ DATA  ====================
     '----- Coordinates
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N
           INPUT #1, X(N, 1), X(N, 2)
     NEXT I
     '----- Connectivity, Material, Mom_Inertia, Dummy, Distributed Load
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), ARIN(N, 1), ARIN(N, 2), UDL(N)
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
ELSTIF:
     '----- Element Stiffness Matrix -----
     I1 = NOC(N, 1): I2 = NOC(N, 2): M = MAT(N)
     X21 = X(I2, 1) - X(I1, 1)
     Y21 = X(I2, 2) - X(I1, 2)
     EL = SQR(X21 * X21 + Y21 * Y21)
     EAL = PM(M, 1) * ARIN(N, 1) / EL
     EIZL = PM(M, 1) * ARIN(N, 2) / EL
     FOR I = 1 TO 6
     FOR J = 1 TO 6
       SEP(I, J) = 0!
     NEXT J: NEXT I
     SEP(1, 1) = EAL: SEP(1, 4) = -EAL: SEP(4, 4) = EAL
     SEP(2, 2) = 12 * EIZL / EL ^ 2: SEP(2, 3) = 6 * EIZL / EL
     SEP(2, 5) = -SEP(2, 2): SEP(2, 6) = SEP(2, 3)
     SEP(3, 3) = 4 * EIZL: SEP(3, 5) = -6 * EIZL / EL: SEP(3, 6) = 2 * EIZL
     SEP(5, 5) = 12 * EIZL / EL ^ 2: SEP(5, 6) = -6 * EIZL / EL
     SEP(6, 6) = 4 * EIZL
     FOR I = 1 TO 6
     FOR J = I TO 6
       SEP(J, I) = SEP(I, J)
     NEXT J: NEXT I
'CONVERT ELEMENT STIFFNESS MATRIX TO GLOBAL SYSTEM
     DCOS(1, 1) = X21 / EL: DCOS(1, 2) = Y21 / EL: DCOS(1, 3) = 0
     DCOS(2, 1) = -DCOS(1, 2): DCOS(2, 2) = DCOS(1, 1): DCOS(2, 3) = 0
     DCOS(3, 1) = 0: DCOS(3, 2) = 0: DCOS(3, 3) = 1
     FOR I = 1 TO 6
     FOR J = 1 TO 6
       ALAMBDA(I, J) = 0!
     NEXT J: NEXT I
     FOR K = 1 TO 2
       IK = 3 * (K - 1)
       FOR I = 1 TO 3
       FOR J = 1 TO 3
         ALAMBDA(I + IK, J + IK) = DCOS(I, J)
       NEXT J: NEXT I
     NEXT K
     IF ISTF = 1 THEN RETURN
     FOR I = 1 TO 6
     FOR J = 1 TO 6
       SE(I, J) = 0
       FOR K = 1 TO 6
         SE(I, J) = SE(I, J) + SEP(I, K) * ALAMBDA(K, J)
       NEXT K
     NEXT J: NEXT I
     FOR I = 1 TO 6: FOR J = 1 TO 6: SEP(I, J) = SE(I, J): NEXT J: NEXT I
     FOR I = 1 TO 6
     FOR J = 1 TO 6
       SE(I, J) = 0
       FOR K = 1 TO 6
         SE(I, J) = SE(I, J) + ALAMBDA(K, I) * SEP(K, J)
       NEXT K
     NEXT J: NEXT I
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

