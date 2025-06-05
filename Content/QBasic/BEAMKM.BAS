     '********    PROGRAM BEAMKM    ********"
     '*      Shaft Vibration Analysis      *"
     '* T.R.Chandrupatla and A.D.Belegundu *"
     '**************************************"
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "SHAFT VIBRATION ANALYSIS"; SPACE$(27);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "File Name for Input Data ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 2      'Material Properties E, Density(Rho)
     INPUT "File Name for Output ", FILE2$
     '----- Total dof is  NQ
     NQ = 2 * NN
     DIM X(NN), NOC(NE, NEN), MAT(NE), PM(NM, NPR), SMI(NE)
     DIM NU(ND), U(ND), F(NQ), MPC(NMPC, 2), BT(NMPC, 3)
     DIM AREA(NE), SE(4, 4), EM(4, 4)
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
     DIM S(NQ, NBW), GM(NQ, NBW)
'----- Global Stiffness Matrix
     FOR N = 1 TO NE
        '--- Element Stiffness and Mass Matrices
        PRINT "Stiffness and Mass Matrices of Element "; N
        GOSUB ELKM
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
                       GM(NR, NC) = GM(NR, NC) + EM(I, J)
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
'----- Modify for Stiffness for Boundary Conditions -----
     '--- Displacement BC ---
     FOR I = 1 TO ND
        N = NU(I)
        S(N, 1) = S(N, 1) + CNST
     NEXT I
     '--- Multi-point Constraints ---
     FOR I = 1 TO NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1)
        S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2)
        IR = I1: IF IR > I2 THEN IR = I2
        IC = ABS(I2 - I1) + 1
        S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2)
     NEXT I
'-----  Additional Springs and Lumped Masses  -----
     CLS : PRINT : PRINT "SPRING SUPPORTS  < DOF# = 0 Exits this mode >"
     PRINT "     DOF#    Spring Const ": VIEW PRINT 5 TO 25
     DO
        IR = CSRLIN
        LOCATE IR, 7: INPUT "", N: IF N = 0 THEN EXIT DO
        LOCATE IR, 16: INPUT "", C
        S(N, 1) = S(N, 1) + C
     LOOP
     VIEW PRINT 2 TO 25: CLS
     CLS : PRINT : PRINT "LUMPED MASSES  < DOF# = 0 Exits this mode >"
     PRINT "     DOF#    Lumped Mass ": VIEW PRINT 5 TO 25
     DO
        IR = CSRLIN
        LOCATE IR, 7: INPUT "", N: IF N = 0 THEN EXIT DO
        LOCATE IR, 16: INPUT "", C
        GM(N, 1) = GM(N, 1) + C
     LOOP
     VIEW PRINT 2 TO 25: CLS : PRINT
     '--- Print Banded Stiffness and Mass Matrices in Output File
     OPEN FILE2$ FOR OUTPUT AS #2
     PRINT #2, "Stiffness and Mass for Data in File "; FILE1$
     PRINT #2, "Num. of DOF    Bandwidth"
     PRINT #2, NQ, NBW
     PRINT #2, "Banded Stiffness Matrix"
     FOR I = 1 TO NQ
        FOR J = 1 TO NBW
           PRINT #2, S(I, J);
        NEXT J: PRINT #2, : NEXT I
     PRINT #2, "Banded Mass Matrix"
     FOR I = 1 TO NQ
        FOR J = 1 TO NBW
           PRINT #2, GM(I, J);
        NEXT J: PRINT #2, : NEXT I
     PRINT #2, "Starting Vector for Inverse Iteration "
     FOR I = 1 TO NQ
         PRINT #2, "1 ";
     NEXT I: PRINT #2,
     CLOSE #2
     PRINT "Global Stiffness and Mass Matrices are in file "; FILE2$
     PRINT "Run INVITR or JACOBI or GENEIGEN to get Eigenvalues and"
     PRINT "Eigenvectors"
     END
GETDATA:
     '===============  READ DATA  ====================
     '----- Coordinates
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N
           INPUT #1, X(N)
     NEXT I
     '----- Connectivity, Material, Mom_Inertia, Area
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), SMI(N), AREA(N)
     NEXT I
     '----- Displacement BC
     LINE INPUT #1, D$
     FOR I = 1 TO ND: INPUT #1, NU(I), U(I): NEXT I
     '----- Component Loads
     LINE INPUT #1, D$
     FOR I = 1 TO NL: INPUT #1, N, C: NEXT I
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
ELKM:
'--------  Element Stiffness and Mass Matrices  --------
      N1 = NOC(N, 1)
      N2 = NOC(N, 2)
      M = MAT(NE)
      EL = ABS(X(N1) - X(N2))
      '--- Element Stiffness
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
      '--- Element Mass
      RHO = PM(M, 2)
      C1 = RHO * AREA(N) * EL / 420
      EM(1, 1) = 156 * C1
      EM(1, 2) = 22 * EL * C1
      EM(1, 3) = 54 * C1
      EM(1, 4) = -13 * EL * C1
         EM(2, 1) = EM(1, 2)
         EM(2, 2) = 4 * EL * EL * C1
         EM(2, 3) = 13 * EL * C1
         EM(2, 4) = -3 * EL * EL * C1
      EM(3, 1) = EM(1, 3)
      EM(3, 2) = EM(2, 3)
      EM(3, 3) = 156 * C1
      EM(3, 4) = -22 * EL * C1
         EM(4, 1) = EM(1, 4)
         EM(4, 2) = EM(2, 4)
         EM(4, 3) = EM(3, 4)
         EM(4, 4) = 4 * EL * EL * C1
     RETURN

