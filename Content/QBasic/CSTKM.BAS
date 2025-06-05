     '***************************************
     '*          PROGRAM  CSTKM             *
     '*    STIFFNESS AND MASS MATRICES      *
     '*   2-D  CONSTANT STRAIN TRIANGLE     *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "STIFFNESS AND MASS FOR CST"; SPACE$(25);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     PRINT "  1) Plane Stress"
     PRINT "  2) Plane Strain"
     INPUT "     Choose 1 or 2 "; LC
     IF LC < 1 OR LC > 2 THEN LC = 1   '--- default is Plane Stress
     PRINT : INPUT "File Name for Input Data ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 4   'Material Properties E, Nu, Alpha, Density(Rho)
     INPUT "File Name for Output ", FILE2$
     '----- Total dof is  NQ
     NQ = NDN * NN
     DIM X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     DIM TH(NE), NU(ND), U(ND), MPC(NMPC, 2), BT(NMPC, 3), D(3, 3)
     DIM B(3, 6), DB(3, 6), SE(6, 6), EM(6, 6), Q(6), STR(3), TL(6)
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
     DIM S(NQ, NBW), GM(NQ, NBW)
'----- Global Stiffness and Mass Matrices
     FOR N = 1 TO NE
        PRINT "Forming Stiffness and Mass Matrices of Element "; N
        GOSUB DBMAT
        GOSUB ELKM             'Element Stiffness and Mass
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
'----- Modify for Stiffness Matrix for Boundary Conditions -----
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
        INPUT #1, MAT(N), TH(N), C
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
ELKM:
     '--- Element Stiffness SE()
     FOR I = 1 TO 6
        FOR J = 1 TO 6
           C = 0
           FOR K = 1 TO 3
              C = C + .5 * ABS(DJ) * B(K, I) * DB(K, J) * TH(N)
           NEXT K
           SE(I, J) = C
        NEXT J
     NEXT I
     '--- Element Mass  EM()
     RHO = PM(M, 4)
     CM = RHO * TH(N) * .5 * ABS(DJ) / 12
     FOR I = 1 TO 6: FOR J = 1 TO 6: EM(I, J) = 0: NEXT J: NEXT I
     '--- Non-zero elements of mass matrix are defined
     EM(1, 1) = 2 * CM: EM(1, 3) = CM: EM(1, 5) = CM
     EM(2, 2) = 2 * CM: EM(2, 4) = CM: EM(2, 6) = CM
     EM(3, 1) = CM: EM(3, 3) = 2 * CM: EM(3, 5) = CM
     EM(4, 2) = CM: EM(4, 4) = 2 * CM: EM(4, 6) = CM
     EM(5, 1) = CM: EM(5, 3) = CM: EM(5, 5) = 2 * CM
     EM(6, 2) = CM: EM(6, 4) = CM: EM(6, 6) = 2 * CM
     RETURN

