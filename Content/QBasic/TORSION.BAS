     '***************************************
     '*         PROGRAM  TORSION            *
     '*   TORSION WITH 3-NODED TRIANGLES    *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "TORSION WITH TRIANGLES  "; SPACE$(23);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     PRINT : INPUT "File Name for Input Data ", FILE1$
     OPEN FILE1$ FOR INPUT AS #2
     LINE INPUT #2, D$: INPUT #2, TITLE$: LINE INPUT #2, D$
     INPUT #2, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #2, D$
     INPUT #2, ND, NL, NMPC
     NPR = 1  'Material Property  ThermalConductivity
     NMPC = 0: NM = 1
     '--- ND = NO. OF SPECIFIED STRESS FUNCTION VALUES
     '--- NL = NO. OF GENERALIZED NODAL FORCES "0" HERE
     '--- NPR =1 (SHEAR MODULUS) AND NMPC = 0 FOR THIS PROGRAM
     '--- ELEMENT CHARACTERISTIC NOT USED
     '--- NO. OF MATERIALS = 1 FOR THIS PROGRAM
     INPUT "File Name for Output ", FILE2$
     DIM X(NN, 2), NOC(NE, 3), MAT(NE), PM(NM, NPR), F(NN)
     DIM NU(ND), U(ND)
     OPEN FILE2$ FOR OUTPUT AS #1
     PRINT : PRINT "      PLOT CHOICE"
     PRINT "  1) No Plot Data"
     PRINT "  2) Create Data File Containing Stress Function Values"
     INPUT "     Choose 1 or 2 "; IPL
     IF IPL < 1 OR IPL > 2 THEN IPL = 1   '--- default is no data
     IF IPL = 2 THEN
       INPUT "File Name for Plot Data "; FILE3$
       OPEN FILE3$ FOR OUTPUT AS #3
     END IF
     '----- Coordinates
     LINE INPUT #2, D$
     FOR I = 1 TO NN
        INPUT #2, N
        FOR J = 1 TO NDIM
           INPUT #2, X(N, J)
        NEXT J
     NEXT I
     '----- Connectivity, Material#, Dummy
     LINE INPUT #2, D$
     FOR I = 1 TO NE
        INPUT #2, N
        FOR J = 1 TO NEN
           INPUT #2, NOC(N, J)
        NEXT J
        INPUT #2, MAT(N)
     NEXT I
     '----- Boundary Conditions
     LINE INPUT #2, D$
     FOR I = 1 TO ND: INPUT #2, NU(I), U(I): NEXT I
     '---- DUMMY READ
     LINE INPUT #2, D$
     '----- Shear Modulus of Material
     LINE INPUT #2, D$
     FOR I = 1 TO NM
        INPUT #2, N
        FOR J = 1 TO NPR
           INPUT #2, PM(N, J)
        NEXT J
     NEXT I
     '--- BAND WIDTH CALCULATION ---
     IDIFF = 0
     FOR K = 1 TO NE
        FOR I = 1 TO 2
           FOR J = I + 1 TO 3
              II = ABS(NOC(K, J) - NOC(K, I))
              IF II > IDIFF THEN IDIFF = II
           NEXT J
        NEXT I
     NEXT K
     NBW = IDIFF + 1
     PRINT #1, "THE BAND WIDTH IS"; NBW
     '--- INITIALIZATION OF STIFFNESS MATRIX
     DIM S(NN, NBW)
     FOR I = 1 TO NN
        F(I) = 0
        FOR J = 1 TO NBW
           S(I, J) = 0
        NEXT J
     NEXT I
     CLOSE #2
     '--- STIFFNESS MATRIX
     DIM BT(2, 3)
     FOR I = 1 TO NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        AREA = .5 * ABS(DETJ)
        '--- LOAD CALCULATION
        C = 2 * AREA / 3
        F(I1) = F(I1) + C: F(I2) = F(I2) + C: F(I3) = F(I3) + C
        '--- STIFFNESS FORMATION
        BT(1, 1) = Y23: BT(1, 2) = Y31: BT(1, 3) = Y12
        BT(2, 1) = X32: BT(2, 2) = X13: BT(2, 3) = X21
        FOR II = 1 TO 3
           FOR JJ = 1 TO 2
              BT(JJ, II) = BT(JJ, II) / DETJ
           NEXT JJ
        NEXT II
        FOR II = 1 TO 3
           FOR JJ = 1 TO 3
              II1 = NOC(I, II): II2 = NOC(I, JJ)
              IF II1 <= II2 THEN
                 SUM = 0
                 FOR J = 1 TO 2
                    SUM = SUM + BT(J, II) * BT(J, JJ)
                 NEXT J
                 S(II1, II2 - II1 + 1) = S(II1, II2 - II1 + 1) + SUM * AREA
              END IF
           NEXT JJ
        NEXT II
     NEXT I
     '--- MODIFY FOR BOUNDARY CONDITIONS
     CNST = S(1, 1)
     FOR I = 2 TO NN
       IF CNST < S(I, 1) THEN CNST = S(I, 1)
     NEXT I
     CNST = CNST * 1000000
     FOR I = 1 TO ND
       N = NU(I)
       S(N, 1) = S(N, 1) + CNST
       F(N) = F(N) + CNST * U(I): NEXT I
     '--- EQUATION SOLVING
     N = NN
     GOSUB BANSOL
     PRINT #1, "NODE#   STRESS FUNCTION VALUE"
     PRINT : PRINT "NODE#   STRESS FUNCTION VALUE"
     FOR I = 1 TO NN
     PRINT #1, I; TAB(10); F(I)
     PRINT I; TAB(10); F(I)
     NEXT I
     IF IPL = 2 THEN
       PRINT #3, "Stress Function Value"
       FOR I = 1 TO NN: PRINT #3, F(I): NEXT I
       CLOSE #3
       PRINT
       PRINT : PRINT "Stress Function Values Data in file "; FILE3$
       PRINT "Run CONTOURA or CONTOURB to plot isotherms"
     END IF
     '---- ANGLE OF TWIST PER UNIT LENGTH
     SUM = 0
     FOR I = 1 TO NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        SUM = SUM + ABS(DETJ) / 3 * (F(I1) + F(I2) + F(I3))
     NEXT I
     PRINT
     INPUT "TORQUE = "; TORQUE
     INPUT "SYMMETRY FACTOR (eg. if 1/4 symmetry, then =4.0) ="; SFAC
     SMOD = PM(1, 1)
     ALPHA = TORQUE / SMOD / SUM / SFAC
     PRINT #1, "TWIST PER UNIT LENGTH = "; ALPHA
     PRINT "TWIST PER UNIT LENGTH = "; ALPHA
     PRINT #1, " ** SHEARING STRESSES TAUYZ, TAUXZ IN EACH ELEMENT ** "
     PRINT #1, "ELEMENT#    TAUYZ           TAUXZ "
     PRINT "ELEMENT#    TAUYZ           TAUXZ "
     FOR I = 1 TO NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        BT(1, 1) = Y23: BT(1, 2) = Y31: BT(1, 3) = Y12
        BT(2, 1) = X32: BT(2, 2) = X13: BT(2, 3) = X21
        FOR II = 1 TO 3
           FOR JJ = 1 TO 2
              BT(JJ, II) = BT(JJ, II) / DETJ
           NEXT JJ
        NEXT II
        TAUYZ = -(BT(1, 1) * F(I1) + BT(1, 2) * F(I2) + BT(1, 3) * F(I3))
        TAUXZ = BT(2, 1) * F(I1) + BT(2, 2) * F(I2) + BT(2, 3) * F(I3)
        TAUYZ = TAUYZ * SMOD * ALPHA
        TAUXZ = TAUXZ * SMOD * ALPHA
        PRINT #1, I; TAB(12); TAUYZ; TAB(28); TAUXZ
        PRINT I; TAB(12); TAUYZ; TAB(28); TAUXZ
     NEXT I
     END
BANSOL:
     '--- FORWARD ELIMINATION
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
     '--- BACK SUBSTITUTION ---
     F(N) = F(N) / S(N, 1)
     FOR II = 1 TO N - 1
        I = N - II
        NBI = N - I + 1
        IF N - I + 1 > NBW THEN NBI = NBW
        SUM = 0
        FOR J = 2 TO NBI
           SUM = SUM + S(I, J) * F(I + J - 1)
        NEXT J
        F(I) = (F(I) - SUM) / S(I, 1)
     NEXT II
     '--- F CONTAINS THE SOLUTION. 'S' IS OVER-WRITTEN
     RETURN

