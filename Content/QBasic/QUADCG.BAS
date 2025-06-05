'********     PROGRAM QuadCG        **********
'*     2-D STRESS ANALYSIS USING 4-NODE      *
'*  QUADRILATERAL ELEMENTS WITH TEMPERATURE  *
'*         CONJUGATE GRADIENT APPROACH       *
'*    T.R.Chandrupatla and A.D.Belegundu     *
'*********************************************
DECLARE SUB ElemStiffness (N%)
DECLARE SUB Stiffness ()
DECLARE SUB ModifyForceForBC ()
DECLARE SUB CGSolve ()
DECLARE SUB StressCalc ()
DECLARE SUB ReactionCalc ()
DECLARE SUB DbMat (N%, ISTR%, IP%)
DECLARE SUB InputData ()
DECLARE SUB IntegPoints ()
DECLARE SUB DMatrix (N%)
DECLARE SUB FileOut ()
DEFINT I-N
DEFDBL A-H, O-Z
COMMON SHARED NN, NE, NM, NDIM, NEN, NDN
COMMON SHARED ND, NL, NCH, NPR, NMPC, NBW
COMMON SHARED X(), NOC(), MAT(), PM()
COMMON SHARED TH(), DT(), NU(), U(), F(), DD(), AD(), GG()
COMMON SHARED D(), B(), DB(), S(), Q(), QE(), STR()
COMMON SHARED TL(), XNI(), A(), G(), MPC(), BT()
COMMON SHARED React(), vonMisesStress(), ShearStress()
COMMON SHARED CNST, LC, IPL, NQ, DJ
COMMON SHARED File1 AS STRING, File2 AS STRING, File3 AS STRING
COMMON SHARED Title AS STRING, Dummy AS STRING
     CLS
     COLOR 1, 3: LOCATE 1, 1
     PRINT "CONJUGATE GRADIENT METHOD - QUAD ELEMENT"; SPACE$(10);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT

     CALL InputData
     CALL Stiffness
     CALL ModifyForceForBC
     CALL CGSolve
     CALL StressCalc
     CALL ReactionCalc
     CALL FileOut
END

SUB CGSolve
     '----- Conjugate Gradient Method starts here
     GG1 = 0
     FOR I = 1 TO NQ
        GG(I) = -F(I): DD(I) = F(I): Q(I) = 0
        GG1 = GG1 + GG(I) * GG(I)
     NEXT I
     '-----  ITERATION LOOP
     ITER = 0
    DO
     ITER = ITER + 1
     PRINT "ITERATION# = "; ITER
     '=====  ELEMENT LOOP  =====
     FOR N = 1 TO NE
        FOR I = 1 TO 4
           IGT = 2 * (NOC(N, I) - 1): ILT = 2 * (I - 1)
           FOR II = 1 TO 2
              IG = IGT + II: IL = ILT + II
              FOR J = 1 TO 4
                 JGT = 2 * (NOC(N, J) - 1): JLT = 2 * (J - 1)
                 FOR JJ = 1 TO 2
                    JG = JGT + JJ: JL = JLT + JJ
                    AD(IG) = AD(IG) + S(N, IL, JL) * DD(JG)
                 NEXT JJ
              NEXT J
           NEXT II
        NEXT I
     NEXT N
     '--- Displacement BC ---
     FOR I = 1 TO ND
        N = NU(I)
        AD(N) = AD(N) + CNST * DD(N)
     NEXT I
     '--- Multi-point Constraints ---
     FOR I = 1 TO NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        C = BT(I, 1) * DD(I1) + BT(I, 2) * DD(I2)
        AD(I1) = AD(I1) + CNST * BT(I, 1) * C
        AD(I2) = AD(I2) + CNST * BT(I, 2) * C
     NEXT I
        DAD = 0
        FOR I = 1 TO NQ
           DAD = DAD + DD(I) * AD(I)
        NEXT I
        AL = GG1 / DAD
        GG2 = 0
        FOR I = 1 TO NQ
           GG(I) = GG(I) + AL * AD(I)
           Q(I) = Q(I) + AL * DD(I)
           GG2 = GG2 + GG(I) * GG(I)
        NEXT I
        IF GG2 < .00000001# THEN EXIT DO
        BTA = GG2 / GG1
        GG1 = GG2
        FOR I = 1 TO NQ
           DD(I) = -GG(I) + BTA * DD(I)
        NEXT I
        FOR I = 1 TO NQ: AD(I) = 0: NEXT I
    LOOP

END SUB

SUB DbMat (N, ISTR, IP)
'-------  DB()  MATRIX  ------
     XI = XNI(IP, 1): ETA = XNI(IP, 2)
     '--- Nodal Coordinates
     THICK = TH(N)
     N1 = NOC(N, 1): N2 = NOC(N, 2)
     N3 = NOC(N, 3): N4 = NOC(N, 4)
     X1 = X(N1, 1): Y1 = X(N1, 2)
     X2 = X(N2, 1): Y2 = X(N2, 2)
     X3 = X(N3, 1): Y3 = X(N3, 2)
     X4 = X(N4, 1): Y4 = X(N4, 2)
     '--- Formation of Jacobian  TJ
     TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4
     TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4
     TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4
     TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4
     '--- Determinant of the JACOBIAN
     DJ = TJ11 * TJ22 - TJ12 * TJ21
     '--- A(3,4) Matrix relates Strains to
     '--- Local Derivatives of u
     A(1, 1) = TJ22 / DJ: A(2, 1) = 0: A(3, 1) = -TJ21 / DJ
     A(1, 2) = -TJ12 / DJ: A(2, 2) = 0: A(3, 2) = TJ11 / DJ
     A(1, 3) = 0: A(2, 3) = -TJ21 / DJ: A(3, 3) = TJ22 / DJ
     A(1, 4) = 0: A(2, 4) = TJ11 / DJ: A(3, 4) = -TJ12 / DJ
     '--- G(4,8) Matrix relates Local Derivatives of u
     '--- to Local Nodal Displacements q(8)
     FOR I = 1 TO 4: FOR J = 1 TO 8
     G(I, J) = 0: NEXT J: NEXT I
     G(1, 1) = -(1 - ETA) / 4: G(2, 1) = -(1 - XI) / 4
     G(3, 2) = -(1 - ETA) / 4: G(4, 2) = -(1 - XI) / 4
     G(1, 3) = (1 - ETA) / 4: G(2, 3) = -(1 + XI) / 4
     G(3, 4) = (1 - ETA) / 4: G(4, 4) = -(1 + XI) / 4
     G(1, 5) = (1 + ETA) / 4: G(2, 5) = (1 + XI) / 4
     G(3, 6) = (1 + ETA) / 4: G(4, 6) = (1 + XI) / 4
     G(1, 7) = -(1 + ETA) / 4: G(2, 7) = (1 - XI) / 4
     G(3, 8) = -(1 + ETA) / 4: G(4, 8) = (1 - XI) / 4
     '--- B(3,8) Matrix Relates Strains to q
     FOR I = 1 TO 3
        FOR J = 1 TO 8
           C = 0
           FOR K = 1 TO 4
              C = C + A(I, K) * G(K, J)
           NEXT K
           B(I, J) = C
        NEXT J
     NEXT I
     '--- DB(3,8) Matrix relates Stresses to q(8)
     FOR I = 1 TO 3
        FOR J = 1 TO 8
           C = 0
           FOR K = 1 TO 3
              C = C + D(I, K) * B(K, J)
           NEXT K:
           DB(I, J) = C
        NEXT J
     NEXT I
     IF ISTR = 2 THEN
           '--- Stress Evaluation
           FOR I = 1 TO NEN
              IIN = NDN * (NOC(N, I) - 1)
              II = NDN * (I - 1)
              FOR J = 1 TO NDN
                 QE(II + J) = Q(IIN + J)
              NEXT J
           NEXT I
           AL = PM(MAT(N), 3)
           C1 = AL * DT(N): IF LC = 2 THEN C1 = C1 * (1 + PNU)
           FOR I = 1 TO 3
              C = 0
              FOR K = 1 TO 8
                 C = C + DB(I, K) * QE(K)
              NEXT K
              STR(I) = C - C1 * (D(I, 1) + D(I, 2))
           NEXT I
     END IF
END SUB

SUB DMatrix (N)
'-----  D() Matrix  -----
     '--- Material Properties
     MATN = MAT(N)
     E = PM(MATN, 1): PNU = PM(MATN, 2)
     AL = PM(MATN, 3)
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
END SUB

SUB ElemStiffness (N)
'--------  Element Stiffness and Temperature Load  -----
     FOR I = 1 TO 8: TL(I) = 0: NEXT I
     DTE = DT(N)
     '--- Weight Factor is ONE
     '--- Loop on Integration Points
     FOR IP = 1 TO 4
        '---  Get DB Matrix at Integration Point IP
        CALL DbMat(N, 1, IP)
        '--- Element Stiffness Matrix  SE
        FOR I = 1 TO 8
           FOR J = 1 TO 8
              C = 0
              FOR K = 1 TO 3
                 C = C + B(K, I) * DB(K, J) * DJ * TH(N)
              NEXT K
              S(N, I, J) = S(N, I, J) + C
           NEXT J
        NEXT I
        '--- Determine Temperature Load TL
        AL = PM(MAT(N), 3)
        C = AL * DTE: IF LC = 2 THEN C = (1 + PNU) * C
        FOR I = 1 TO 8
           TL(I) = TL(I) + TH(N) * DJ * C * (DB(1, I) + DB(2, I))
        NEXT I
     NEXT IP
END SUB

SUB FileOut
     '===== Print Displacements, Stresses, and Reactions
     INPUT "Output File d:\dir\fileName.ext ", File2
     OPEN File2 FOR OUTPUT AS #2
     PRINT #2, "Program Quad - CHANDRUPATLA & BELEGUNDU"
     PRINT #2, "Output for Data from " + File1
     PRINT #2, Title
     IF LC = 1 THEN PRINT #2, "Plane Stress Analysis"
     IF LC = 2 THEN PRINT #2, "Plane Strain Analysis"
     '----- Displacements -----
     PRINT "NODE#   X-Displ     Y-Displ"
     PRINT #2, "NODE#   X-Displ     Y-Displ"
     FOR I = 1 TO NN
        PRINT USING " ###"; I;
        PRINT USING "  ##.####^^^^"; Q(2 * I - 1); Q(2 * I)
        PRINT #2, USING " ###"; I;
        PRINT #2, USING "  ##.####^^^^"; Q(2 * I - 1); Q(2 * I)
     NEXT I
     '----- Stress Output -----
     PRINT #2, "ELEM#    vonMises Stresses at 4 Integration points"
     FOR N = 1 TO NE
        PRINT #2, USING " ###"; N;
        FOR IP = 1 TO 4
        PRINT #2, USING "  ##.###^^^^"; vonMisesStress(N, IP);
        NEXT IP
        PRINT #2,
     NEXT N
     PRINT "Complete results are in file "; File2
  IF IPL > 1 THEN
     INPUT "Plot Data File d:\dir\fn ", File3
     OPEN File3 FOR OUTPUT AS #3
     IF IPL = 2 THEN PRINT #3, "Max. in-plane Shear Stress "
     IF IPL = 3 THEN PRINT #3, "Von Mises Stress "
     FOR N = 1 TO NE
        FOR IP = 1 TO 4
           IF IPL = 2 THEN PRINT #3, ShearStress(N, IP);
           IF IPL = 3 THEN PRINT #3, vonMisesStress(N, IP);
        NEXT IP
        PRINT #3,
     NEXT N
     PRINT "Element Stess Data in "; File3
     PRINT "Run BESTFITQ and then CONTOURA or CONTOURB to plot stresses"
     CLOSE #3
  END IF
     '----- Reactions -----
     PRINT #2, "DOF#  REACTION"
     FOR I = 1 TO ND
        N = NU(I)
        PRINT #2, USING " ###"; N;
        PRINT #2, "##.###^^^^ "; React(I)
     NEXT I
     CLOSE #2
     PRINT "Complete results are in file "; File2
END SUB

SUB InputData
     INPUT "Input File d:\dir\fileName.ext ", File1
     PRINT "1. Plane Stress Analysis"
     PRINT "2. Plane Strain Analysis"
     INPUT "   Your Choice <1 or 2> "; LC
     IF LC < 1 OR LC > 2 THEN LC = 1   '--- default is Plane Stress
     NPR = 3   'Material Properties E, Nu, Alpha
     PRINT " 1) No Plot Data"
     PRINT " 2) Create Data File for in-plane Shear Stress"
     PRINT " 3) Create Data File for Von Mises Stress"
     INPUT "    Choose 1 or 2 or 3"; IPL
     OPEN File1 FOR INPUT AS #1
     LINE INPUT #1, Dummy: INPUT #1, Title
     LINE INPUT #1, Dummy: INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, Dummy: INPUT #1, ND, NL, NMPC
     NPR = 3      'Material Properties E, Nu, Alpha
     '----- Total dof is  NQ
     NQ = NDN * NN
     DIM X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     DIM TH(NE), DT(NE), NU(ND), U(ND), F(NQ), DD(NQ), AD(NQ), GG(NQ)
     DIM D(3, 3), B(3, 8), DB(3, 8), S(NE, 8, 8), Q(NQ), QE(8), STR(3)
     DIM TL(8), XNI(4, 2), A(3, 4), G(4, 8), MPC(NMPC, 2), BT(NMPC, 3)
     '=============  READ DATA  ===============
     '----- Coordinates
     LINE INPUT #1, Dummy
     FOR I = 1 TO NN
        INPUT #1, N
        FOR J = 1 TO NDIM
           INPUT #1, X(N, J)
        NEXT J
     NEXT I
     '----- Connectivity, Material, Thickness, Temp-change
     LINE INPUT #1, Dummy
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), TH(N), DT(N)
     NEXT I
     '----- Displacement BC
     LINE INPUT #1, Dummy
     FOR I = 1 TO ND
        INPUT #1, NU(I), U(I)
     NEXT I
     '----- Component Loads
     LINE INPUT #1, Dummy
     FOR I = 1 TO NL
        INPUT #1, N
        INPUT #1, F(N)
     NEXT I
     '----- Material Properties
     LINE INPUT #1, Dummy
     FOR I = 1 TO NM
        INPUT #1, N
        FOR J = 1 TO NPR
           INPUT #1, PM(N, J)
        NEXT J
     NEXT I
     IF NMPC > 0 THEN
        '-----  Multi-point Constraints
        LINE INPUT #1, Dummy
        FOR I = 1 TO NMPC
           INPUT #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        NEXT I
     END IF
     CLOSE #1
END SUB

SUB IntegPoints
'------- Integration Points XNI() --------
     C = .57735026919#
     XNI(1, 1) = -C: XNI(1, 2) = -C
     XNI(2, 1) = C: XNI(2, 2) = -C
     XNI(3, 1) = C: XNI(3, 2) = C
     XNI(4, 1) = -C: XNI(4, 2) = C
END SUB

SUB ModifyForceForBC
'----- GG() diagonal stiffness summation
'----- Decide Penalty Parameter CNST -----
     CNST = 0
     FOR I = 1 TO NQ
        IF CNST < GG(I) THEN CNST = GG(I)
     NEXT I
     CNST = CNST * 10000
'----- Modify right hand side F() for Boundary Conditions -----
     '--- Displacement BC ---
     FOR I = 1 TO ND
        N = NU(I)
        F(N) = F(N) + CNST * U(I)
     NEXT I
     '--- Multi-point Constraints ---
     FOR I = 1 TO NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
     NEXT I
END SUB

SUB ReactionCalc
     DIM React(ND)
     '----- Reaction Calculation -----
     FOR I = 1 TO ND
        N = NU(I)
        React(I) = CNST * (U(I) - Q(N))
     NEXT I
END SUB

SUB Stiffness
     '----- Stiffness Matrix of Each Element -----
     CALL IntegPoints
     FOR N = 1 TO NE
        PRINT "Forming Stiffness Matrix of Element "; N
        CALL DMatrix(N)
        CALL ElemStiffness(N)
        PRINT ".... Placing in Global Locations"
        '--- Add temperature loads
        FOR II = 1 TO NEN
           NRT = NDN * (NOC(N, II) - 1)
           FOR IT = 1 TO NDN
              NR = NRT + IT
              I = NDN * (II - 1) + IT
              F(NR) = F(NR) + TL(I)
              GG(NR) = GG(NR) + S(N, I, I)
           NEXT IT
        NEXT II
     NEXT N
END SUB

SUB StressCalc
     DIM vonMisesStress(NE, 4), ShearStress(NE, 4)
     '-----  Stress Calculations
     FOR N = 1 TO NE
        CALL DMatrix(N)
        FOR IP = 1 TO 4
           CALL DbMat(N, 2, IP) '--- Get DB Matrix with Stress calculation
        '--- Von Mises Stress at Integration Point
           C = 0: IF LC = 2 THEN C = PNU * (STR(1) + STR(2))
           C1 = (STR(1) - STR(2)) ^ 2 + (STR(2) - C) ^ 2 + (C - STR(1)) ^ 2
           SV = SQR(.5 * C1 + 3 * STR(3) ^ 2)
        '--- Maximum Shear Stress R
           R = SQR(.25 * (STR(1) - STR(2)) ^ 2 + (STR(3)) ^ 2)
           ShearStress(N, IP) = R
           vonMisesStress(N, IP) = SV
        NEXT IP
     NEXT N
END SUB

