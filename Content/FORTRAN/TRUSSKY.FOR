C    ******		TRUSSKY		******
C    **** Two Dimensional Truss Analysis using Skyline Storage ****
      DIMENSION X(100,2),NOC(200,4),PM(10,2),AREA(20)
      DIMENSION NU(50),U(50),DT(200),MAT(200),MPC(20,2),BT(20,3)
      DIMENSION A(10000),F(200),SE(4,4),TL(4),ID(200)
      CHARACTER*16 FILE1, FILE2
      CHARACTER*81 DUMMY, TITLE
      IMAX = 200
      PRINT *, ' ======================================='
      PRINT *, '         PROGRAM TRUSSKY                '
      PRINT *, '   T.R.Chandrupatla and A.D.Belegundu   '
      PRINT *, ' ======================================='
      PRINT *, 'Data File Name '
      READ '(A)', FILE1
      LINP = 10
      OPEN (UNIT = 10, FILE = FILE1, STATUS = 'OLD')
      PRINT *, 'Output File Name'
      READ '(A)', FILE2
      READ(LINP,'(A)') DUMMY
      READ(LINP,'(A)') TITLE
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) NN, NE, NM, NDIM, NEN, NDN
      NQ = 2 * NN
      READ(LINP,'(A)') DUMMY      
      READ(LINP,*) ND, NL, NMPC
C     Material Properties E, Alpha
      NPR = 2
C   ============  READ DATA FROM FILE  ===========
C     ----- Coordinates -----
      READ(LINP,'(A)') DUMMY
      DO 900 I = 1, NN
  900 READ (LINP, *) N, X(N, 1), X(N, 2)
C     ----- Connectivity -----
      READ(LINP,'(A)') DUMMY
      DO 1000 I=1,NE
      READ (LINP, *) N,NOC(N,1),NOC(N,2),MAT(N),AREA(N),DT(N)
 1000 CONTINUE
C    ----- Specified Displacements -----
      READ(LINP,'(A)') DUMMY
      DO 1010  I = 1, ND
 1010 READ (LINP, *) NU(I), U(I)
C     ----- Component Loads -----
      READ(LINP,'(A)') DUMMY      
C     ----- First initialize and then read loads -----
      DO 1015 I=1, NQ
 1015 F(I) = 0
      DO 1020  I = 1, NL
 1020 READ (LINP, *) N, F(N)
C     ----- Material Properties -----
      READ(LINP,'(A)') DUMMY      
      DO 1030  I = 1, NM
 1030	READ (LINP, *) N, (PM(N, J), J = 1, NPR)
C     ----- Multi-point Constraints B1*Qi+B2*Qj=B0
      IF (NMPC .GT. 0) THEN
	READ(LINP,'(A)') DUMMY
        DO 1040  I = 1, NMPC
 1040	READ (LINP,*) BT(I,1),MPC(I,1),BT(I,2),MPC(I,2),BT(I,3)
      END IF
      CLOSE(LINP)
C     ----- Preparation for Skyline Storage -----
      DO 1050 I = 1, NE
        II = NOC(I, 1)
	IF (II .GT. NOC(I, 2)) II = NOC(I, 2)
        II1 = NDN * (II - 1) + 1
        DO 1050 J = 1, NEN
           NIJ = NDN * (NOC(I, J) - 1)
           DO 1050 JJ = 1, NDN
              NDG = NIJ + JJ
	      NHT = NDG - II1 + 1
              IF (NHT .GT. ID(NDG)) ID(NDG) = NHT
 1050 CONTINUE
C     ----- Skyline Height adjustment for MPC
      DO 1060 I = 1, NMPC
        I1 = MPC(I, 1)
	I2 = MPC(I, 2)
        NDG = I1
	IF (NDG .LT. I2) NDG = I2
        NHT = ABS(I2 - I1) + 1
	IF (NHT .GT. ID(NDG)) ID(NDG) = NHT
 1060 CONTINUE 
      DO 1070 I = 2, NQ
 1070 ID(I) = ID(I) + ID(I - 1)
      NSUM = ID(NQ)
      DO 1080 I = 1, NSUM
 1080 A(I) = 0 
C    -----  Stiffness Matrix -----
      DO 1100 N = 1, NE
        I1 = NOC(N, 1)
	I2 = NOC(N, 2)
        I3 = MAT(N)
        X21 = X(I2, 1) - X(I1, 1)
        Y21 = X(I2, 2) - X(I1, 2)
        EL = SQRT(X21 * X21 + Y21 * Y21)
        EAL = PM(I3, 1) * AREA(N) / EL
        CS = X21 / EL
	SN = Y21 / EL
C     ----------- Element Stiffness Matrix SE() -----------
        SE(1, 1) = CS * CS * EAL
        SE(1, 2) = CS * SN * EAL
	SE(2, 1) = SE(1, 2)
        SE(1, 3) = -CS * CS * EAL
	SE(3, 1) = SE(1, 3)
        SE(1, 4) = -CS * SN * EAL
	SE(4, 1) = SE(1, 4)
        SE(2, 2) = SN * SN * EAL
        SE(2, 3) = -CS * SN * EAL
	SE(3, 2) = SE(2, 3)
        SE(2, 4) = -SN * SN * EAL
	SE(4, 2) = SE(2, 4)
        SE(3, 3) = CS * CS * EAL
        SE(3, 4) = CS * SN * EAL
	SE(4, 3) = SE(3, 4)
        SE(4, 4) = SN * SN * EAL
C     -------------- Temperature Load TL() ---------------
        EE0 = PM(I3, 2) * DT(N) * PM(I3, 1) * AREA(N)
        TL(1) = -EE0 * CS
	TL(2) = -EE0 * SN
        TL(3) = EE0 * CS
	TL(4) = EE0 * SN
        PRINT *, '..... Adding Stiffness to Global Locations'
        DO 1100 II = 1, NEN
           NCT = NDN * (NOC(N, II) - 1)
           DO 1100 IT = 1, NDN
              NC = NCT + IT
	      IID = ID(NC)
              I = NDN * (II - 1) + IT
              DO 1090 JJ = 1, NEN
                 NRT = NDN * (NOC(N, JJ) - 1)
                 DO 1090 JT = 1, 2
                    J = NDN * (JJ - 1) + JT
                    NR = NRT + JT
                    IF (NR .LE. NC) THEN
                       NLOC = IID - (NC - NR)
                       A(NLOC) = A(NLOC) + SE(I, J)
                    END IF
 1090         CONTINUE 
              F(NC) = F(NC) + TL(I)
 1100      CONTINUE
C    ----- Decide Penalty Parameter CNST -----
      CNST = 0
      DO 1110 I = 1, NQ
        II = ID(I)
        IF (CNST .LT. A(II)) CNST = A(II)
 1110 CONTINUE
      CNST = CNST * 10000
C    ----- Modify for Boundary Conditions -----
C     --- Displacement BC ---
      DO 1120 I = 1, ND
        N = NU(I)
	II = ID(N)
        A(II) = A(II) + CNST
 1120 F(N) = F(N) + CNST * U(I)
C    --- Multi-point Constraints ---
      DO 1130 I = 1, NMPC
        I1 = MPC(I, 1)
	I2 = MPC(I, 2)
        A(ID(I1)) = A(ID(I1)) + CNST * BT(I, 1) * BT(I, 1)
        A(ID(I2)) = A(ID(I2)) + CNST * BT(I, 2) * BT(I, 2)
        II = I1
	IF (II .LT. I2) II = I2
        IL = ID(II) - ABS(I2 - I1)
        A(IL) = A(IL) + CNST * BT(I, 1) * BT(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
 1130 CONTINUE 
C   ----- Equation Solving using Skyline Solver -----
      CALL SKYLIN(A, F, ID, NSUM, NQ)
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)      
      WRITE(LOUT,'(1X,2A)') 'Output for Input Data from File ' ,FILE1
      PRINT *, 'Output for Input Data from File ' ,FILE1
      WRITE(LOUT,'(1X,A)')  TITLE
      PRINT *, TITLE
      PRINT *, 'Node#    X-Displ       Y-Displ'
      WRITE(LOUT,'(1X,A)') 'Node#    X-Displ       Y-Displ'
      WRITE(LOUT,'(1X,I4,2E15.4)') (I,F(2*I-1),F(2*I),I=1,NN)
      WRITE(*,'(1X,I4,2E15.4)') (I,F(2*I-1),F(2*I),I=1,NN)
C     ----- Stress Calculation -----
      PRINT *, 'Elem#    Stress'
      WRITE(LOUT,'(1X,A)') ' Elem#    Stress'
      DO 1140  I = 1, NE
        I1 = NOC(I, 1)
        I2 = NOC(I, 2)
        I3 = MAT(I)
        X21 = X(I2, 1) - X(I1, 1)
        Y21 = X(I2, 2) - X(I1, 2)
        EL = SQRT(X21 * X21 + Y21 * Y21)
        CS = X21 / EL
        SN = Y21 / EL
        J2 = 2 * I1
        J1 = J2 - 1
        K2 = 2 * I2
        K1 = K2 - 1
        DLT = (F(K1) - F(J1)) * CS + (F(K2) - F(J2)) * SN
        STRESS = PM(I3, 1) * (DLT / EL - PM(I3, 2) * DT(I))
	    WRITE(*, '(1X,I4,E15.4)') I, STRESS
        WRITE(LOUT, '(1X,I4,E15.4)') I, STRESS
 1140   CONTINUE
C     ----- Reaction Calculation -----
      PRINT *, 'DOF#      Reaction'
      WRITE(LOUT,'(1X,A)')  'DOF#      Reaction'
      DO 1150  I = 1, ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
	    WRITE(*,'(1X,I4,E15.4)') N, R
        WRITE(LOUT,'(1X,I4,E15.4)') N, R
 1150 CONTINUE
      CLOSE(LOUT)
      PRINT *, 'Results are now available in the text file ', FILE2
      END

      SUBROUTINE SKYLIN(A, F, ID, NSUM, NQ)
C     -----  Sky Line Solver -----
      DIMENSION A(NSUM), F(NQ), ID(NQ)
C     --- Forward Elimination ---
      DO 2010 J = 2, NQ
        NJ = ID(J) - ID(J - 1)
        IF (NJ .NE. 1) THEN
           K1 = 0
	   NJ = J - NJ + 1
           DO 2000 K = NJ, J - 1
              K1 = K1 + 1
	      KJ = ID(J - 1) + K1
	      KK = ID(K)
              C = A(KJ) / A(KK)
              DO 2000 I = K + 1, J
                 NI = ID(I) - ID(I - 1)
                 IF ((I - K + 1) .LE. NI) THEN
                    IJ = ID(J) - J + I
                    KI = ID(I) - I + K
                    A(IJ) = A(IJ) - C * A(KI)
                 END IF
 2000      CONTINUE
        END IF
 2010 CONTINUE
      DO 2020 K = 1, NQ - 1
	KK = ID(K)
        C = F(K) / A(KK)
        DO 2020 I = K + 1, NQ
           NI = ID(I) - ID(I - 1)
           IF ((I - K + 1) .LE. NI) THEN
              KI = ID(I) - I + K
              F(I) = F(I) - C * A(KI)
           END IF
 2020 CONTINUE
C    --- Back-substitution ---
      NS = ID(NQ)
      F(NQ) = F(NQ) / A(NS)
      DO 2030 I1 = 1, NQ - 1
        I = NQ - I1
        II = ID(I)
        C = 1 / A(II)
	F(I) = C * F(I)
        DO 2030 J = I + 1, NQ
           J1 = J - I + 1
	   NJ = ID(J) - ID(J - 1)
           IF (J1 .LE. NJ) THEN
              IJ = ID(J) - J + I
              F(I) = F(I) - C * A(IJ) * F(J)
           END IF
 2030 CONTINUE 
      RETURN
      END
