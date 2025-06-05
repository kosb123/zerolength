C    ******		TRUSS2D		******
C    **** Two Dimensional Truss Analysis ****
      DIMENSION X(100,2),NOC(200,4),PM(10,2),AREA(20)
      DIMENSION NU(50),U(50),DT(200),MAT(200),MPC(20,2),BT(20,3)
      DIMENSION S(200,50),F(200),SE(4,4),TL(4)
      CHARACTER*16 FILE1, FILE2
      CHARACTER*81 DUMMY, TITLE
      IMAX = 200
      PRINT *, ' ======================================='
      PRINT *, '           PROGRAM TRUSS2D               '
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
C     ----- Bandwidth Evaluation -----
      NBW = 0
      DO 1050  N = 1, NE
        NABS = NDN * (ABS(NOC(N, 1) - NOC(N, 2)) + 1)
        IF (NBW .LT. NABS) NBW = NABS
 1050 CONTINUE 
      DO 1060  I = 1, NMPC
        NABS = ABS(MPC(I, 1) - MPC(I, 2)) + 1
        IF (NBW .LT. NABS) NBW = NABS
 1060 CONTINUE
      PRINT *, 'Bandwidth = ', NBW
C    ----- First initialize stiffness matrix
      DO 1070 I = 1, NQ
      DO 1070 J = 1, NBW
 1070 S(I,J) = 0
C    -----  Stiffness Matrix -----
      DO 1090  N = 1, NE
        I1 = NOC(N, 1)
        I2 = NOC(N, 2)
        I3 = MAT(N)
        X21 = X(I2, 1) - X(I1, 1)
        Y21 = X(I2, 2) - X(I1, 2)
        EL = SQRT(X21 * X21 + Y21 * Y21)
        EAL = PM(I3, 1) * AREA(N) / EL
        CS = X21 / EL
        SN = Y21 / EL
C    ----------- Element Stiffness Matrix SE() -----------
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
C    -------------- Temperature Load TL() ---------------
        EE0 = PM(I3, 2) * DT(N) * PM(I3, 1) * AREA(N)
        TL(1) = -EE0 * CS
        TL(2) = -EE0 * SN
        TL(3) = EE0 * CS
        TL(4) = EE0 * SN
        PRINT *,  '..... Adding Stiffness to Global Locations'
	DO 1090  II = 1, NEN
           NRT = NDN * (NOC(N, II) - 1)
	   DO 1090  IT = 1, NDN
              NR = NRT + IT
              I = NDN * (II - 1) + IT
              DO 1080  JJ = 1, NEN
                 NCT = NDN * (NOC(N, JJ) - 1)
                 DO 1080  JT = 1, NDN
                    J = NDN * (JJ - 1) + JT
                    NC = NCT + JT - NR + 1
                    IF (NC .GT. 0) THEN
                       S(NR, NC) = S(NR, NC) + SE(I, J)
                    END IF
 1080           CONTINUE
              F(NR) = F(NR) + TL(I)
 1090	     CONTINUE
C     ----- Decide Penalty Parameter CNST -----
      CNST = 0
      DO 1100  I = 1, NQ
        IF (CNST .LT. S(I, 1)) CNST = S(I, 1)
 1100 CONTINUE 
      CNST = CNST * 10000
C    ----- Modify for Boundary Conditions -----
C    --- Displacement BC ---
      DO 1110  I = 1, ND
        N = NU(I)
        S(N, 1) = S(N, 1) + CNST
        F(N) = F(N) + CNST * U(I)
 1110 CONTINUE
C    --- Multi-point Constraints ---
      DO 1120  I = 1, NMPC
        I1 = MPC(I, 1)
        I2 = MPC(I, 2)
        S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1)
        S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2)
        IR = I1
        IF (IR .GT.I2) IR = I2
        IC = ABS(I2 - I1) + 1
        S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
 1120 CONTINUE
C    ----- Equation Solving using Band Solver -----
      CALL BAND(S, F, IMAX, NBW, NQ)
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
      DO 1130  I = 1, NE
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
 1130   CONTINUE
C     ----- Reaction Calculation -----
      PRINT *, 'DOF#      Reaction'
      WRITE(LOUT,'(A)')  'DOF#      Reaction'
      DO 1140  I = 1, ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
	WRITE(*,'(1X,I4,E15.4)') N, R
        WRITE(LOUT,'(1X,I4,E15.4)') N, R
 1140 CONTINUE
      CLOSE(LOUT)
      PRINT *, 'Results are now available in the text file ', FILE2
      END
      
      SUBROUTINE BAND(A, B, IMAX, NBW, N)
      DIMENSION A(IMAX,NBW), B(IMAX)
      N1 = N - 1
      PRINT *, ' *** FORWARD ELIMINATION ***'
      DO 2100 K = 1, N1
      NK = N - K + 1
      IF(NK .GT. NBW) NK = NBW
      DO 2100 I = 2, NK
      C1 = A(K, I) / A(K, 1)
      I1 = K + I - 1
      DO 2000 J = I, NK
      J1 = J - I + 1
 2000 A(I1, J1) = A(I1, J1) - C1 * A(K, J)
 2100 B(I1) = B(I1) - C1 * B(K)
      PRINT *, ' *** BACK SUBSTITUTION ***'
      B(N) = B(N) / A(N, 1)
      DO 2300 KK = 1, N1
      K = N - KK
      C1 = 1 / A(K, 1)
      B(K) = C1 * B(K)
      NK = N - K + 1
      IF (NK .GT. NBW) NK = NBW
      DO 2200 J = 2, NK
      B(K) = B(K) - C1 * A(K, J) * B(K + J - 1)
 2200 CONTINUE
 2300 CONTINUE
      RETURN
      END
