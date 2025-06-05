C   ********           FRAME3D           ********
C   *****     3-D FRAME ANALYSIS BY FEM     *****
      DIMENSION X(150,3),NOC(100,3),MAT(100),PM(10,5),ARIN(100,4)
      DIMENSION NU(40),U(40),S(600,90),F(600),SEP(12,12),SE(12,12)
      DIMENSION MPC(30,2),BT(30,3),UDL(100,2),ALMBDA(12,12)
      DIMENSION ED(12),EDP(12),EF(12)
      CHARACTER*16 FILE1, FILE2
      CHARACTER*81 DUMMY, TITLE
      IMAX = 600
      PRINT *, '********    PROGRAM FRAME3D      ********'
      PRINT *, '*      3-D   FRAME ANALYSIS BY FEM      *'
      PRINT *, '*  T.R.Chandrupatla and A.D.Belegundu   *'
      PRINT *, '*****************************************'
      PRINT *, 'Data File Name '
      READ '(A)', FILE1
      LINP = 10
      OPEN (UNIT = 10, FILE = FILE1, STATUS = 'OLD')
      PRINT *, 'Output File Name'
      READ '(A)', FILE2
      READ(LINP,'(A)') DUMMY
      READ(LINP,'(A)') TITLE
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) NN, NE, NM, NDIM, NEN, NDN, NNREF
      NNT = NN + NNREF
      NQ = NDN * NN
      READ(LINP,'(A)') DUMMY      
      READ(LINP,*) ND, NL, NMPC
C    --- Material properties E, G (shear modulus)
      NPR = 2
C   ============  READ DATA FROM FILE  ===========
C     ----- Coordinates -----
      READ(LINP,'(A)') DUMMY
      DO 900 I = 1, NNT
  900 READ (LINP, *) N, X(N, 1), X(N, 2), X(N, 3)
C     ----- Connectivity  Ref_Pt Mat# Iy  Iz  J  UDLy'  UDLz'  -----
      READ(LINP,'(A)') DUMMY
      DO 1000 I=1,NE
      READ (LINP, *) N,NOC(N,1),NOC(N,2),NOC(N,3),MAT(N),ARIN(N,1),
     1ARIN(N,2),ARIN(N,3),ARIN(N,4),UDL(N,1),UDL(N,2)
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
C     ----- Global Stiffness Matrix
      DO 1090 N = 1, NE
        PRINT *, 'Forming Stiffness Matrix of Element ', N
        ISTF = 2
	  CALL ELSTIF(N,NOC,X,MAT,PM,ARIN,ALMBDA,SE,SEP,EL,ISTF)
	   PRINT *, '.... Placing in Global Locations'
	   DO 1090  II = 1, NEN
           NRT = NDN * (NOC(N, II) - 1)
	       DO 1090  IT = 1, NDN
              NR = NRT + IT
              I = NDN * (II - 1) + IT
              DO 1090  JJ = 1, NEN
                 NCT = NDN * (NOC(N, JJ) - 1)
                 DO 1090  JT = 1, NDN
                    J = NDN * (JJ - 1) + JT
                    NC = NCT + JT - NR + 1
                    IF (NC .GT. 0) THEN
                       S(NR, NC) = S(NR, NC) + SE(I, J)
                    END IF
 1090	     CONTINUE
C     ----- Decide Penalty Parameter CNST -----
      CNST = 0
      DO 1100  I = 1, NQ
        IF (CNST .LT. S(I, 1)) CNST = S(I, 1)
 1100 CONTINUE
      CNST = CNST * 10000
C    ----- Loads due to uniformly distributed load on element
      DO 1130 N = 1, NE
      IF ((ABS(UDL(N, 1)) .GT. 0) .OR. (ABS(UDL(N, 2)) .GT. 0)) THEN
       ISTF = 1
       CALL ELSTIF(N,NOC,X,MAT,PM,ARIN,ALMBDA,SE,SEP,EL,ISTF)
       I1 = NOC(N, 1)
       I2 = NOC(N, 2)
       ED(1) = 0
       ED(4) = 0
       ED(7) = 0
       ED(10) = 0
       ED(2) = UDL(N, 1) * EL / 2
       ED(8) = ED(2)
       ED(6) = UDL(N, 1) * EL**2 / 12
       ED(12) = -ED(6)
       ED(3) = UDL(N, 2) * EL / 2
       ED(9) = ED(3)
       ED(5) = -UDL(N, 2) * EL**2 / 12
       ED(11) = -ED(5)
       DO 1110 I = 1, 12
         EDP(I) = 0
         DO 1110 K = 1, 12
	   EDP(I) = EDP(I) + ALMBDA(K, I) * ED(K)
 1110  CONTINUE
       DO 1120 I = 1, 6
         F(6 * I1 - 6 + I) = F(6 * I1 - 6 + I) + EDP(I)
         F(6 * I2 - 6 + I) = F(6 * I2 - 6 + I) + EDP(I + 6)
 1120  CONTINUE
      END IF
 1130 CONTINUE
C    ----- Modify for Boundary Conditions -----
C    --- Displacement BC ---
      DO 1140  I = 1, ND
        N = NU(I)
        S(N, 1) = S(N, 1) + CNST
        F(N) = F(N) + CNST * U(I)
 1140 CONTINUE
C    --- Multi-point Constraints ---
      DO 1150  I = 1, NMPC
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
 1150 CONTINUE
C    ----- Equation Solving using Band Solver -----
      CALL BAND(S, F, IMAX, NBW, NQ)
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)      
      WRITE(LOUT,'(1X,2A)') 'Output for Input Data from File ' ,FILE1
      PRINT *, 'Output for Input Data from File ' ,FILE1
      WRITE(LOUT,'(1X,A)')  TITLE
      PRINT *, TITLE
      PRINT *, 'NODE#   X-Displ.    Y-Displ.    Z-Rot.    X-Rot.      Y-
     1Rot.      Z-Rot.'
      WRITE(LOUT,'(A)') 'NODE#   X-Displ.    Y-Displ.    Z-Rot.    X-Rot
     1.      Y-Rot.      Z-Rot.'
      DO 1160 I = 1, NN       
       I1 = 6 * I - 5
       I2 = I1 + 1
       I3 = I1 + 2
       I4 = I1 + 3
       I5 = I1 + 4
       I6 = I1 + 5
       WRITE (*, '(I4,6E12.4)') I,F(I1),F(I2),F(I3),F(I4),F(I5),F(I6)
       WRITE(LOUT, '(I4,6E12.4)') I,F(I1),F(I2),F(I3),F(I4),F(I5),F(I6)
 1160 CONTINUE
C     ----- Reaction Calculation -----
      PRINT *, 'DOF#      Reaction'
      WRITE(LOUT,'(1X,A)')  'DOF#      Reaction'
      DO 1170  I = 1, ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
	    WRITE(*,'(1X,I4,E15.4)') N, R
        WRITE(LOUT,'(1X,I4,E15.4)') N, R
 1170 CONTINUE
C ---- Member End-Actions
      WRITE(LOUT, '(1X,A)') ' Member End-Forces '
      DO 1230 N = 1, NE
       ISTF = 1
       CALL ELSTIF(N,NOC,X,MAT,PM,ARIN,ALMBDA,SE,SEP,EL,ISTF)
       I1 = NOC(N, 1)
       I2 = NOC(N, 2)
       DO 1180 I = 1, 6
         ED(I) = F(6 * I1 - 6 + I)
         ED(I + 6) = F(6 * I2 - 6 + I)
 1180  CONTINUE
       DO 1190 I = 1, 12
         EDP(I) = 0
         DO 1190 K = 1, 12
	   EDP(I) = EDP(I) + ALMBDA(I, K) * ED(K)
 1190 CONTINUE
C  END FORCES DUE TO DISTRIBUTED LOADS
       IF ((ABS(UDL(N, 1)) .GT. 0) .OR. (ABS(UDL(N, 2)) .GT. 0)) THEN
         ED(1) = 0
         ED(4) = 0
         ED(7) = 0
         ED(10) = 0
         ED(2) = -UDL(N, 1) * EL / 2
         ED(8) = ED(2)
         ED(6) = -UDL(N, 1) * EL**2 / 12
         ED(12) = -ED(6)
         ED(3) = -UDL(N, 2) * EL / 2
         ED(9) = ED(3)
         ED(5) = UDL(N, 2) * EL**2 / 12
         ED(11) = -ED(5)
       ELSE
         DO 1200 K = 1, 12
 1200    ED(K) = 0
       END IF
       DO 1210 I = 1, 12
         EF(I) = ED(I)
         DO 1210 K = 1, 12
           EF(I) = EF(I) + SEP(I, K) * EDP(K)
 1210  CONTINUE
       WRITE(LOUT,'(1X,A,I4)') ' Member #', N
       DO 1220 I = 1, 2
	 II = (I - 1) * 6
         WRITE(LOUT, '(1X,6E12.4)') EF(II + 1), EF(II + 2), EF(II + 3), 
     1EF(II + 4), EF(II + 5), EF(II + 6)
 1220  CONTINUE
 1230  CONTINUE
      CLOSE(LOUT)
      END

      SUBROUTINE ELSTIF(N,NOC,X,MAT,PM,ARIN,ALMBDA,SE,SEP,EL,ISTF)
C   ===== SUBROUTINE ELEMENT STIFFNESS =====
      DIMENSION X(150,3),NOC(100,3),MAT(100),PM(10,5),ARIN(100,4)
      DIMENSION SE(12,12),SEP(12,12),DCOS(3,3),ALMBDA(12,12)
C----- Element Stiffness Matrix -----

      I1 = NOC(N, 1)
      I2 = NOC(N, 2)
      I3 = NOC(N, 3)
      M = MAT(N)
      X21 = X(I2, 1) - X(I1, 1)
      Y21 = X(I2, 2) - X(I1, 2)
      Z21 = X(I2, 3) - X(I1, 3)
      EL = SQRT(X21 * X21 + Y21 * Y21 + Z21 * Z21)
      EAL = PM(M, 1) * ARIN(N, 1) / EL
      EIYL = PM(M, 1) * ARIN(N, 2) / EL
      EIZL = PM(M, 1) * ARIN(N, 3) / EL
      GJL = PM(M, 2) * ARIN(N, 4) / EL
      DO 3000 I = 1, 12
      DO 3000 J = 1, 12
       SEP(I, J) = 0
 3000 CONTINUE
      SEP(1, 1) = EAL
      SEP(1, 7) = -EAL
      SEP(7, 7) = EAL
      SEP(4, 4) = GJL
      SEP(4, 10) = -GJL
      SEP(10, 10) = GJL
      SEP(2, 2) = 12 * EIZL / EL**2
      SEP(2, 6) = 6 * EIZL / EL
      SEP(2, 8) = -SEP(2, 2)
      SEP(2, 12) = SEP(2, 6)
      SEP(3, 3) = 12 * EIYL / EL**2
      SEP(3, 5) = -6 * EIYL / EL
      SEP(3, 9) = -SEP(3, 3)
      SEP(3, 11) = SEP(3, 5)
      SEP(5, 5) = 4 * EIYL
      SEP(5, 9) = 6 * EIYL / EL
      SEP(5, 11) = 2 * EIYL
      SEP(6, 6) = 4 * EIZL
      SEP(6, 8) = -6 * EIZL / EL
      SEP(6, 12) = 2 * EIZL
      SEP(8, 8) = 12 * EIZL / EL**2
      SEP(8, 12) = -6 * EIZL / EL
      SEP(9, 9) = 12 * EIYL / EL**2
      SEP(9, 11) = 6 * EIYL / EL
      SEP(11, 11) = 4 * EIYL
      SEP(12, 12) = 4 * EIZL
      DO 3010 I = 1, 12
      DO 3010 J = I, 12
       SEP(J, I) = SEP(I, J)
 3010 CONTINUE
C   CONVERT ELEMENT STIFFNESS MATRIX TO GLOBAL SYSTEM
      DCOS(1, 1) = X21 / EL
      DCOS(1, 2) = Y21 / EL
      DCOS(1, 3) = Z21 / EL
      EIP1 = X(I3, 1) - X(I1, 1)
      EIP2 = X(I3, 2) - X(I1, 2)
      EIP3 = X(I3, 3) - X(I1, 3)
      C1 = DCOS(1, 2) * EIP3 - DCOS(1, 3) * EIP2
      C2 = DCOS(1, 3) * EIP1 - DCOS(1, 1) * EIP3
      C3 = DCOS(1, 1) * EIP2 - DCOS(1, 2) * EIP1
      CC = SQRT(C1 * C1 + C2 * C2 + C3 * C3)
      DCOS(3, 1) = C1 / CC
      DCOS(3, 2) = C2 / CC
      DCOS(3, 3) = C3 / CC
      DCOS(2, 1) = DCOS(3, 2) * DCOS(1, 3) - DCOS(1, 2) * DCOS(3, 3)
      DCOS(2, 2) = DCOS(1, 1) * DCOS(3, 3) - DCOS(3, 1) * DCOS(1, 3)
      DCOS(2, 3) = DCOS(3, 1) * DCOS(1, 2) - DCOS(1, 1) * DCOS(3, 2)
      DO 3020 I = 1, 12
      DO 3020 J = 1, 12
       ALMBDA(I, J) = 0
 3020 CONTINUE
      DO 3030 K = 1, 4
       IK = 3 * (K - 1)
       DO 3030 I = 1, 3
       DO 3030 J = 1, 3
	 ALMBDA(I + IK, J + IK) = DCOS(I, J)
 3030 CONTINUE
      IF (ISTF .EQ. 1) GO TO 3070
      DO 3040 I = 1, 12
      DO 3040 J = 1, 12
       SE(I, J) = 0
       DO 3040 K = 1, 12
	 SE(I, J) = SE(I, J) + SEP(I, K) * ALMBDA(K, J)
 3040 CONTINUE
      DO 3050 I = 1, 12
      DO 3050 J = 1, 12
         SEP(I, J) = SE(I, J)
 3050 CONTINUE
      DO 3060 I = 1, 12
      DO 3060 J = 1, 12
       SE(I, J) = 0
       DO 3060 K = 1, 12
 3060  SE(I, J) = SE(I, J) + ALMBDA(K, I) * SEP(K, J)
 3070 CONTINUE
      RETURN
      END

      SUBROUTINE BAND(A, B, IMAX, NBW, N)
C ----- Equation Solving Using Banded Storage
      DIMENSION A(IMAX,NBW), B(IMAX)
      N1 = N - 1
      PRINT *, '*** FORWARD ELIMINATION ***'
      DO 2100 K = 1, N1
      NK = N - K + 1
      IF (NK .GT. NBW) NK = NBW
      DO 2100 I = 2, NK
      C1 = A(K, I) / A(K, 1)
      I1 = K + I - 1
      DO 2000 J = I, NK
      J1 = J - I + 1
 2000 A(I1, J1) = A(I1, J1) - C1 * A(K, J)
 2100 B(I1) = B(I1) - C1 * B(K)
      PRINT *, '*** BACK SUBSTITUTION ***'
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
