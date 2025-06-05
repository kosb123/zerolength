
C   ------------------   PROGRAM TORSION2   ----------------
      DIMENSION X(200,2),NOC(300,3),MAT(300),PM(10,3),
     $ NU(100),U(100),S(200,50),F(200),BT(2,3)
      CHARACTER*16 FILE1,FILE2,FILE3
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '***************************************'
      PRINT *, '*         PROGRAM  TORSION2           *'
      PRINT *, '*    TORSION WITH 3-NODED TRIANGLES   *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'     
C     IMAX = FIRST DIMENSION OF THE S-MATRIX, IX=1ST DIM. OF X-MATRIX, Etc.
      IMAX = 200                                 
      IX = 200
      INOC = 300
      IPM = 10
      PRINT *,'Input Data File Name <DOS file name>'
      READ '(A)', FILE1
      LINP = 10
      OPEN(UNIT = 10, FILE = FILE1, STATUS = 'UNKNOWN')
      PRINT *,'Output Data File Name <DOS file name>'
      READ '(A)', FILE2
      LOUT = 11
      OPEN(UNIT = 11,FILE = FILE2,STATUS = 'UNKNOWN')
      READ(LINP,'(A)')DUMMY
      READ(LINP,'(A)')TITLE
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) NN, NE, NM, NDIM, NEN, NDN
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) ND, NL, NMPC
C  --- Material Property ThermalConductivity
      NPR = 1
      NMPC = 0
      NM = 1
C  --- ND = NO. OF SPECIFIED STRESS FUNCTION VALUES
C  --- NL = NO. OF GENERALIZED NODAL FORCES "0" HERE
C  --- NPR =1 (SHEAR MODULUS) AND NMPC = 0 FOR THIS PROGRAM
C  --- ELEMENT CHARACTERISTIC NOT USED
C  --- NO. OF MATERIALS = 1 FOR THIS PROGRAM
      PRINT *,'PLOT CHOICE'
      PRINT *,'  1) No Plot Data'
      PRINT *,'  2) Create Data File Containing Stress Function Values'
      PRINT *,'Choose 1 or 2 '
      READ(5,*) IPL
      IF (IPL .LT. 1 .OR. IPL .GT. 2) IPL = 1  
C     --- default is no data
      IF (IPL .GT. 1) THEN
        PRINT *,'Give File Name for Plot Data'
            READ '(A)',FILE3
        LOUT2 = 12
            OPEN(UNIT = 12,FILE = FILE3,STATUS = 'UNKNOWN')
      ENDIF
C     ----- Coordinates
      READ(LINP,'(A)')DUMMY
      DO 20 I = 1, NN
         READ(LINP,*) N ,(X(N,J),J = 1,NDIM)
 20   CONTINUE
C     ----- Connectivity, Material#, Dummy
      READ(LINP,'(A)')DUMMY
      DO 30 I = 1,NE
         READ(LINP,*) N, (NOC(N,J),J = 1,NEN),MAT(N)
 30   CONTINUE
C     ----- Boundary Conditions
      READ(LINP,'(A)') DUMMY
      IF (ND.GT.0)THEN
        DO 40 I = 1,ND
          READ(LINP,*) NU(I), U(I)
 40     CONTINUE
      ENDIF
C     ----- DUMMY READ
      READ(LINP,'(A)') DUMMY
C     ----- Shear Modulus of Material
      READ(LINP,'(A)')DUMMY
      READ (LINP,*)(N,(PM(N,J),J = 1,NPR),I = 1,NM)
C     --- BAND WIDTH CALCULATION ---
      NBW = 0
      DO 60 I = 1,NE
         NMIN = NOC(I, 1)
         NMAX = NOC(I, 1)
         DO 50 J = 2,3
            IF (NMIN .GT. NOC(I, J)) NMIN = NOC(I, J)
            IF (NMAX .LT. NOC(I, J)) NMAX = NOC(I, J)
 50      CONTINUE
         NTMP = NDN * (NMAX - NMIN + 1)
         IF (NBW .LT. NTMP) NBW = NTMP
 60     CONTINUE
      PRINT *, 'THE BAND WIDTH IS ',NBW
C     --- INITIALIZATION OF STIFFNESS MATRIX
      DO 70 I = 1, NN
        F(I) = 0.
      DO 70 J = 1,NBW
        S(I, J) = 0.
 70   CONTINUE
      CLOSE(LINP)
C     --- STIFFNESS MATRIX
      DO 110 I = 1, NE
         I1 = NOC(I, 1)
         I2 = NOC(I, 2)
         I3 = NOC(I, 3)
         X32 = X(I3, 1) - X(I2, 1)
         X13 = X(I1, 1) - X(I3, 1)
         X21 = X(I2, 1) - X(I1, 1)
         Y23 = X(I2, 2) - X(I3, 2)
         Y31 = X(I3, 2) - X(I1, 2)
         Y12 = X(I1, 2) - X(I2, 2)
         DETJ = X13 * Y23 - X32 * Y31
         AREA = .5 * ABS(DETJ)
C        --- LOAD CALCULATION
         C = 2 * AREA / 3
         F(I1) = F(I1) + C
         F(I2) = F(I2) + C
         F(I3) = F(I3) + C
C        --- STIFFNESS FORMATION
         BT(1, 1) = Y23
         BT(1, 2) = Y31
         BT(1, 3) = Y12
         BT(2, 1) = X32
         BT(2, 2) = X13
         BT(2, 3) = X21
         DO 80 II = 1, 3
            DO 80 JJ = 1,2
               BT(JJ, II) = BT(JJ, II) / DETJ
80       CONTINUE
         DO 100 II = 1, 3
            DO 100 JJ = 1, 3
               II1 = NOC(I, II)
               II2 = NOC(I, JJ)
               IF (II1 .LE.II2) THEN
                  SUM = 0.
                  DO 90 J = 1, 2
                     SUM = SUM + BT(J, II) * BT(J, JJ)
 90               CONTINUE
                  IC = II2 - II1 + 1
                  S(II1, IC) = S(II1, IC) + SUM * AREA
               END IF
 100         CONTINUE
 110    CONTINUE
      IF (ND .GT. 0) THEN
C    ----- Decide Penalty Parameter CNST -----
        CNST = S(1, 1)
        DO 120 I = 2, NN
          IF (CNST .LT. S(I, 1)) CNST = S(I, 1)
 120    CONTINUE
        CNST = CNST * 1000000.
C    ----- Modify for Boundary Conditions -----
        DO 130 I = 1, ND
          N = NU(I)
          S(N, 1) = S(N, 1) + CNST
          F(N) = F(N) + CNST * U(I)
 130    CONTINUE
      END IF
C     --- EQUATION SOLVING
      CALL BANSOL(NN,NBW,IMAX,S,F)
      WRITE (LOUT,'(2(1X,A))') 'Output for Input Data File: ',FILE1
      WRITE(LOUT,'(A)') TITLE
      PRINT *, TITLE
      WRITE(LOUT,*) 'NODE#   Stress Function Value'
      PRINT *,'NODE#   Stress Function Value'
      DO 140 I = 1,NN
	 WRITE(*,'(1X,I5,E15.5)') I, F(I)
	 WRITE(LOUT,'(1X,I5,E15.5)') I, F(I)
 140  CONTINUE
      IF (IPL .GT. 1) THEN
            WRITE(LOUT2,*) 'Stress Function Value'
        DO 150 I=1,NN
              WRITE(LOUT2,*)  F(I)
 150    CONTINUE
        PRINT *, 'Stress Function Values in Data file ',FILE3
        PRINT *, 'Run CONTOURA or CONTOURB'
        CLOSE (LOUT2)
      END IF
C    ----- ANGLE OF TWIST PER UNIT LENGTH
      SUM = 0.
      DO 160 I = 1,NE
         I1 = NOC(I, 1)
         I2 = NOC(I, 2)
         I3 = NOC(I, 3)
         X32 = X(I3, 1) - X(I2, 1)
         X13 = X(I1, 1) - X(I3, 1)
         X21 = X(I2, 1) - X(I1, 1)
         Y23 = X(I2, 2) - X(I3, 2)
         Y31 = X(I3, 2) - X(I1, 2)
         Y12 = X(I1, 2) - X(I2, 2)
         DETJ = X13 * Y23 - X32 * Y31
         SUM = SUM + ABS(DETJ) / 3 * (F(I1) + F(I2) + F(I3))        
 160  CONTINUE
      PRINT *, 'TORQUE'
      READ *, TORQUE
      PRINT *, 'SYMMETRY FACTOR (eg. if 1/4 symmetry, then=4.0)= '
      READ *, SFAC
      SMOD = PM(1, 1)
      ALPHA = TORQUE / SMOD / SUM / SFAC
      WRITE(LOUT, '(A,E15.5)') 'TWIST PER UNIT LENGTH = ', ALPHA
      PRINT  *, 'TWIST PER UNIT LENGTH = ', ALPHA
      WRITE(LOUT,*) '-- SHEARING STRESSES TAUYZ, TAUXZ IN EACH ELEMENT'
      WRITE(LOUT,*) 'ELEMENT#    TAUYZ       TAUXZ'
      DO 180 I = 1, NE
        I1 = NOC(I, 1)
        I2 = NOC(I, 2)
        I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1)
        X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2)
        Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        BT(1, 1) = Y23
        BT(1, 2) = Y31
        BT(1, 3) = Y12
        BT(2, 1) = X32
        BT(2, 2) = X13
        BT(2, 3) = X21
        DO 170 II = 1, 3
        DO 170 JJ = 1, 2
 170    BT(JJ, II) = BT(JJ, II) / DETJ
        TAUYZ = -(BT(1, 1)*F(I1) + BT(1, 2)*F(I2) + BT(1, 3)*F(I3))
        TAUXZ = BT(2, 1)*F(I1) + BT(2, 2)*F(I2) + BT(2, 3)*F(I3)
        TAUYZ = TAUYZ * SMOD * ALPHA
        TAUXZ = TAUXZ * SMOD * ALPHA
        WRITE(LOUT,'(I4,2E16.5)') I, TAUYZ, TAUXZ
        WRITE(*,'(1X,I4,2E16.5)') I, TAUYZ, TAUXZ
 180    CONTINUE
      CLOSE(LOUT)
      PRINT *,'Complete results are in file ',FILE2
      END
           
      SUBROUTINE BANSOL(NQ,NBW,IMAX,S,F)
      DIMENSION S(IMAX,1),F(1)
      N = NQ
C     ----- Forward Elimination -----
      DO 1020 K = 1, N-1
         NBK = N - K + 1
         IF ((N - K + 1) .GT. NBW) NBK = NBW
         DO 1010 I = K+1, NBK+K-1
            I1 = I - K + 1
            C = S(K, I1) / S(K, 1)
            DO 1000 J = I, NBK+K-1
               J1 = J - I + 1
               J2 = J - K + 1
               S(I, J1) = S(I, J1) - C * S(K, J2)
 1000       CONTINUE
            F(I) = F(I) - C * F(K)
 1010   CONTINUE
 1020 CONTINUE
C     ----- Back Substitution -----
      F(N) = F(N) / S(N, 1)
      DO 1040 II = 1, N-1
         I = N - II
         NBI = N - I + 1
         IF ((N - I + 1) .GT. NBW) NBI = NBW
         SUM = 0.
         DO 1030 J = 2, NBI
           SUM = SUM + S(I, J) * F(I + J - 1)
 1030    CONTINUE
         F(I) = (F(I) - SUM) / S(I, 1)
 1040 CONTINUE
      RETURN 
      END
