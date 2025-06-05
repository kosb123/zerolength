C     *****		      TETRA3D             *****
C     ***  Three Dimensional Stress Analysis  ***
C     ***        Tetrahedral Elements         ***
      DIMENSION X(100,3),NOC(100,4),MAT(100),PM(10,3),NU(50),U(50)
      DIMENSION S(300,50),F(300),D(6,6),B(6,12),DB(6,12),QT(12)
      DIMENSION SE(12,12),STR(6),MPC(20,2),BT(20,3),DT(100)
      CHARACTER*16 FILE1, FILE2
      CHARACTER*81 DUMMY,TITLE      
      IMAX = 300
      PRINT *, '*******      PROGRAM TETRA3D        *******'
      PRINT *, '*    Three Dimensional Stress Analysis    *'
      PRINT *, '*          Tetrahedral Elements           *'
      PRINT *, '*   T.R.Chandrupatla and A.D.Belegundu    *'
      PRINT *, '*******************************************'
      PRINT *, 'File Name for Input Data '
      READ '(A)', FILE1
      LINP=10
      OPEN (UNIT = 10, FILE = FILE1, STATUS = 'OLD')
      PRINT *, 'Give Name of Output File '
      READ '(A)', FILE2
      READ(LINP,'(A)') DUMMY
      READ(LINP,'(A)') TITLE
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) NN, NE, NM, NDIM, NEN, NDN
      READ(LINP,'(A)') DUMMY      
      READ(LINP,*) ND, NL, NMPC
C    --- Material properties E, Nu, Alpha
      NPR = 3 
C     ----- Total dof is  NQ
      NQ = NDN * NN
C   ============  READ DATA FROM FILE  ===========
C     ----- Coordinates -----
      READ(LINP,'(A)') DUMMY
      DO 900 I = 1, NN
  900 READ (LINP, *) N, X(N,1), X(N,2), X(N,3)
C     ----- Connectivity, Material, Temp-change -----
      READ(LINP,'(A)') DUMMY
      DO 1000 I=1,NE
      READ (LINP, *) N,NOC(N,1),NOC(N,2),NOC(N,3),NOC(N,4),MAT(N),DT(N)
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
      DO 1050 I = 1, NE
        NMIN = NOC(I, 1)
        NMAX = NOC(I, 1)
        DO 1045 J = 2, NEN
           IF (NMIN .GT. NOC(I, J)) NMIN = NOC(I, J)
           IF (NMAX .LT. NOC(I, J)) NMAX = NOC(I, J)
 1045   CONTINUE
        NTMP = NDN * (NMAX - NMIN + 1)
        IF (NBW .LT. NTMP) NBW = NTMP
 1050 CONTINUE
      DO 1060  I = 1, NMPC
        NABS = ABS(MPC(I, 1) - MPC(I, 2)) + 1
        IF (NBW .LT. NABS) NBW = NABS
 1060 CONTINUE
      PRINT *, 'Bandwidth = ', NBW
C    ----- First initialize stiffness and Mass matrices
      DO 1070 I = 1, NQ
      DO 1070 J = 1, NBW
 1070 S(I,J) = 0      
C  ----- Global Stiffness Matrix -----
      DO 1090 N = 1, NE
	     PRINT *, 'Forming Stiffness Matrix of Element ', N
	 CALL DBMAT(N,B,D,DB,MAT,PM,X,NOC,DJ,AL,PNU)
C     --- Element Stiffness
         DO 1072 I = 1, 12
         DO 1072 J = 1, 12
         SE(I, J) = 0
         DO 1072 K = 1, 6
            SE(I, J) = SE(I, J) + B(K, I) * DB(K, J) * ABS(DJ) / 6
 1072    CONTINUE
C      --- Temperature Load Vector QT()
         C = AL * DT(N)
         DO 1074 I = 1, 12
            DSUM = DB(1, I) + DB(2, I) + DB(3, I)
            QT(I) = C * ABS(DJ) * DSUM / 6
 1074    CONTINUE 	  
	  PRINT *, '.... Placing in Global Locations'
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
	      F(NR) = F(NR) + QT(I)
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
C   ----- Equation Solving
      PRINT *, '.... Solving Equations'
      CALL BAND(S, F, IMAX, NBW, NQ)      
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)
      WRITE(LOUT,'(1X,2A)') 'Output for Input Data from File ', FILE1
      PRINT *, 'Output for Input Data from File ', FILE1
      WRITE(LOUT,'(1X,A)')  TITLE
      PRINT *, TITLE     
      PRINT *, 'Node#  X-Displ        Y-Displ        Z-Displ'
      WRITE(LOUT,'(1X,A)') 'Node#  X-Displ       Y-Displ       Z-Displ'
      WRITE(LOUT,'(1X,I4,3E15.4)') (I,F(3*I-2),F(3*I-1),F(3*I),I=1,NN)
      WRITE(*,'(1X,I4,3E15.4)') (I,F(3*I-2),F(3*I-1),F(3*I),I=1,NN)
C     ----- Reaction Calculation -----
      PRINT *, 'DOF#      Reaction'
      WRITE(LOUT,'(1X,A)')  'DOF#      Reaction'
      DO 1140  I = 1, ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
        WRITE(*,'(1X,I4,E15.4)') N, R
        WRITE(LOUT,'(1X,I4,E15.4)') N, R
 1140 CONTINUE            
C    --- Stress Calculations
      PI = 3.141593
      DO 1150 N = 1, NE
	 CALL DBMAT(N,B,D,DB,MAT,PM,X,NOC,DJ,AL,PNU)
         CALL STRESS(F,AL,NOC,DT,N,D,DB,STR)
C    --- Principal Stress Calculations
         AI1 = STR(1) + STR(2) + STR(3)
         AI21 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1)
         AI22 = STR(4) * STR(4) + STR(5) * STR(5) + STR(6) * STR(6)
         AI2 = AI21 - AI22
         AI31 = STR(1) * STR(2) * STR(3) + 2 * STR(4) * STR(5) * STR(6)
         AI32 = STR(1)*STR(4)**2 + STR(2)*STR(5)**2 + STR(3)*STR(6)**2
         AI3 = AI31 - AI32
         C1 = AI2 - AI1**2 / 3
         C2 = -2 * (AI1 / 3)**3 + AI1 * AI2 / 3 - AI3
         C3 = 2 * SQRT(-C1 / 3)
         TH = -3 * C2 / (C1 * C3)
         TH2 = ABS(1 - TH * TH)
         IF (TH .EQ. 0) TH = PI / 2
         IF (TH .GT. 0) TH = ATAN(SQRT(TH2) / TH)
         IF (TH .LT. 0) TH = PI - ATAN(SQRT(TH2) / TH)
         TH = TH / 3
C      --- Principal Stresses
         P1 = AI1 / 3 + C3 * COS(TH)
         P2 = AI1 / 3 + C3 * COS(TH + 2 * PI / 3)
         P3 = AI1 / 3 + C3 * COS(TH + 4 * PI / 3)
         WRITE(LOUT,*) ' ELEMENT NO. '
         WRITE(LOUT,'(I5)') N
         WRITE(LOUT,*) '      SX             SY             SZ'
         WRITE(LOUT,'(3E15.4)') STR(1), STR(2), STR(3)
         WRITE(LOUT,*) '     TYZ            TXZ            TXY'
         WRITE(LOUT,'(3E15.4)') STR(4), STR(5), STR(6)
         WRITE(LOUT,*) ' Principal Stresses'
         WRITE(LOUT,*) '      P1             P2             P3'
         WRITE(LOUT,'(3E15.4)') P1, P2, P3
 1150  CONTINUE
      CLOSE(LOUT)
      PRINT *, 'Results are in the file ',FILE2
      END

      SUBROUTINE DBMAT(N,B,D,DB,MAT,PM,X,NOC,DJ,AL,PNU)
      DIMENSION X(100,3),NOC(100,4),MAT(100),PM(10,3)
      DIMENSION D(6,6),B(6,12),DB(6,12)      
C   ----- D(), B() and DB() matrices
C    --- First the D-Matrix
      M = MAT(N)
      E = PM(M, 1)
      PNU = PM(M, 2)
      AL = PM(M, 3)
      C4 = E / ((1 + PNU) * (1 - 2 * PNU))
      C1 = C4 * (1 - PNU)
      C2 = C4 * PNU
      C3 = .5 * E / (1 + PNU)
      DO 3000 I = 1, 6
      DO 3000 J = 1, 6
      D(I, J) = 0
 3000 CONTINUE
      D(1, 1) = C1
      D(1, 2) = C2
      D(1, 3) = C2
      D(2, 1) = C2
      D(2, 2) = C1
      D(2, 3) = C2
      D(3, 1) = C2
      D(3, 2) = C2
      D(3, 3) = C1
      D(4, 4) = C3
      D(5, 5) = C3
      D(6, 6) = C3
C   --- Strain-Displacement Matrix B()
      I1 = NOC(N, 1)
      I2 = NOC(N, 2)
      I3 = NOC(N, 3)
      I4 = NOC(N, 4)
      X14 = X(I1, 1) - X(I4, 1)
      X24 = X(I2, 1) - X(I4, 1)
      X34 = X(I3, 1) - X(I4, 1)
      Y14 = X(I1, 2) - X(I4, 2)
      Y24 = X(I2, 2) - X(I4, 2)
      Y34 = X(I3, 2) - X(I4, 2)
      Z14 = X(I1, 3) - X(I4, 3)
      Z24 = X(I2, 3) - X(I4, 3)
      Z34 = X(I3, 3) - X(I4, 3)
      DJ1 = X14 * (Y24 * Z34 - Z24 * Y34)
      DJ2 = Y14 * (Z24 * X34 - X24 * Z34)
      DJ3 = Z14 * (X24 * Y34 - Y24 * X34)
      DJ = DJ1 + DJ2 + DJ3
      A11 = (Y24 * Z34 - Z24 * Y34) / DJ
      A21 = (Z24 * X34 - X24 * Z34) / DJ
      A31 = (X24 * Y34 - Y24 * X34) / DJ
      A12 = (Y34 * Z14 - Z34 * Y14) / DJ
      A22 = (Z34 * X14 - X34 * Z14) / DJ
      A32 = (X34 * Y14 - Y34 * X14) / DJ
      A13 = (Y14 * Z24 - Z14 * Y24) / DJ
      A23 = (Z14 * X24 - X14 * Z24) / DJ
      A33 = (X14 * Y24 - Y14 * X24) / DJ
C    ---  B Matrix
      DO 3010 I = 1, 6
      DO 3010 J = 1, 12
 3010 B(I, J) = 0
      B(1, 1) = A11
      B(1, 4) = A12
      B(1, 7) = A13
      B(1, 10) = -A11 - A12 - A13
      B(2, 2) = A21
      B(2, 5) = A22
      B(2, 8) = A23
      B(2, 11) = -A21 - A22 - A23
      B(3, 3) = A31
      B(3, 6) = A32
      B(3, 9) = A33
      B(3, 12) = -A31 - A32 - A33
      B(4, 2) = A31
      B(4, 3) = A21
      B(4, 5) = A32
      B(4, 6) = A22
      B(4, 8) = A33
      B(4, 9) = A23
      B(4, 11) = B(3, 12)
      B(4, 12) = B(2, 11)
      B(5, 1) = A31
      B(5, 3) = A11
      B(5, 4) = A32
      B(5, 6) = A12
      B(5, 7) = A33
      B(5, 9) = A13
      B(5, 10) = B(3, 12)
      B(5, 12) = B(1, 10)
      B(6, 1) = A21
      B(6, 2) = A11
      B(6, 4) = A22
      B(6, 5) = A12
      B(6, 7) = A23
      B(6, 8) = A13
      B(6, 10) = B(2, 11)
      B(6, 11) = B(1, 10)
C   --- DB Matrix DB = D*B
      DO 3020 I = 1, 6
      DO 3020 J = 1, 12
      DB(I, J) = 0
      DO 3020 K = 1, 6
 3020 DB(I, J) = DB(I, J) + D(I, K) * B(K, J)
      RETURN
      END
      
      SUBROUTINE STRESS(F,AL,NOC,DT,N,D,DB,STR)
      DIMENSION F(300),NOC(100,4),D(6,6),DB(6,12),QT(12),STR(6),DT(100)      
C --- Stress Evaluation (Element Nodal Displacements stored in QT() )
      DO 4000 I = 1, 4
         IN = 3 * (NOC(N, I) - 1)
         II = 3 * (I - 1)
         DO 4000 J = 1, 3
            QT(II + J) = F(IN + J)
 4000 CONTINUE
            C1 = AL * DT(N)
            DO 4020 I = 1, 6
               STR(I) = 0
               DO 4010K = 1, 12
 4010          STR(I) = STR(I) + DB(I, K) * QT(K)
 4020       STR(I) = STR(I) - C1 * (D(I, 1) + D(I, 2) + D(I, 3))
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
 
