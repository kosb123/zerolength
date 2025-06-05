C    -----------------------  FEM1D  ------------------------
      DIMENSION X(50),NOC(25,2),F(50),AREA(25),MAT(25),DT(25),
     $  PM(10,2),NU(20),U(20),MPC(15,2),BT(15,3)
      DIMENSION S(100,55)
C     IMAX = FIRST DIMENSION OF THE S-MATRIX
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '***************************************'
      PRINT *, '*          PROGRAM FEM1D              *'
      PRINT *, '*    WITH MULTI-POINT CONSTRAINTS     *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'      
      IMAX=100
      PRINT *,'Input Data File Name <DOS file name>'
      READ '(A)',FILE1
      LINP=10
      OPEN(UNIT=10,FILE=FILE1,STATUS='UNKNOWN')
      PRINT *,'Output Data File Name <DOS file name>'
      READ '(A)',FILE2
      LOUT=11
      OPEN(UNIT=11,FILE=FILE2,STATUS='UNKNOWN')
      READ(LINP,'(A)')DUMMY
      READ(LINP,'(A)')TITLE
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) NN, NE, NM, NDIM, NEN, NDN
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) ND, NL, NMPC
C     Material Properties E, Alpha
      NPR = 2   
C      ----- Coordinates -----
      READ(LINP,'(A)')DUMMY   
      READ(LINP,*) (N,X(N),I=1,NN)
C      ----- Connectivity -----
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) (N,NOC(N,1),NOC(N,2),MAT(N),AREA(N),DT(N),I=1,NE)
C      ----- Specified Displacements -----
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) (NU(I),U(I),I=1,ND)
C      ----- Component Loads -----
      READ(LINP,'(A)')DUMMY
      DO 101 I=1,NN
101   F(I)=0.
      READ(LINP,*) (N,F(N),I=1,NL)
C      ----- Material Properties -----
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) (N,(PM(N,J),J=1,NPR),I=1,NM)
C      ----- Multi-point Constraints B1*Qi+B2*Qj=B0
      IF (NMPC .GT. 0) THEN
         READ(LINP,'(A)')DUMMY
         DO 9 I=1,NMPC
          READ(LINP,*)BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
9        CONTINUE
      ENDIF
      CLOSE (LINP)
C      ----- Bandwidth Evaluation -----
      NBW = 0
      DO 11 N=1,NE
         NABS = IABS(NOC(N, 1) - NOC(N, 2)) + 1
         IF (NBW .LT. NABS) NBW = NABS
11    CONTINUE
      DO 13 I=1,NMPC
         NABS = IABS(MPC(I, 1) - MPC(I, 2)) + 1
         IF (NBW .LT. NABS) NBW = NABS
13    CONTINUE
      DO 102 I=1,NN
      DO 102 J=1,NBW
102   S(I,J)=0.      
C      ----- Stiffness Matrix -----
      DO 25 N=1,NE
         N1 = NOC(N, 1)
         N2 = NOC(N, 2)
         N3 = MAT(N)
         X21 = X(N2) - X(N1)
         EL = ABS(X21)
         EAL = PM(N3, 1) * AREA(N) / EL
         IF (NPR .GT. 1) C = PM(N3, 2)
         TL = PM(N3, 1) * C * DT(N) * AREA(N) * EL / X21
C	  ----- Temperature Loads -----
         F(N1) = F(N1) - TL
         F(N2) = F(N2) + TL
C	  ----- Element Stiffness in Global Locations -----
         S(N1, 1) = S(N1, 1) + EAL
         S(N2, 1) = S(N2, 1) + EAL
         IR = N1
         IF (IR .GT. N2) IR = N2
         IC = IABS(N2 - N1) + 1
         S(IR, IC) = S(IR, IC) - EAL
25    CONTINUE
C      ----- Decide Penalty Parameter CNST -----
      CNST = 0
      DO 27 I=1,NN
        IF (CNST .LT. S(I, 1)) CNST = S(I, 1)
27    CONTINUE
      CNST = CNST * 10000
C      ----- Modify for Boundary Conditions -----
C	  --- Displacement BC ---
         DO 29 I=1,ND
            N = NU(I)
            S(N, 1) = S(N, 1) + CNST
            F(N) = F(N) + CNST * U(I)
29       CONTINUE
C	  --- Multi-point Constraints ---
         DO 31 I=1,NMPC
            I1 = MPC(I, 1)
            I2 = MPC(I, 2)
            S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1)
            S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2)
            IR = I1
            IF (IR .GT. I2) IR = I2
            IC = IABS(I2 - I1) + 1
            S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2)
            F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
            F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
31       CONTINUE
C      ----- Equation Solving using Band Solver -----
      CALL BANSOL(NN,NBW,IMAX,S,F)
      WRITE(*,'(1X,A)') TITLE
      WRITE(LOUT,'(1X,A)') TITLE
      PRINT *,'NODE#   DISPLACEMENT'
      WRITE(LOUT,*)'NODE#   DISPLACEMENT'
      DO 33 I=1,NN
	 WRITE(*,'(1X,I5,E15.5)') I, F(I)
	 WRITE(LOUT,'(1X,I5,E15.5)') I, F(I)
33    CONTINUE
C      ----- Stress Calculation -----
      PRINT *, 'ELEM#    STRESS'
      WRITE(LOUT,*) 'ELEM#    STRESS'
      DO 35 N=1,NE
         N1 = NOC(N, 1)
	 N2 = NOC(N, 2)
         N3 = MAT(N)
         EPS = (F(N2) - F(N1)) / (X(N2) - X(N1))
         IF (NPR .GT. 1) C = PM(N3, 2)
         STRESS = PM(N3, 1) * (EPS - C * DT(N))
	 WRITE(*,'(1X,I5,E15.5)') N, STRESS
	 WRITE(LOUT,'(1X,I5,E15.5)') N, STRESS
35    CONTINUE
C      ----- Reaction Calculation -----
      PRINT *, 'NODE#    REACTION'
      WRITE(LOUT,*) 'NODE#    REACTION'
      DO 37 I=1,ND
         N = NU(I)
         R = CNST * (U(I) - F(N))
	 WRITE(*,'(1X,I5,E15.5)') N, R
	 WRITE(LOUT,'(1X,I5,E15.5)') N, R
37    CONTINUE
      CLOSE(LOUT)
      PRINT *, 'RESULTS ARE IN FILE ', FILE2
      END

      SUBROUTINE BANSOL(NN,NBW,IMAX,S,F)
      DIMENSION S(IMAX,1),F(1)
      N = NN
C      ----- Forward Elimination -----
      DO 39 K=1,N-1
         NBK = N - K + 1
         IF ((N - K + 1) .GT. NBW) NBK = NBW
         DO 43 I=K+1, NBK+K-1
            I1 = I - K + 1
            C = S(K, I1) / S(K, 1)
            DO 41 J=I, NBK+K-1
               J1 = J - I + 1
               J2 = J - K + 1
               S(I, J1) = S(I, J1) - C * S(K, J2)
41          CONTINUE
            F(I) = F(I) - C * F(K)
43        CONTINUE
39    CONTINUE
C      ----- Back Substitution -----
      F(N) = F(N) / S(N, 1)
      DO 47 II=1,N-1
         I = N - II
         NBI = N - I + 1
         IF ((N - I + 1) .GT. NBW) NBI = NBW
         SUM = 0.
         DO 45 J=2,NBI
           SUM = SUM + S(I, J) * F(I + J - 1)
45       CONTINUE
         F(I) = (F(I) - SUM) / S(I, 1)
47    CONTINUE
      RETURN  
      END
