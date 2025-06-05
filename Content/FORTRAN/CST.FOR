C    ------------------------  CST  ---------------------------
      DIMENSION X(200,2),NOC(300,3),MAT(300),PM(10,3),TH(300),
     $ DT(300),NU(100),U(100),F(400),MPC(50,2),BT(50,3),
     $ D(3,3),B(3,6),DB(3,6),SE(6,6),Q(6),STR(3),TL(6),S(400,95)
      CHARACTER*16 FILE1,FILE2,FILE3
      CHARACTER*81 DUMMY,TITLE      
      PRINT *, '***************************************'
      PRINT *, '*          PROGRAM  CST              *'
      PRINT *, '*   2-D  CONSTANT STRAIN TRIANGLE     *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'
C     IMAX = FIRST DIMENSION OF THE S-MATRIX, IX=1ST DIM. OF X-MATRIX, Etc.
      IMAX = 400                                  
      IX=200
      INOC=300
      IPM=10
      IBT=50
      IMPC=50
      PRINT *,'  1) Plane Stress Analysis'
      PRINT *,'  2) Plane Strain Analysis'
      PRINT *,'Choose 1 or 2'
      READ (5,*) LC
      IF (LC .LT. 1 .OR. LC .GT. 2)LC = 1   
C     --- default is Plane Stress
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
      NPR = 3   
C     Dimensioned for 3 properties E, Nu, Alpha
      PRINT *,'PLOT CHOICE'
      PRINT *,'  1) No Plot Data'
      PRINT *,'  2) Create Data File for in-plane Shear Stress'
      PRINT *,'  3) Create Data File for Von Mises Stress'
      PRINT *,'Choose 1 or 2 or 3'
      READ(5,*) IPL
      IF (IPL .LT. 1 .OR. IPL .GT. 3)IPL = 1   
C     --- default is no data
      IF (IPL .GT. 1) THEN
        PRINT *,'Give File Name for Plot Data'
	READ '(A)',FILE3
        LOUT2=12
	OPEN(UNIT=12,FILE=FILE3,STATUS='UNKNOWN')
      ENDIF
C     ----- Total dof is  NQ
      NQ = NDN * NN 
      CALL GETDAT(NQ,NN,NDIM,NE,NEN,ND,NL,NPR,NM,NMPC,IX,INOC,IPM,
     $ IBT,IMPC,LINP,X,NOC,MAT,TH,DT,NU,U,F,PM,BT,MPC)
C'---- - Bandwidth NBW from Connectivity NOC() and MPC
      NBW = 0
      DO 5 I = 1,NE
         NMIN = NOC(I, 1)
         NMAX = NOC(I, 1)
         DO 3 J = 2,3
            IF (NMIN .GT. NOC(I, J)) NMIN = NOC(I, J)
            IF (NMAX .LT. NOC(I, J)) NMAX = NOC(I, J)
3        CONTINUE
         NTMP = NDN * (NMAX - NMIN + 1)
         IF (NBW .LT. NTMP) NBW = NTMP
5     CONTINUE
      DO 7 I = 1,NMPC
         NABS = IABS(MPC(I, 1) - MPC(I, 2)) + 1
         IF (NBW .LT. NABS) NBW = NABS
7     CONTINUE
      PRINT *,'The Bandwidth is', NBW
      DO 101 I=1,NQ
      DO 101 J=1,NBW
101   S(I,J)=0.
C'---- - Global Stiffness Matrix
      DO 25 N = 1,NE
         PRINT *,'Forming Stiffness Matrix of Element', N
         CALL DBMAT(N,IPM,INOC,IX,LC,DJ,MAT,PM,D,NOC,X,B,DB)
C     --- Element Stiffness
         DO 11 I = 1,6
            DO 9 J = 1,6
               C = 0
               DO 8 K = 1,3
                  C = C + .5 * ABS(DJ) * B(K, I) * DB(K, J) * TH(N)
8              CONTINUE
               SE(I, J) = C
9           CONTINUE
11       CONTINUE
C     --- Temperature Load Vector
         M=MAT(N)
         PNU = PM(M, 2)
         AL = PM(M, 3)
         C = AL * DT(N)
         IF (LC .EQ. 2) C = C * (1 + PNU)
         DO 15 I = 1,6
            TL(I) = .5 * C * TH(N) * ABS(DJ) * (DB(1, I) + DB(2, I))
15       CONTINUE
         PRINT *,'.... Placing in Global Locations'
         DO 23 II = 1, NEN
            NRT = NDN * (NOC(N, II) - 1)
            DO 21 IT = 1, NDN
               NR = NRT + IT
               I = NDN * (II - 1) + IT
               DO 19 JJ = 1, NEN
                  NCT = NDN * (NOC(N, JJ) - 1)
                  DO 17 JT = 1, NDN
                     J = NDN * (JJ - 1) + JT
                     NC = NCT + JT - NR + 1
                     IF (NC .GT. 0) THEN
                        S(NR, NC) = S(NR, NC) + SE(I, J)
                     END IF
17                CONTINUE
19              CONTINUE
               F(NR) = F(NR) + TL(I)
21          CONTINUE
23        CONTINUE
25    CONTINUE
C'---- - Decide Penalty Parameter CNST -----
      CNST = 0
      DO 27 I = 1, NQ
         IF (CNST .LT. S(I, 1)) CNST = S(I, 1)
27    CONTINUE
      CNST = CNST * 10000
C'---- - Modify for Boundary Conditions -----
C	 --- Displacement BC ---
      DO 29 I = 1, ND
         N = NU(I)
         S(N, 1) = S(N, 1) + CNST
         F(N) = F(N) + CNST * U(I)
29    CONTINUE
C	  --- Multi-point Constraints ---
         DO 31 I = 1, NMPC
            I1 = MPC(I, 1)
            I2 = MPC(I, 2)
            S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1)
            S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2)
            IR = I1
            IF (IR .GT. I2) IR = I2
            IC = ABS(I2 - I1) + 1
            S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2)
            F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
            F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
31       CONTINUE
C   ----- Equation Solving
      CALL BANSOL(NQ,NBW,IMAX,S,F)
      WRITE (LOUT,'(1X,2A)') 'Output for Input Data File :',FILE1
      WRITE(LOUT,'(1X,A)') TITLE
      WRITE(*,'(1X,A)') TITLE
      IF (LC .EQ. 1) WRITE(LOUT,*)'Plane Stress Analysis'
      IF (LC .EQ. 2) WRITE(LOUT,*)'PRINT #2, "Plane Strain Analysis'
      PRINT *, 'Node#      X-Displ        Y-Displ'
      WRITE(LOUT,'(A)') 'Node#      X-Displ        Y-Displ'
      WRITE(LOUT,'(1X,I4,2E15.4)') (I,F(2*I-1),F(2*I),I=1,NN)
      WRITE(*,'(1X,I4,2E15.4)') (I,F(2*I-1),F(2*I),I=1,NN)
C   ----- Reaction Calculation
      PRINT *,'DOF#        Reaction'
      WRITE(LOUT,*) 'DOF#        Reaction'
      DO 35 I = 1, ND
         N = NU(I)
         R = CNST * (U(I) - F(N))
	 WRITE(*,'(1X,I4,E15.4)') N, R
	 WRITE(LOUT,'(1X,I4,E15.4)') N, R
35    CONTINUE
      IF (IPL .GT. 1) THEN
         IF (IPL.EQ. 2)WRITE(LOUT2,*)'Max. in-plane Shear Stress'
         IF (IPL .EQ. 3)WRITE(LOUT2,*)'Von Mises Stress'
      END IF
C     -----  Stress Calculations
      WRITE(LOUT,*)'ELEM#   SX     SY   TXY   S1    S2   ANGLE SX-->S1'
      DO 37 N = 1, NE
         CALL DBMAT(N,IPM,INOC,IX,LC,DJ,MAT,PM,D,NOC,X,B,DB)
         CALL STRESS(N,LC,INOC,IPM,PM,MAT,NOC,Q,STR,D,DB,DT,F)  
C      --- Principal Stress Calculations
         IF (STR(3) .EQ. 0) THEN
            S1 = STR(1)
            S2 = STR(2)
            ANG = 0
            IF (S2 .GT. S1) THEN
               S1 = STR(2)
               S2 = STR(1)
               ANG = 90
            END IF
         ELSE
            C = .5 * (STR(1) + STR(2))
            R = SQRT(.25 * (STR(1) - STR(2))** 2 + (STR(3))** 2)
            S1 = C + R
            S2 = C - R
            IF (C .GT. STR(1)) THEN
               ANG = 57.2957795 * ATAN(STR(3) / (S1 - STR(1)))
               IF (STR(3) .GT. 0) ANG = 90 - ANG
               IF (STR(3) .GT. 0) ANG = -90 - ANG
            ELSE
               ANG = 57.29577951 * ATAN(STR(3) / (STR(1) - S2))
            END IF
         END IF
	 WRITE(LOUT,'(1X,A,I4)') 'ELEM# ', N
	 WRITE(LOUT,*)'    SX            SY            TXY'
	 WRITE(LOUT,'(1X,3E14.5)') STR(1), STR(2), STR(3)
	 WRITE(LOUT,*)'    S1            S2            ANGLE SR-->S1'
	 WRITE(LOUT,'(1X,3E14.5)') S1, S2, ANG
         IF (IPL .EQ. 2) WRITE(LOUT2,*) .5 * (S1 - S2)
         IF (IPL .EQ. 3) THEN
            S3 = 0
            IF (LC .EQ. 2) S3 = PNU * (S1 + S2)
            C = (S1 - S2) ** 2 + (S2 - S3)** 2 + (S3 - S1) ** 2
            WRITE(LOUT2,*) SQRT(.5 * C)
         END IF
37    CONTINUE
      CLOSE (LOUT)
      PRINT *,'Complete results are in file ',FILE2
      IF (IPL .GT. 1) THEN
       CLOSE (LOUT2)
       PRINT *,'Element Stress Data in file ',FILE3
       PRINT *,'Run BESTFIT and then CONTOURA or CONTOURB to plot stress
     1es'
      END IF
      END 
      
      
      SUBROUTINE GETDAT(NQ,NN,NDIM,NE,NEN,ND,NL,NPR,NM,NMPC,IX,INOC,IPM
     $ ,IBT,IMPC,LINP,X,NOC,MAT,TH,DT,NU,U,F,PM,BT,MPC)
      DIMENSION X(IX,1),NOC(INOC,1),MAT(1),TH(1),DT(1),NU(1),U(1),F(1),
     $ PM(IPM,1),BT(IBT,1),MPC(IMPC,1)
C     ===============  READ DATA  ====================
C     ----- Coordinates
      READ(LINP,'(A)')DUMMY
      DO 39 I = 1, NN
         READ(LINP,*) N ,(X(N,J),J=1,NDIM)
39    CONTINUE
C     ----- Connectivity, Material, Thickness, Temp-change
      READ(LINP,'(A)')DUMMY
      DO 41 I = 1,NE
         READ(LINP,*) N, (NOC(N,J),J=1,NEN),MAT(N),TH(N),DT(N)
41    CONTINUE
C     ----- Displacement BC
      READ(LINP,'(A)')DUMMY
      READ(LINP,*)(NU(I), U(I),I=1,ND)
C     ----- Component Loads
      READ(LINP,'(A)')DUMMY
      DO 29 I=1,NQ
        F(I)=0.
29    CONTINUE      
      READ(LINP,*) (N, F(N),I=1,NL)
C     ----- Material Properties
      READ(LINP,'(A)')DUMMY
         READ (LINP,*)(N,(PM(N,J),J=1,NPR),I=1,NM)
      IF (NMPC .GT. 0) THEN
C	 -----	Multi-point Constraints
         READ(LINP,'(A)')DUMMY
         DO 43 I = 1,NMPC
            READ(LINP,*)BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2),BT(I,3)
43       CONTINUE
      END IF
      CLOSE (LINP)
      RETURN      
      END
      
      
      SUBROUTINE DBMAT(N,IPM,INOC,IX,LC,DJ,MAT,PM,D,NOC,X,B,DB)
      DIMENSION MAT(1),PM(IPM,1),D(3,1),NOC(INOC,1),X(IX,1),B(3,1),
     $ DB(3,1)
C     ----- D(), B() and DB() matrices
C     --- First the D-Matrix
      M = MAT(N)
      E = PM(M, 1)
      PNU = PM(M, 2)
      AL = PM(M, 3)
C     --- D() Matrix
      IF (LC .EQ. 1) THEN
C	 --- Plane Stress
         C1 = E / (1 - PNU ** 2)
         C2 = C1 * PNU
      ELSE
C	 --- Plane Strain
         C = E / ((1 + PNU) * (1 - 2 * PNU))
         C1 = C * (1 - PNU)
         C2 = C * PNU
      END IF
      C3 = .5 * E / (1 + PNU)
      D(1, 1) = C1
      D(1, 2) = C2
      D(1, 3) = 0
      D(2, 1) = C2
      D(2, 2) = C1
      D(2, 3) = 0
      D(3, 1) = 0
      D(3, 2) = 0
      D(3, 3) = C3
C     --- Strain-Displacement Matrix B()
      I1 = NOC(N, 1)
      I2 = NOC(N, 2)
      I3 = NOC(N, 3)
      X1 = X(I1, 1)
      Y1 = X(I1, 2)
      X2 = X(I2, 1)
      Y2 = X(I2, 2)
      X3 = X(I3, 1)
      Y3 = X(I3, 2)
      X21 = X2 - X1
      X32 = X3 - X2
      X13 = X1 - X3
      Y12 = Y1 - Y2
      Y23 = Y2 - Y3
      Y31 = Y3 - Y1
      DJ = X13 * Y23 - X32 * Y31   
C     --- Definition of B() Matrix
      B(1, 1) = Y23 / DJ
      B(2, 1) = 0
      B(3, 1) = X32 / DJ
      B(1, 2) = 0
      B(2, 2) = X32 / DJ
      B(3, 2) = Y23 / DJ
      B(1, 3) = Y31 / DJ
      B(2, 3) = 0
      B(3, 3) = X13 / DJ
      B(1, 4) = 0
      B(2, 4) = X13 / DJ
      B(3, 4) = Y31 / DJ
      B(1, 5) = Y12 / DJ
      B(2, 5) = 0
      B(3, 5) = X21 / DJ
      B(1, 6) = 0
      B(2, 6) = X21 / DJ
      B(3, 6) = Y12 / DJ
C     --- DB Matrix DB = D*B
      DO 49 I = 1, 3
         DO 47 J = 1,6
            C = 0
            DO 45 K = 1,3
               C = C + D(I, K) * B(K, J)
45          CONTINUE
            DB(I, J) = C
47       CONTINUE
49    CONTINUE
      RETURN  
      END
      
      
      SUBROUTINE STRESS(N,LC,INOC,IPM,PM,MAT,NOC,Q,STR,D,DB,DT,F)  
      DIMENSION PM(IPM,1),MAT(1),NOC(INOC,1),Q(1),STR(1),D(3,1),DB(3,1),
     $ DT(1),F(1)
C     ----- Stress Evaluation
      M = MAT(N)
      PNU = PM(M, 2)
      AL = PM(M, 3)
      I1 = NOC(N, 1)
      I2 = NOC(N, 2)
      I3 = NOC(N, 3)
      Q(1) = F(2 * I1 - 1)
      Q(2) = F(2 * I1)
      Q(3) = F(2 * I2 - 1)
      Q(4) = F(2 * I2)
      Q(5) = F(2 * I3 - 1)
      Q(6) = F(2 * I3)
      C1 = AL * DT(N)
      IF (LC .EQ.2) C1 = C1 * (1 + PNU)
      DO 53 I = 1,3
         C = 0
         DO 51 K = 1,6
            C = C + DB(I, K) * Q(K)
51       CONTINUE
         STR(I) = C - C1 * (D(I, 1) + D(I, 2))
53    CONTINUE
      RETURN  
      END
      
      
      SUBROUTINE BANSOL(NQ,NBW,IMAX,S,F)
            DIMENSION S(IMAX,1),F(1)
      N = NQ
C     ----- Forward Elimination -----
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
C     ----- Back Substitution -----
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
