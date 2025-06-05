C    -----------------  HEAT2D  -----------------------
      DIMENSION X(200,2),NOC(300,3),MAT(300),PM(10,3),
     $ EHS(300),NU(100),U(100),F(400),BT(2,3),S(400,95)
      CHARACTER*16 FILE1,FILE2,FILE3
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '***************************************'
      PRINT *, '*          PROGRAM  HEAT2D            *'
      PRINT *, '*   HEAT 2-D  WITH 3-NODED TRIANGLES  *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'
C     IMAX = FIRST DIMENSION OF THE S-MATRIX, IX=1ST DIM. OF X-MATRIX, Etc.
      IMAX = 400                                  
      IX=200
      INOC=300
      IPM=10
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
C    --- Material property Thermal Conductivity
      NPR = 1
      NMPC = 0
      NDN = 1
      NDIM=2
      NEN=3
C      ELEMENT HEAT SOURCE, EHS(I),I=1,...,NE
C     --- ND = NO. OF SPECIFIED TEMPERATURES
C     --- NL = NO. OF NODAL HEAT SOURCES
C     ---  NPR =1 (THERMAL CONDUCTIVITY) AND NMPC = 0
      PRINT *,'PLOT CHOICE'
      PRINT *,'  1) No Plot Data'
      PRINT *,'  2) Create Data File for Temperatures'
      PRINT *,'Choose 1 or 2 '
      READ(5,*) IPL
      IF (IPL .LT. 1 .OR. IPL .GT. 2)IPL = 1   
C     --- default is no data
      IF (IPL .GT. 1) THEN
        PRINT *,'Give File Name for Plot Data'
	READ '(A)',FILE3
        LOUT2=12
	OPEN(UNIT=12,FILE=FILE3,STATUS='UNKNOWN')
      ENDIF
C     ----- Coordinates
      READ(LINP,'(A)')DUMMY
      DO 39 I = 1, NN
         READ(LINP,*) N ,(X(N,J),J=1,NDIM)
39    CONTINUE
C     ----- Connectivity, Material, Element Heat Source
      READ(LINP,'(A)')DUMMY
      DO 41 I = 1,NE
         READ(LINP,*) N, (NOC(N,J),J=1,NEN),MAT(N),EHS(N)
41    CONTINUE
C     ----- Temperature BC
      READ(LINP,'(A)')DUMMY
      IF (ND.GT.0)THEN
        DO 50 I=1,ND
          READ(LINP,*)NU(I), U(I)
50      CONTINUE
      ENDIF
C     ----- Nodal Heat Sources
      DO 42 I=1,NN
        F(I)=0.
42    CONTINUE
      READ(LINP,'(A)')DUMMY      
      IF (NL.GT.0)THEN
        DO 51 I=1,NL
          READ(LINP,*) N, F(N)
51      CONTINUE
      ENDIF
C     ----- Thermal Conductivity
      READ(LINP,'(A)')DUMMY
      READ (LINP,*)(N,(PM(N,J),J=1,NPR),I=1,NM)
C     --- BAND WIDTH CALCULATION ---
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
      WRITE(6,*)'THE BAND WIDTH IS',NBW
C     --- INITIALIZATION OF CONDUCTIVITY MATRIX
      DO 8 I = 1, NN
      DO 8 J = 1,NBW
        S(I, J) = 0.
8     CONTINUE
      READ(LINP,'(A)')DUMMY
      READ (LINP,*) NHF
      IF (NHF .GT. 0) THEN
        DO 11 I=1,NHF
           READ(LINP,*) N1, N2, V
           ELEN = SQRT((X(N1, 1) - X(N2, 1)) ** 2 + (X(N1, 2) - 
     $            X(N2, 2)) ** 2)
           F(N1) = F(N1) - ELEN * V / 2
           F(N2) = F(N2) - ELEN * V / 2
11    CONTINUE
      END IF
      READ(LINP,'(A)')DUMMY
      READ(LINP,*)NCONV
      IF (NCONV .GT. 0) THEN
        DO 13 I = 1,NCONV
           READ(LINP,*)N1, N2, H, TINF
           ELEN = SQRT((X(N1, 1) - X(N2, 1))**2 + 
     $           (X(N1, 2) - X(N2, 2)) ** 2)
           F(N1) = F(N1) + ELEN * H * TINF / 2
           F(N2) = F(N2) + ELEN * H * TINF / 2
           S(N1, 1) = S(N1, 1) + H * ELEN / 3
           S(N2, 1) = S(N2, 1) + H * ELEN / 3
           IF (N1 .GE. N2) THEN
              N3 = N1
              N1 = N2
              N2 = N3
           END IF
           S(N1, N2 - N1 + 1) = S(N1, N2 - N1 + 1) + H * ELEN / 6
13    CONTINUE
      END IF
      CLOSE(LINP)
C     --- CONDUCTIVITY MATRIX
      DO 15 I = 1, NE
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
C         --- ELEMENT HEAT SOURCES
         IF (EHS(I) .NE.0.) THEN
            C = EHS(I) * AREA / 3
            F(I1) = F(I1) + C
            F(I2) = F(I2) + C
            F(I3) = F(I3) + C
         END IF
         BT(1, 1) = Y23
         BT(1, 2) = Y31
         BT(1, 3) = Y12
         BT(2, 1) = X32
         BT(2, 2) = X13
         BT(2, 3) = X21
         DO 19 II = 1, 3
            DO 17 JJ = 1,2
               BT(JJ, II) = BT(JJ, II) / DETJ
17          CONTINUE
19       CONTINUE
         DO 21 II = 1,3
            DO 23 JJ = 1,3
               II1 = NOC(I, II)
               II2 = NOC(I, JJ)
               IF (II1 .LE.II2) THEN
                  SUM = 0.
                  DO 25 J = 1,2
                     SUM = SUM + BT(J, II) * BT(J, JJ)
25                CONTINUE
                  IC = II2 - II1 + 1
                  S(II1, IC) = S(II1, IC) + SUM * AREA * PM(MAT(I), 1)
               END IF
23          CONTINUE
21       CONTINUE
15    CONTINUE
      IF (ND .GT. 0) THEN
C'---- - Decide Penalty Parameter CNST -----
        CNST = 0.
        DO 27 I = 1, NN
          IF (CNST .LT. S(I, 1)) CNST = S(I, 1)
27      CONTINUE
        CNST = CNST * 10000.
C'---- - Modify for Boundary Conditions -----
C        --- Temperature BC ---
        DO 29 I = 1, ND
          N = NU(I)
          S(N, 1) = S(N, 1) + CNST
          F(N) = F(N) + CNST * U(I)
29      CONTINUE
      END IF
C     --- EQUATION SOLVING
      CALL BANSOL(NN,NBW,IMAX,S,F)
      WRITE (LOUT,'(1X,2A)') 'Output for Input Data File: ',FILE1
      WRITE(LOUT,'(1X,A)') TITLE
      PRINT *, TITLE
      WRITE(LOUT,*)'NODE#   Temperature'
      PRINT *,'NODE#   Temperature'
      DO 33 I = 1,NN
	 WRITE(*,'(I5,E15.5)') I, F(I)
	 WRITE(LOUT,'(I5,E15.5)') I, F(I)
33    CONTINUE
      IF (IPL .GT. 1) THEN
	WRITE(LOUT2,*)'(Node Temp.) for Data in File ',FILE1
        DO 60 I=1,NN
	  WRITE(LOUT2,*) I,F(I)
60      CONTINUE
      END IF
      CLOSE (LOUT2)
      WRITE(LOUT,*)' -- CONDUCTION HEAT FLOW PER UNIT AREA IN EACH ELEME
     $NT -- '
      WRITE(LOUT,*)'ELEMENT#   QX= -K*DT/DX    QY= -K*DT/DY '
      DO 45 I = 1,NE
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
         DO 35 II = 1,3
            DO 37 JJ = 1,2
               BT(JJ, II) = BT(JJ, II) / DETJ
37          CONTINUE
35       CONTINUE
         QX = BT(1, 1) * F(I1) + BT(1, 2) * F(I2) + BT(1, 3) * F(I3)
         QX = -QX * PM(MAT(I), 1)
         QY = BT(2, 1) * F(I1) + BT(2, 2) * F(I2) + BT(2, 3) * F(I3)
         QY = -QY * PM(MAT(I), 1)
	 WRITE(LOUT, '(I5,2E15.5)') I, QX, QY
45    CONTINUE 
      CLOSE(LOUT)
      PRINT *,'Complete results are in file ',FILE2
      IF (IPL .GT. 1) THEN
       CLOSE (LOUT2)
       PRINT *,'Element Stress Data in file ',FILE3
       PRINT *,'Run CONTOURA or CONTOURB to plot ISOTHERMS'
      END IF
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
