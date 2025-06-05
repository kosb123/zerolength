C     **************           CSTKM           **************
C     ----- Stiffness and Mass Matrices for CST Element -----
      DIMENSION X(100,2),NOC(100,3),MAT(100),PM(10,4)
      DIMENSION NU(50),U(50),MPC(20,2),BT(20,3),TH(100)
      DIMENSION SE(6,6),EM(6,6),S(200,50), GM(200,50)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '***************************************'
      PRINT *, '*          PROGRAM  CSTKM             *'
      PRINT *, '*    STIFFNESS AND MASS MATRICES      *'
      PRINT *, '*   2-D  CONSTANT STRAIN TRIANGLE     *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'
      PRINT *, 'PROBLEM TYPE'
      PRINT *, '  1) Plane Stress'
      PRINT *, '  2) Plane Strain'
      PRINT *, '     Choose 1 or 2'
      READ *, LC
C ------- Default is Plane Stress     
      IF ((LC .LT. 1) .OR. (LC .GT. 2)) LC = 1
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
C    --- Material Properties E, Nu, Alpha, Density(Rho)
      NPR = 4
C    ----- Total dof is  NQ
      NQ = NDN * NN
C   ============  READ DATA FROM FILE  ===========
C     ----- Coordinates -----
      READ(LINP,'(A)') DUMMY
      DO 900 I = 1, NN
  900 READ (LINP, *) N, X(N,1), X(N,2)
C     ----- Connectivity, Material, Mom_Inertia, Area -----
      READ(LINP,'(A)') DUMMY
      DO 1000 I=1,NE
      READ (LINP, *) N,NOC(N,1),NOC(N,2),NOC(N,3),MAT(N),TH(N),C
 1000 CONTINUE
C    ----- Specified Displacements -----
      READ(LINP,'(A)') DUMMY
      DO 1010  I = 1, ND
 1010 READ (LINP, *) NU(I), U(I)
C     ----- Component Loads read as dummy
      READ(LINP,'(A)') DUMMY     
      DO 1020  I = 1, NL
 1020 READ (LINP, *) N, C
C     ----- Material Properties -----
      READ(LINP,'(A)') DUMMY     
      DO 1030  I = 1, NM
 1030   READ (LINP, *) N, (PM(N, J), J = 1, NPR)
C     ----- Multi-point Constraints B1*Qi+B2*Qj=B0
      IF (NMPC .GT. 0) THEN
        READ(LINP,'(A)') DUMMY
        DO 1040  I = 1, NMPC
 1040   READ (LINP,*) BT(I,1),MPC(I,1),BT(I,2),MPC(I,2),BT(I,3)
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
      S(I,J) = 0
 1070 GM(I,J) = 0
C  ----- Global Stiffness and Mass Matrices -----
      DO 1090 N = 1, NE
          PRINT *, 'Forming Stiffness and Mass Matrices of Element ', N
          CALL ELKM(LC,N,MAT,PM,X,NOC,TH,SE,EM)
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
                       GM(NR, NC) = GM(NR, NC) + EM(I, J)                      
                    END IF
 1090        CONTINUE
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
 1120 CONTINUE
C   -----  Additional Springs and Lumped Masses  -----
      PRINT *,  'SPRING SUPPORTS  < dof# = 0 Exits this mode >'
 1130 PRINT *,  'dof#'
      READ *, N
      IF (N .EQ. 0) GO TO 1140
      PRINT *, 'Spring Const'
      READ *, C
      S(N, 1) = S(N, 1) + C
      GOTO 1130
 1140 PRINT *,  'LUMPED MASSES  < dof# = 0 Exits this mode >'
      PRINT *,  'dof#'
      READ *, N
      IF (N .EQ. 0) GO TO 1150
      PRINT *, 'Lumped Mass'
      READ *, C
      GM(N, 1) = GM(N, 1) + C
      GOTO 1140
 1150 CONTINUE     
C --- Print Banded Stiffness and Mass Matrices in Output File
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)
      WRITE(LOUT,'(1X,2A)') 'Stiffness and Mass for Data in File ', FILE1
      WRITE(LOUT,'(1X,A)') 'Num. of DOF    Bandwidth'
      WRITE(LOUT,*) NQ, NBW
      WRITE(LOUT,*) 'BANDED STIFFNESS MATRIX'
      WRITE(LOUT,*) ((S(I, J),J=1,NBW),I=1,NQ)
      WRITE(LOUT,*) 'BANDED MASS MATRIX'
      WRITE(LOUT,*) ((GM(I, J),J=1,NBW),I=1,NQ)
      WRITE(LOUT,*) 'STARTING VECTOR FOR INVERSE ITERATION'
      WRITE(LOUT,*) ('1 ',I=1,NQ)
      CLOSE (LOUT)
      PRINT *, 'Output saved in file ', FILE2
      PRINT *, 'Global Stiffness and Mass Matrices are in file ', FILE2
      PRINT *, 'Run INVITR or JACOBI or GENEIGEN program to get '
      PRINT *, 'Eigenvalues andEigenvectors'
      END

      SUBROUTINE ELKM(LC,N,MAT,PM,X,NOC,TH,SE,EM)
      DIMENSION X(100,2),NOC(100,3),MAT(100),PM(10,4),TH(100)
      DIMENSION D(3,3),B(3,6),DB(3,6),SE(6,6),EM(6,6)
C ----- D(), B() and DB() matrices
C --- First the D-Matrix
      M = MAT(N)
      E = PM(M, 1)
      PNU = PM(M, 2)
      AL = PM(M, 3)
C     --- D() Matrix
      IF (LC .EQ. 1) THEN
C      --- Plane Stress
        C1 = E / (1 - PNU**2)
        C2 = C1 * PNU
      ELSE
C     --- Plane Strain
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
      DO 2010 I = 1, 3
        DO 2010 J = 1, 6
           C = 0
           DO 2000 K = 1, 3
              C = C + D(I, K) * B(K, J)
 2000      CONTINUE
           DB(I, J) = C
 2010  CONTINUE
C     --- Element Stiffness SE()
      DO 2030 I = 1, 6
        DO 2030 J = 1, 6
           C = 0
           DO 2020 K = 1, 3
              C = C + .5 * ABS(DJ) * B(K, I) * DB(K, J) * TH(N)
 2020      CONTINUE
           SE(I, J) = C
 2030   CONTINUE
C     --- Element Mass  EM()
      RHO = PM(M, 4)
      CM = RHO * TH(N) * .5 * ABS(DJ) / 12
      DO 2040 I = 1, 6
      DO 2040 J = 1, 6
      EM(I, J) = 0
 2040 CONTINUE
C     --- Non-zero elements of mass matrix are defined
      EM(1, 1) = 2 * CM
      EM(1, 3) = CM
      EM(1, 5) = CM
      EM(2, 2) = 2 * CM
      EM(2, 4) = CM
      EM(2, 6) = CM
      EM(3, 1) = CM
      EM(3, 3) = 2 * CM
      EM(3, 5) = CM
      EM(4, 2) = CM
      EM(4, 4) = 2 * CM
      EM(4, 6) = CM
      EM(5, 1) = CM
      EM(5, 3) = CM
      EM(5, 5) = 2 * CM
      EM(6, 2) = CM
      EM(6, 4) = CM
      EM(6, 6) = 2 * CM
      RETURN
      END
