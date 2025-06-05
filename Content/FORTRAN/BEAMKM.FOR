C         *******      PROGRAM BEAMKM    ********
C         *****   Shaft Vibration Analysis  *****
      DIMENSION X(51),NOC(50,2),NU(20),U(20),MAT(50),PM(5,3)
      DIMENSION SMI(50),S(102,15),GM(102,15),F(102)
      DIMENSION MPC(10,2),BT(10,3),SE(4,4),EM(4,4),AREA(50)
      CHARACTER*16 FILE1, FILE2
      CHARACTER*81 DUMMY, TITLE
      IMAX = 102
      PRINT *,'***************************************'
      PRINT *,'*            PROGRAM BEAMKM           *'
      PRINT *,'*         Beam Bending Analysis       *'
      PRINT *,'*  T.R.Chandrupatla and A.D.Belegundu *'
      PRINT *,'***************************************'
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
C   --- Material Properties E, Density(Rho)
      NPR = 2
C   ============  READ DATA FROM FILE  ===========
C     ----- Coordinates -----
      READ(LINP,'(A)') DUMMY
      DO 900 I = 1, NN
  900 READ (LINP, *) N, X(N)
C     ----- Connectivity, Material, Mom_Inertia, Area -----
      READ(LINP,'(A)') DUMMY
      DO 1000 I=1,NE
      READ (LINP, *) N,NOC(N,1),NOC(N,2),MAT(N),SMI(N),AREA(N)
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
 1030	READ (LINP, *) N, (PM(N, J), J = 1, NPR)
C     ----- Multi-point Constraints B1*Qi+B2*Qj=B0
      IF (NMPC .GT. 0) THEN
	READ(LINP,'(A)') DUMMY
        DO 1040  I = 1, NMPC
 1040	READ (LINP, *) BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
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
C  ----- Global Stiffness Matrix
      DO 1090 N = 1, NE
	 PRINT *, 'Forming Stiffness and Mass Matrices of Element ', N
	 CALL ELKM(N,X,NOC,MAT,PM,SMI,SE,EM,AREA)
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
C   -----  Additional Springs and Lumped Masses  -----
      PRINT *,	'SPRING SUPPORTS  < dof# = 0 Exits this mode >'
 1130 PRINT *,  'dof#'
      READ *, N
      IF (N .EQ. 0) GO TO 1140
      PRINT *, 'Spring Const'
      READ *, C
      S(N, 1) = S(N, 1) + C
      GOTO 1130
 1140 PRINT *,	'LUMPED MASSES  < dof# = 0 Exits this mode >'
 1145 PRINT *,	'dof#'
      READ *, N
      IF (N .EQ. 0) GO TO 1150
      PRINT *, 'Lumped Mass'
      READ *, C
      GM(N, 1) = GM(N, 1) + C
      GOTO 1145
 1150 CONTINUE      
C --- Print Banded Stiffness and Mass Matrices in Output File
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)      
      WRITE(LOUT,'(1X,2A))') 'Stiffness and Mass for Data in File ', FILE1
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
      PRINT *, 'Run INVITR or JACOBI or GENEIGEN program to get'
      PRINT *, 'Eigenvalues and Eigenvectors'
      END
      
      SUBROUTINE ELKM(N,X,NOC,MAT,PM,SMI,SE,EM,AREA)
      DIMENSION X(51),NOC(50,2),MAT(50),PM(5,3)
      DIMENSION SMI(50),SE(4,4),EM(4,4),AREA(50)      
C  --------  Element Stiffness and Mass Matrices  --------
         N1 = NOC(N, 1)
         N2 = NOC(N, 2)
	 M = MAT(N)
         EL = ABS(X(N1) - X(N2))
C  ----- Element Stiffness -----         
         EIL = PM(M, 1) * SMI(N) / EL**3
         SE(1, 1) = 12 * EIL
         SE(1, 2) = EIL * 6 * EL
         SE(1, 3) = -12 * EIL
         SE(1, 4) = EIL * 6 * EL
            SE(2, 1) = SE(1, 2)
            SE(2, 2) = EIL * 4 * EL * EL
            SE(2, 3) = -EIL * 6 * EL
            SE(2, 4) = EIL * 2 * EL * EL
         SE(3, 1) = SE(1, 3)
         SE(3, 2) = SE(2, 3)
         SE(3, 3) = EIL * 12
         SE(3, 4) = -EIL * 6 * EL
            SE(4, 1) = SE(1, 4)
            SE(4, 2) = SE(2, 4)
            SE(4, 3) = SE(3, 4)
            SE(4, 4) = EIL * 4 * EL * EL
C   --- Element Mass
      RHO = PM(M, 2)
      C1 = RHO * AREA(N) * EL / 420
      EM(1, 1) = 156 * C1
      EM(1, 2) = 22 * EL * C1
      EM(1, 3) = 54 * C1
      EM(1, 4) = -13 * EL * C1
         EM(2, 1) = EM(1, 2)
         EM(2, 2) = 4 * EL * EL * C1
         EM(2, 3) = 13 * EL * C1
         EM(2, 4) = -3 * EL * EL * C1
      EM(3, 1) = EM(1, 3)
      EM(3, 2) = EM(2, 3)
      EM(3, 3) = 156 * C1
      EM(3, 4) = -22 * EL * C1
         EM(4, 1) = EM(1, 4)
         EM(4, 2) = EM(2, 4)
         EM(4, 3) = EM(3, 4)
         EM(4, 4) = 4 * EL * EL * C1                        
      RETURN
      END
