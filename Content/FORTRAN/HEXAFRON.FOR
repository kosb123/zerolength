C *****                    HEXAFRON                    *****
C *****  3-D STRESS ANALYSIS USING HEXAHEDRAL ELEMENT  *****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200,3),NOC(100,8),MAT(100),PM(10,3),NU(50),U(50)
      DIMENSION F(600),D(6,6),B(6,24),DB(6,24),QT(24),MPC(20,2)
      DIMENSION BT(20,3),SE(24,24),XI(3,8),XNI(3,8),STR(6),DT(100)
      DIMENSION S(200,200),ISBL(200),INDX(200)
      CHARACTER*16 FILE1, FILE2
      CHARACTER*81 DUMMY,TITLE      
C --- If IBL > IBMAX increase IBMAX and change dimension of S() ---      
      IBMAX = 200
      PRINT *, '*****        PROGRAM HEXAFRON           *****'
      PRINT *,'*      3-D STRESS ANALYSIS USING  8-NODE     *'
      PRINT *,'*       ISOPARAMETRIC HEXAHEDRAL ELEMENT     *'
      PRINT *,'*             USING FRONTAL SOLVER           *'
      PRINT *,'*       T.R.Chandrupatla and A.DBelegundu    *'
      PRINT *,'**********************************************'
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
      READ (LINP, *) N,(NOC(N,J),J=1,8),MAT(N),DT(N)
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
      NEDF = NEN * NDN
      CALL PREFRN(NN,NE,NEN,NDN,NQ,NMPC,NOC,MPC,IBL)
      IF (IBL .GT. IBMAX) THEN
         PRINT *, 'Increase IBMAX to ', IBL
         PRINT *, 'and change dimension of S(),ISBL(), and INDX()'
         PRINT *, 'in the main program and the subroutines.'
         STOP
      END IF
      DO 1045 I = 1, IBL
 1045 INDX(I) = I
      NFRON = 0
      NTOGO = 0
      NDCNT = 0
C  --- Use same format for writing and reading say '(I4,E16.8)'
      LSCR = 12
      OPEN(UNIT = 12, STATUS = 'SCRATCH', 
     $ACCESS = 'DIRECT', FORM = 'FORMATTED', RECL = 20)
      ICOUNT = 0
C   =====  FRONTAL ASSEMBLY & ELIMINATON ETC.  =====
C   ----- Corner Nodes and Integration Points
      CALL INTEG(XI,XNI)
      MTN1 = 0
      DO 1070 N = 1, NE
         PRINT *, '... Forming Stiffness Matrix of Element ', N
         MTN = MAT(N)
         IF (MTN .NE. MTN1) CALL DMAT(MTN,AL,PM,D)
         MTN1 = MTN
	 CALL ELSTIF(N,SE,QT,XI,XNI,D,B,DB,DT,X,AL,NOC)
	 IF (N .EQ. 1) THEN
            CNST = 0
            DO 1050 I = 1, NEDF
 1050       CNST = CNST + SE(I, I)
            CNST = 1E+11 * CNST
	    CALL MPCFRN(INDX,ISBL,MPC,NMPC,NFRON,S,F,BT,CNST)
         END IF         
C     ----- Account for temperature loads QT()
         DO 1060 I = 1, NEN
           IL = 3 * (I - 1)
           IG = 3 * (ABS(NOC(N, I)) - 1)
           DO 1060 J = 1, 3
              IL = IL + 1
              IG = IG + 1
              F(IG) = F(IG) + QT(IL)
 1060	 CONTINUE
C      ----- Frontal assembly  and Forward Elimination
	 CALL FRONT(N,NOC,NEN,NDN,ND,ICOUNT,INDX,ISBL,S,F,NFRON,
     1NTOGO,NDCNT,SE,NU,CNST,U,LSCR,NEDF)
 1070 CONTINUE
C    ----- Assembly and reduction are complete
C    ----- Now Backsubstitute
      CALL BACSUB(ICOUNT, F, LSCR)
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)
      WRITE(LOUT,'(1X,2A)') 'Output for Input Data from File ', FILE1
      PRINT *, 'Output for Input Data from File ', FILE1
      WRITE(LOUT,'(1X,A)')  TITLE
      PRINT *, TITLE     
      PRINT *, 'Node#        X-Displ        Y-Displ        Z-Displ'
      WRITE(LOUT,'(1X,A)') 'Node#       X-Displ       Y-Displ       Z-Di
     1spl'
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
C   -----  Stress Calculations -----
      MTN1 = 0
      DO 1190 N = 1, NE
         WRITE(LOUT, '(A,I4)') 'von Mises Stress at 8 Integation Pts. in
     1 Elem#	', N
         MTN = MAT(N)
         IF (MTN .NE. MTN1) CALL DMAT(MTN,AL,PM,D)
         MTN1 = MTN
         CAL = AL * DT(N)
         DO 1180 IP = 1, 8
C        --- Von Mises Stress at Integration Points
            CALL DBMAT(N,IP,X,NOC,D,B,DB,DJ,XI,XNI)
C         --- Element Nodal Displacements stored in QT()
            DO 1150 I = 1, 8
               IN = 3 * (ABS(NOC(N, I)) - 1)
               II = 3 * (I - 1)
               DO 1150 J = 1, 3
                  QT(II + J) = F(IN + J)
 1150       CONTINUE
C         --- Stress Calculation STR = DB * Q
            DO 1170 I = 1, 6
               STR(I) = 0
               DO 1160 J = 1, 24
 1160          STR(I) = STR(I) + DB(I, J) * QT(J)          
 1170       STR(I) = STR(I) - CAL * (D(I, 1) + D(I, 2) + D(I, 3))
C         --- Calculation of von Mises Stress at IP
            SIV1 = STR(1) + STR(2) + STR(3)
            SIV2 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1)
            SIV2 = SIV2 - STR(4) ** 2 - STR(5) ** 2 - STR(6) ** 2
            VM = SQRT(SIV1 * SIV1 - 3 * SIV2)
            WRITE (LOUT, '(E15.5)') VM
 1180     CONTINUE
 1190  CONTINUE
       CLOSE(LOUT)
	CLOSE(12)
       PRINT *, 'Results are in the file ',FILE2
      END
      
      SUBROUTINE PREFRN(NN,NE,NEN,NDN,NQ,NMPC,NOC,MPC,IBL)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NOC(100,8),MPC(20,2),IDE(600)
C    ----- Mark Last Appearance of Node / Make it negative in NOC()
C          Last appearance is first appearance for reverse element order
         DO 2010  I = 1, NN
            DO 2000 J = NE, 1, -1
               DO 2000 K = 1, NEN
                  IF (I .EQ. NOC(J, K)) GO TO 2010
 2000    CONTINUE
 2010    NOC(J, K) = -I
C     ===== Block Size Determination
         NQ = NN * NDN
         DO 2020 I = 1, NQ
 2020    IDE(I) = 0
         DO 2030 I = 1, NMPC
         DO 2030  J = 1, 2
 2030    IDE(MPC(I, J)) = 1
         IFRON = 0
	 DO 2035 I = 1, NQ
 2035	 IFRON = IFRON + IDE(I)
         IBL = IFRON
         DO 2060 N = 1, NE
            INEG = 0
            DO 2050 I = 1, NEN
               I1 = NOC(N, I)
               IA = NDN * (ABS(I1) - 1)
               DO 2040 J = 1, NDN
                  IA = IA + 1
                  IF (IDE(IA) .EQ. 0) THEN
                     IFRON = IFRON + 1
                     IDE(IA) = 1
                  END IF
 2040          CONTINUE
               IF (I1 .LT. 0) INEG = INEG + 1
 2050       CONTINUE
            IF (IBL .LT. IFRON) IBL = IFRON
            IFRON = IFRON - NDN * INEG
 2060    CONTINUE
	     RETURN
         END

      SUBROUTINE INTEG(XI,XNI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XI(3,8),XNI(3,8)
C   ------- Integration Points XNI() --------
      C = .57735026919
      XI(1, 1) = -1
      XI(2, 1) = -1
      XI(3, 1) = -1
      XI(1, 2) = 1
      XI(2, 2) = -1
      XI(3, 2) = -1
      XI(1, 3) = 1
      XI(2, 3) = 1
      XI(3, 3) = -1
      XI(1, 4) = -1
      XI(2, 4) = 1
      XI(3, 4) = -1
      XI(1, 5) = -1
      XI(2, 5) = -1
      XI(3, 5) = 1
      XI(1, 6) = 1
      XI(2, 6) = -1
      XI(3, 6) = 1
      XI(1, 7) = 1
      XI(2, 7) = 1
      XI(3, 7) = 1
      XI(1, 8) = -1
      XI(2, 8) = 1
      XI(3, 8) = 1
      DO 2200 I = 1, 8
         XNI(1, I) = C * XI(1, I)
         XNI(2, I) = C * XI(2, I)
         XNI(3, I) = C * XI(3, I)
 2200 CONTINUE
      RETURN
      END
      
      SUBROUTINE DMAT(MTN,AL,PM,D)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PM(10,3),D(6,6)
C   --- D() Matrix relating Stresses to Strains
      E = PM(MTN, 1)
      PNU = PM(MTN, 2)
      AL = PM(MTN, 3)
      C1 = E / ((1 + PNU) * (1 - 2 * PNU))
      C2 = .5 * E / (1 + PNU)
      DO 2300 I = 1, 6
      DO 2300 J = 1, 6
 2300 D(I, J) = 0
      D(1, 1) = C1 * (1 - PNU)
      D(1, 2) = C1 * PNU
      D(1, 3) = D(1, 2)
      D(2, 1) = D(1, 2)
      D(2, 2) = D(1, 1)
      D(2, 3) = D(1, 2)
      D(3, 1) = D(1, 3)
      D(3, 2) = D(2, 3)
      D(3, 3) = D(1, 1)
      D(4, 4) = C2
      D(5, 5) = C2
      D(6, 6) = C2
      RETURN
      END
      
      SUBROUTINE ELSTIF(N,SE,QT,XI,XNI,D,B,DB,DT,X,AL,NOC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200,3),NOC(100,8),D(6,6),B(6,24),DB(6,24)
      DIMENSION SE(24,24),XI(3,8),XNI(3,8),DT(100),QT(24)
C  --------  Element Stiffness  -----
      DO 2400 I = 1, 24
      QT(I) = 0
      DO 2400 J = 1, 24
 2400 SE(I, J) = 0
      DTE = DT(N)
C  --- Weight Factor is ONE
C  --- Loop on Integration Points
      DO 2430 IP = 1, 8
         PRINT *, 'Integration Point = ', IP
         CALL DBMAT(N,IP,X,NOC,D,B,DB,DJ,XI,XNI)
C        --- Element Stiffness Matrix  SE
         DO 2410 I = 1, 24
         DO 2410 J = 1, 24
         DO 2410 K = 1, 6
 2410    SE(I, J) = SE(I, J) + B(K, I) * DB(K, J) * DJ
C    --- Determine Temperature Load QT()
         C = AL * DTE
         DO 2420 I = 1, 24
         DSUM = DB(1, I) + DB(2, I) + DB(3, I)
 2420    QT(I) = QT(I) + C * ABS(DJ) * DSUM / 6
 2430 CONTINUE
      RETURN
      END

      SUBROUTINE DBMAT(N,IP,X,NOC,D,B,DB,DJ,XI,XNI)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(200,3),NOC(100,8),D(6,6),B(6,24),DB(6,24),XI(3,8)
      DIMENSION XNI(3,8),GN(3,8),AJ(3,3),TJ(3,3),H(9,24),G(6,9)    
C    -------  DB()  MATRIX  ------
C    --- Gradient of Shape Functions - The GN() Matrix
      DO 2510 I = 1, 3
         DO 2510 J = 1, 8
            C = 1
            DO 2500 K = 1, 3
               IF (K .NE. I) THEN
                  C = C * (1 + XI(K, J) * XNI(K, IP))
               END IF
 2500       CONTINUE
 2510       GN(I, J) = .125 * XI(I, J) * C
C    --- Formation of Jacobian  TJ
      DO 2520 I = 1, 3
         DO 2520 J = 1, 3
            TJ(I, J) = 0
            DO 2520 K = 1, 8
               KN = ABS(NOC(N, K))
               TJ(I, J) = TJ(I, J) + GN(I, K) * X(KN, J)
 2520 CONTINUE
C    --- Determinant of the JACOBIAN
      DJ1 = TJ(1, 1) * (TJ(2, 2) * TJ(3, 3) - TJ(3, 2) * TJ(2, 3))
      DJ2 = TJ(1, 2) * (TJ(2, 3) * TJ(3, 1) - TJ(3, 3) * TJ(2, 1))
      DJ3 = TJ(1, 3) * (TJ(2, 1) * TJ(3, 2) - TJ(3, 1) * TJ(2, 2))
      DJ = DJ1 + DJ2 + DJ3
C    --- Inverse of the Jacobian AJ()
      AJ(1, 1) = (TJ(2, 2) * TJ(3, 3) - TJ(2, 3) * TJ(3, 2)) / DJ
      AJ(1, 2) = (TJ(3, 2) * TJ(1, 3) - TJ(3, 3) * TJ(1, 2)) / DJ
      AJ(1, 3) = (TJ(1, 2) * TJ(2, 3) - TJ(1, 3) * TJ(2, 2)) / DJ
      AJ(2, 1) = (TJ(2, 3) * TJ(3, 1) - TJ(2, 1) * TJ(3, 3)) / DJ
      AJ(2, 2) = (TJ(1, 1) * TJ(3, 3) - TJ(1, 3) * TJ(3, 1)) / DJ
      AJ(2, 3) = (TJ(1, 3) * TJ(2, 1) - TJ(1, 1) * TJ(2, 3)) / DJ
      AJ(3, 1) = (TJ(2, 1) * TJ(3, 2) - TJ(2, 2) * TJ(3, 1)) / DJ
      AJ(3, 2) = (TJ(1, 2) * TJ(3, 1) - TJ(1, 1) * TJ(3, 2)) / DJ
      AJ(3, 3) = (TJ(1, 1) * TJ(2, 2) - TJ(1, 2) * TJ(2, 1)) / DJ
C     --- H() Matrix relates local derivatives of  u  to local
C         displacements  q
      DO 2530 I = 1, 9
      DO 2530 J = 1, 24
 2530 H(I, J) = 0
      DO 2540 I = 1, 3
      DO 2540 J = 1, 3
        IR = 3 * (I - 1) + J
      DO 2540 K = 1, 8
        IC = 3 * (K - 1) + I
        H(IR, IC) = GN(J, K)
 2540 CONTINUE
C     --- G() Matrix relates strains to local derivatives of  u
      DO 2550 I = 1, 6
      DO 2550 J = 1, 9
 2550 G(I, J) = 0
      G(1, 1) = AJ(1, 1)
      G(1, 2) = AJ(1, 2)
      G(1, 3) = AJ(1, 3)
      G(2, 4) = AJ(2, 1)
      G(2, 5) = AJ(2, 2)
      G(2, 6) = AJ(2, 3)
      G(3, 7) = AJ(3, 1)
      G(3, 8) = AJ(3, 2)
      G(3, 9) = AJ(3, 3)
      G(4, 4) = AJ(3, 1)
      G(4, 5) = AJ(3, 2)
      G(4, 6) = AJ(3, 3)
      G(4, 7) = AJ(2, 1)
      G(4, 8) = AJ(2, 2)
      G(4, 9) = AJ(2, 3)
      G(5, 1) = AJ(3, 1)
      G(5, 2) = AJ(3, 2)
      G(5, 3) = AJ(3, 3)
      G(5, 7) = AJ(1, 1)
      G(5, 8) = AJ(1, 2)
      G(5, 9) = AJ(1, 3)
      G(6, 1) = AJ(2, 1)
      G(6, 2) = AJ(2, 2)
      G(6, 3) = AJ(2, 3)
      G(6, 4) = AJ(1, 1)
      G(6, 5) = AJ(1, 2)
      G(6, 6) = AJ(1, 3)
C    --- B() Matrix relates strains to  q
      DO 2560 I = 1, 6
      DO 2560 J = 1, 24
         B(I, J) = 0
      DO 2560 K = 1, 9
         B(I, J) = B(I, J) + G(I, K) * H(K, J)
 2560 CONTINUE
C     --- DB() Matrix relates stresses to  q
      DO 2570 I = 1, 6
      DO 2570 J = 1, 24
         DB(I, J) = 0
      DO 2570 K = 1, 6
         DB(I, J) = DB(I, J) + D(I, K) * B(K, J)
 2570 CONTINUE
      RETURN
      END

      SUBROUTINE MPCFRN(INDX,ISBL,MPC,NMPC,NFRON,S,F,BT,CNST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(600),MPC(20,2),S(200,200),ISBL(200),INDX(200),
     1BT(20,3)
C     ----- Modifications for Multipoint Constraints by Penalty Method
         DO 2640 I = 1, NMPC
            I1 = MPC(I, 1)
            IFL = 0
            DO 2600 J = 1, NFRON
               J1 = INDX(J)
               IF (I1 .EQ. ISBL(J1)) THEN
                  IFL = 1
                  GO TO 2610
               END IF
 2600       CONTINUE
 2610       IF (IFL .EQ. 0) THEN
               NFRON = NFRON + 1
               J1 = INDX(NFRON)
               ISBL(J1) = I1
            END IF
            I2 = MPC(I, 2)
            IFL = 0
            DO 2620 K = 1, NFRON
               K1 = INDX(K)
               IF (K1 .EQ. ISBL(K1)) THEN
                  IFL = 1
                  GO TO 2630
               END IF
 2620	    CONTINUE
 2630       IF (IFL .EQ. 0) THEN
               NFRON = NFRON + 1
               K1 = INDX(NFRON)
               ISBL(K1) = I2
            END IF
C      ----- Stiffness Modification
            S(J1, J1) = S(J1, J1) + CNST * BT(I, 1) ** 2
            S(K1, K1) = S(K1, K1) + CNST * BT(I, 2) ** 2
            S(J1, K1) = S(J1, K1) + CNST * BT(I, 1) * BT(I, 2)
            S(K1, J1) = S(J1, K1)
C      ----- Force Modification
            F(I1) = F(I1) + CNST * BT(I, 3) * BT(I, 1)
            F(I2) = F(I2) + CNST * BT(I, 3) * BT(I, 2)
 2640    CONTINUE
         RETURN
         END
                  
      SUBROUTINE FRONT(N,NOC,NEN,NDN,ND,ICOUNT,INDX,ISBL,S,F,NFRON,
     1NTOGO,NDCNT,SE,NU,CNST,U,LSCR,NEDF)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION NOC(100,8),NU(50),U(50),F(600),SE(24,24),IEBL(24)
      DIMENSION S(200,200),ISBL(200),INDX(200)      
C     ----- Frontal Method Assembly and Elimination
C     ----------------  Assembly of Element N  --------------------
         DO 3010 I = 1, NEN
            I1 = NOC(N, I)
            IA = ABS(I1)
	    IS1 = 1
	    IF (I1 .LT. 0) IS1 = -1
            IDF = NDN * (IA - 1)
            IE1 = NDN * (I - 1)
            DO 3010 J = 1, NDN
               IDF = IDF + 1
               IE1 = IE1 + 1
               IFL = 0
               IF (NFRON .GT. NTOGO) THEN
                  DO 3000 II = NTOGO + 1, NFRON
                     IX = INDX(II)
                     IF (IDF .EQ. ISBL(IX)) THEN
                        IFL = 1
			GO TO 3005
                     END IF
 3000		  CONTINUE
               END IF
 3005	    CONTINUE
               IF (IFL .EQ. 0) THEN
                  NFRON = NFRON + 1
                  II = NFRON
                  IX = INDX(II)
               END IF
               ISBL(IX) = IDF
               IEBL(IE1) = IX
               IF (IS1 .EQ. -1) THEN
                  NTOGO = NTOGO + 1
                  ITEMP = INDX(NTOGO)
                  INDX(NTOGO) = INDX(II)
                  INDX(II) = ITEMP
               END IF
 3010       CONTINUE
         DO 3020 I = 1, NEDF
            I1 = IEBL(I)
            DO 3020 J = 1, NEDF
               J1 = IEBL(J)
               S(I1, J1) = S(I1, J1) + SE(I, J)
 3020    CONTINUE
C  ------------------------------------------------------------------
      IF (NDCNT .LT. ND) THEN
C  -----  Modification for displacement BCs / Penalty Approach  -----
	 DO 3035 I = 1, NTOGO
            I1 = INDX(I)
            IG = ISBL(I1)
               DO 3030 J = 1, ND
                  IF (IG .EQ. NU(J)) THEN
                     S(I1, I1) = S(I1, I1) + CNST
                     F(IG) = F(IG) + CNST * U(J)
C   --------------  Counter for check
                     NDCNT = NDCNT + 1
		     GO TO 3035
                  END IF
 3030	  CONTINUE
 3035  CONTINUE
      END IF
C    ------------   Elimination of completed variables   -----------
         NTG1 = NTOGO
         DO 3070 II = 1, NTG1
            IPV = INDX(1)
            IPG = ISBL(IPV)
            PIVOT = S(IPV, IPV)
C        -----  Write separator "0" and PIVOT value to disk  -----
            IBA = 0
	      ICOUNT = ICOUNT + 1
		  WRITE(LSCR,'(I4,E16.8)',REC=ICOUNT) IBA, PIVOT
            S(IPV, IPV) = 0
            DO 3050 I = 2, NFRON
               I1 = INDX(I)
               IG = ISBL(I1)
               IF (S(I1, IPV) .NE. 0) THEN
                   C = S(I1, IPV) / PIVOT
                   S(I1, IPV) = 0
                   DO 3040 J = 2, NFRON
                      J1 = INDX(J)
                      IF (S(IPV, J1) .NE. 0) THEN
                         S(I1, J1) = S(I1, J1) - C * S(IPV, J1)
                      END IF
 3040              CONTINUE
                   F(IG) = F(IG) - C * F(IPG)
               END IF
 3050       CONTINUE
            DO 3060 J = 2, NFRON
C       -----  Write Variable# and Reduced Coeff/PIVOT to disk  -----
               J1 = INDX(J)
               IF (S(IPV, J1) .NE. 0) THEN
                  ICOUNT = ICOUNT + 1
                  IBA = ISBL(J1)
		        WRITE(LSCR,'(I4,E16.8)',REC=ICOUNT) IBA, S(IPV, J1) / PIVOT
                  S(IPV, J1) = 0
               END IF
 3060       CONTINUE
            ICOUNT = ICOUNT + 1
C     -----  Write Eliminated Variable# and RHS/PIVOT to disk  -----
	       WRITE(LSCR,'(I4,E16.8)',REC=ICOUNT) IPG, F(IPG) / PIVOT
            F(IPG) = 0
C      ----- (NTOGO) into (1); (NFRON) into (NTOGO)
C      ----- IPV into (NFRON) and reduce front & NTOGO sizes by 1
            IF (NTOGO .GT. 1) THEN
               INDX(1) = INDX(NTOGO)
            END IF
            INDX(NTOGO) = INDX(NFRON)
            INDX(NFRON) = IPV
            NFRON = NFRON - 1
            NTOGO = NTOGO - 1
 3070    CONTINUE 
         RETURN
         END
      SUBROUTINE BACSUB(ICOUNT, F, LSCR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(600)
C     ===== Backsubstitution
 4000    IF (ICOUNT .LE. 0) GO TO 4020
	     READ(LSCR,'(I4,E16.8)', REC = ICOUNT) N1, F(N1)
         ICOUNT = ICOUNT - 1
 4010    CONTINUE
	   READ(LSCR,'(I4,E16.8)', REC = ICOUNT) N2, C
         ICOUNT = ICOUNT - 1
         IF (N2 .EQ. 0) GO TO 4000
         F(N1) = F(N1) - C * F(N2)
         GO TO 4010
 4020    CONTINUE
         RETURN
         END 