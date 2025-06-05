C     ********           JACOBI           **********
      DIMENSION S(100,100),GM(100,100),EVL(100)
      DIMENSION EVC(100,100),NORD(100)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY
      IMAX = 100
      PRINT *, '*****       PROGRAM JACOBI            *'
      PRINT *, '*     Generalized Jacobi Method       *'
      PRINT *, '*       for symmetric matrices        *'
      PRINT *, '*  T.R.Chandrupatla and A.D.Belegundu *'
      PRINT *, '***************************************'
      PRINT *,'Give Name of Input File '
      READ '(A)', FILE1
      LINP = 10
      OPEN (UNIT = 10, FILE = FILE1, STATUS = 'OLD')     
      PRINT *,'Give Name of Output File '
      READ '(A)', FILE2
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)
C   --- Read in Number of Equations and Bandwidth ---
      READ(LINP,'(A)') DUMMY
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) NQ, NBW
C   NORD( ) is for ascending order of eigenvalues
      DO 100 I = 1, NQ
  100 NORD(I) = I
      TOL = .000001
      PI = 3.14159
      PRINT *, 'DEFAULT TOLERANCE IS  1E-6'
C   --- Banded Stiffness Matrix into S(NQ,NQ) ---
      READ(LINP,'(A)') DUMMY
C  Eigenvector array EVC() is used for reading in the Data
      READ(LINP,*)((EVC(I,J),J=1,NBW),I=1,NQ)
      DO 110 I = 1, NQ
      DO 110 JN = 1, NBW
      J = I + JN - 1
      IF(J .GT. NQ) GO TO 110
      S(I, J) = EVC(I,JN)
      S(J, I) = EVC(I,JN)
  110 CONTINUE
C   --- Read in Mass Matrix into GM(NQ,NQ) ---
      READ(LINP,'(A)') DUMMY
      READ(LINP,*)((EVC(I,J),J=1,NBW),I=1,NQ)
      DO 120 I = 1, NQ
      DO 120 JN = 1, NBW
      J = I + JN - 1
      IF(J .GT. NQ) GO TO 120
      GM(I, J) = EVC(I,JN)
      GM(J, I) = EVC(I,JN)
  120 CONTINUE
      CLOSE(LINP)
C   --- Initialize Eigenvector Matrix ---
      DO 140 I = 1, NQ
      DO 130 J = 1, NQ
  130 EVC(I, J) = 0.
  140 EVC(I, I) = 1.
      C1 = S(1, 1)
      C2 = GM(1, 1)
      DO 150 I = 1, NQ
      IF(C1 .GT. S(I, I)) C1 = S(I, I)
      IF(C2 .LT. GM(I, I)) C2 = GM(I, I)
  150 CONTINUE
      TOLS = TOL * C1
      TOLM = TOL * C2
      PRINT *, 'Max. Number of Sweeps .. enter number >= 50'
      READ *, NSWMAX
C   ----- Generalixed Jacobi's Method -----
      K1 = 1
      I1 = 1
      NSW = 0
  160 NSW = NSW + 1
      IF(NSW .GT. NSWMAX) GO TO 330
      PRINT *, '==== SWEEP NUMBER ==> ', NSW
      NQM1 = NQ - 1
      DO 260 K = K1, NQM1
      DO 260 I = I1, K
      J = NQ - K + I
      IF((ABS(S(I,J)).LE.TOLS).AND.(ABS(GM(I,J)).LE.TOLM)) GO TO 260
      AA = S(I, I) * GM(I, J) - GM(I, I) * S(I, J)
      BB = S(J, J) * GM(I, J) - GM(J, J) * S(I, J)
      CC = S(I, I) * GM(J, J) - GM(I, I) * S(J, J)
      CAB = .25 * CC * CC + AA * BB
      IF(CAB .LT. 0.) GO TO 340
      IF(AA .EQ. 0.) GO TO 170
      IF(BB .EQ. 0.) GO TO 180
      SQC = SQRT(CAB)
      IF(CC .LT. 0) SQC = -SQC
      ALP = (-.5 * CC + SQC) / AA
      BET = -AA * ALP / BB
      GO TO 190
  170 BET = 0.
      ALP = -S(I, J) / S(I, I)
      GO TO 190
  180 ALP = 0.
      BET = -S(I, J) / S(J, J)
C   --- Only Upper Triangular Part is used in Diagonalization ---
  190 IF(I .EQ. 1) GO TO 205
      IM1 = I - 1
      DO 200 N = 1, IM1
      SI = S(N, I)
      SJ = S(N, J)
      EMI = GM(N, I)
      EMJ = GM(N, J)
      S(N, I) = SI + BET * SJ
      S(N, J) = SJ + ALP * SI
      GM(N, I) = EMI + BET * EMJ
      GM(N, J) = EMJ + ALP * EMI
  200 CONTINUE
  205 IF(J .EQ. NQ) GO TO 220
      JP1 = J + 1
      DO 210 N = JP1, NQ
      SI = S(I, N)
      SJ = S(J, N)
      EMI = GM(I, N)
      EMJ = GM(J, N)
      S(I, N) = SI + BET * SJ
      S(J, N) = SJ + ALP * SI
      GM(I, N) = EMI + BET * EMJ
      GM(J, N) = EMJ + ALP * EMI
  210 CONTINUE
  220 IF(J .EQ. I + 1) GO TO 240
      IP1 = I + 1
      JM1 = J - 1
      DO 230 N = IP1, JM1
      SI = S(I, N)
      SJ = S(N, J)
      EMI = GM(I, N)
      EMJ = GM(N, J)
      S(I, N) = SI + BET * SJ
      S(N, J) = SJ + ALP * SI
      GM(I, N) = EMI + BET * EMJ
      GM(N, J) = EMJ + ALP * EMI
  230 CONTINUE
  240 SII = S(I, I)
      SIJ = S(I, J)
      SJJ = S(J, J)
      S(I, J) = 0
      S(I, I) = SII + 2 * BET * SIJ + BET * BET * SJJ
      S(J, J) = SJJ + 2 * ALP * SIJ + ALP * ALP * SII
      EII = GM(I, I)
      EIJ = GM(I, J)
      EJJ = GM(J, J)
      GM(I, J) = 0
      GM(I, I) = EII + 2 * BET * EIJ + BET * BET * EJJ
      GM(J, J) = EJJ + 2 * ALP * EIJ + ALP * ALP * EII
C   --- Eigenvectors ---
      DO 250 N = 1, NQ
      EVI = EVC(N, I)
      EVJ = EVC(N, J)
      EVC(N, I) = EVI + BET * EVJ
      EVC(N, J) = EVJ + ALP * EVI
  250 CONTINUE
  260 CONTINUE
      DO 270 K = 1, NQM1
      DO 270 I = 1, K
      J = NQ - K + I
      IF((ABS(S(I,J)).LE.TOLS).AND.(ABS(GM(I,J)).LE.TOLM)) GO TO 270
      K1 = K
      I1 = I
      GO TO 160
  270 CONTINUE
C   --- Calculation of Eigenvalues ---
      DO 280 I = 1, NQ
      IF(ABS(GM(I, I)) .LT. TOLM) GM(I, I) = TOLM
  280 EVL(I) = S(I, I) / GM(I, I)
C   --- Scaling of Eigenvectors ---
      DO 290 I = 1, NQ
      EM2 = SQRT(ABS(GM(I, I)))
      DO 290 J = 1, NQ
  290 EVC(J, I) = EVC(J, I) / EM2
C   -----   RESULTS   -----
C   Ascending Order of Eigenvalues
      DO 310 I = 1, NQ-1
      II = NORD(I)
      I1 = II
      C1 = EVL(II)
      J1 = I
      DO 300 J = I+1, NQ
      IJ = NORD(J)
      IF(C1 .LE. EVL(IJ)) GO TO 300
      C1 = EVL(IJ)
      I1 = IJ
      J1 = J
  300 CONTINUE
      NORD(I) = I1
      NORD(J1) = II
  310 CONTINUE
      DO 320 I = 1, NQ
      II = NORD(I)
      OMEGA = SQRT(EVL(II))
      FREQ = .5 * OMEGA / PI
      PRINT '(1X,A,I4)' , ' Eigenvalue Number', I
      WRITE(LOUT,'(1X,A,I4)') 'Eigenvalue Number', I
      PRINT '(3(1X,A,E12.5),A)', ' Eigenvalue = ', EVL(II), '    Omega =
     1 ', OMEGA, '    Freq = ', FREQ, ' Hz'
      WRITE(LOUT,'(3(1X,A,E12.5),A)') 'Eigenvalue = ', EVL(II), '    Ome
     1ga= ', OMEGA, '	 Freq = ', FREQ, ' Hz'
      PRINT *, 'Eigenvector '
      WRITE(LOUT,'(1X,A)') 'Eigenvector'
      PRINT '(1X,6E12.4)', (EVC(J, II),J=1,NQ)
      WRITE(LOUT,'(1X,6E12.4)') (EVC(J, II),J=1,NQ)
  320 CONTINUE
      GO TO 350
  330 PRINT *, 'NO CONVERGENCE'
      WRITE(LOUT,'()') 'No Convergence'
      GO TO 350
  340 PRINT *, 'SQUARE ROOT OF NEGATIVE TERM -- CHECK MATRICES'
      WRITE(LOUT,'()') 'Square Root of Negative Term -- Check Matrices'
  350 CLOSE(LOUT)
      END
