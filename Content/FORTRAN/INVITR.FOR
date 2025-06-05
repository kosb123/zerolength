
C     *********            INVITR             **********
      DIMENSION S(200,50),GM(200,50),EVC(200,50),EV1(200)
      DIMENSION EV2(200),EVT(200),EVS(200),ST(200),EVL(50)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY
      IMAX = 200
      PRINT *, '*****      PROGRAM INVITR         *****'
      PRINT *, '*      Inverse Iteration Method       *'
      PRINT *, '*  for eigenvalues and eigenvectors   *'
      PRINT *, '*        Searching in Subspace        *'
      PRINT *, '*         for Banded Matrices         *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'
      PRINT *, 'Give Name of Input File '
      READ '(A)', FILE1
      LINP = 10
      OPEN (UNIT = 10, FILE = FILE1, STATUS = 'OLD')      
      PRINT *, 'Give Name of Output File '
      READ '(A)', FILE2
      LOUT = 11
      OPEN (UNIT = 11, FILE = FILE2)
C  --- Read in Number of Equations and Bandwidth ---
      READ(LINP,'(A)') DUMMY
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) NQ, NBW
      TOL = .000001
      ITMAX = 50
      SH = 0
      NEV = 0
      PI = 3.14159
      PRINT *, '       TOLERANCE'
      PRINT *, ' 1. Default < 1E-6>'
      PRINT *, ' 2. Specify Tolerance'
      PRINT *, '  Your choice <1 or 2> '
      READ *,IFL1
      IF(IFL1.EQ.1) GO TO 100
      PRINT *, ' Tolerance desired '
      READ *, TOL
  100 PRINT '(A)',' Number of Eigenvalues Desired'
      READ *, NEV
C  --- Read in Stiffness Matrix ---
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) ((S(I, J),J=1,NBW),I=1,NQ)
C  === Read in Mass Matrix ---
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) ((GM(I, J),J=1,NBW),I=1,NQ)
      READ(LINP,'(A)') DUMMY
      READ(LINP,*) (ST(I),I=1,NQ)
      CLOSE(LINP)
      PRINT *, ' SHIFT VALUE FOR EIGENVALUES <Type 0 for zero shift>'
      READ *, SH
      If (SH .NE. 0) THEN
         DO 110 I = 1, NQ
           DO 110 J = 1, NBW
             S(I, J) = S(I, J) - SH * GM(I, J)
  110 CONTINUE
      END IF
C  --- Reduce Stiffness to Upper Triangle ---
      CALL REDLHS(S,NQ,NBW,IMAX)
      DO 250 NV = 1, NEV
C  --- Starting Value for Eigenvector ---
      DO 130 I = 1, NQ
  130 EV1(I) = ST(I)
      EL2 = 0.
      ITER = 0
  132 EL1 = EL2
      ITER = ITER + 1
      IF(ITER .GT. ITMAX )GO TO 260
      IF(NV .EQ. 1 )GO TO 162
C  ----  Starting Vector Orthogonal to  ----
C  ----       Evaluated Vectors         ----
      NV1 = NV - 1
      DO 160 I = 1, NV1
      CV = 0.
      DO 140 K = 1, NQ
      KA = K - NBW + 1
      KZ = K + NBW - 1
      IF(KA .LT. 1 )KA = 1
      IF(KZ .GT. NQ )KZ = NQ
      DO 140 L = KA, KZ
      IF(L .LT. K )GO TO 135
      K1 = K
      L1 = L - K + 1
      GO TO 140
  135 K1 = L
      L1 = K - L + 1
  140 CV = CV + EVS(K) * GM(K1, L1) * EVC(L, I)
      DO 150 K = 1, NQ
  150 EV1(K) = EV1(K) - CV * EVC(K, I)
  160 CONTINUE
  162 DO 180 I = 1, NQ
      IA = I - NBW + 1
      IZ = I + NBW - 1
      EVT(I) = 0.
      IF(IA .LT. 1 )IA = 1
      IF(IZ .GT. NQ )IZ = NQ
      DO 170 K = IA, IZ
      IF(K .LT. I )GO TO 165
      I1 = I
      K1 = K - I + 1
      GO TO 170
  165 I1 = K
      K1 = I - K + 1
  170 EVT(I) = EVT(I) + GM(I1, K1) * EV1(K)
      EV2(I) = EVT(I)
  180 CONTINUE
C --- Reduce Right Side and Solve ---
      CALL BACSUB(S,NQ,NBW,EV2,IMAX)
      C1 = 0.
      C2 = 0.
      DO 190 I = 1, NQ
  190 C1 = C1 + EV2(I) * EVT(I)
      DO 200 I = 1, NQ
      IA = I - NBW + 1
      IZ = I + NBW - 1
      EVT(I) = 0.
      IF(IA .LT. 1 )IA = 1
      IF(IZ .GT. NQ )IZ = NQ
      DO 200 K = IA, IZ
      IF(K .LT. I )GO TO 195
      I1 = I
      K1 = K - I + 1
      GO TO 200
  195 I1 = K
      K1 = I - K + 1
  200 EVT(I) = EVT(I) + GM(I1, K1) * EV2(K)
      DO 210 I = 1, NQ
  210 C2 = C2 + EV2(I) * EVT(I)
      EL2 = C1 / C2
      C2 = SQRT(C2)
      DO 220 I = 1, NQ
      EV1(I) = EV2(I) / C2
  220 EVS(I) = EV1(I)
      IF(ABS(EL2 - EL1) / ABS(EL2) .GT. TOL) GO TO 132
      DO 230 I = 1, NQ
  230 EVC(I, NV) = EV1(I)
      PRINT '(2(1X,A,I4))' , 'Eigenvalue Number', NV , '    Iteration Numb
     1er', ITER
      WRITE(LOUT,'(2(1X,A,I5))') 'Eigenvalue Number', NV, '    Iteration Nu
     1mber', ITER
      EL2 = EL2 + SH
      EVL(NV) = EL2
      OMEGA = SQRT(EL2)
      FREQ = .5 * OMEGA / PI
      PRINT '(3(1X,A,E12.5),A)', 'Eigenvalue = ', EL2, '    Omega = ', O
     1MEGA, '    Freq = ', FREQ, ' Hz'
      WRITE(LOUT,'(3(1X,A,E12.5),A)') 'Eigenvalue = ', EL2, '    Omega =
     1', OMEGA, '	  Freq = ', FREQ, ' Hz'
      PRINT *, 'Eigenvector '
      WRITE(LOUT,'(1X,A)') 'Eigenvector'
      PRINT '(1X,6E12.4)', (EVC(I, NV),I=1,NQ)
      WRITE(LOUT,'(1X,6E12.4)') (EVC(I, NV),I=1,NQ)
  250 CONTINUE
      GO TO 270
  260 PRINT *, 'No Convergence for ', ITER, ' Iterations'
  270 CLOSE(LOUT)
      END
      SUBROUTINE REDLHS(A,N,NBW,IMAX)
      DIMENSION A(IMAX,NBW)
C     ---- Gauss Elimination LDU Approach ----
C   --- REDUCTION TO UPPER TRIANGULAR MATRIX ---
      N1 = N - 1 
      DO 2000 K = 1, N1
      NK = N - K + 1
      IF(NK .GT. NBW )NK = NBW
      DO 2000 I = 2, NK
      C1 = A(K, I) / A(K, 1)
      I1 = K + I - 1
      DO 2000 J = I, NK
      J1 = J - I + 1
      A(I1, J1) = A(I1, J1) - C1 * A(K, J)
 2000 CONTINUE
      RETURN
      END
      SUBROUTINE BACSUB(A,N,NBW,B,IMAX)
      DIMENSION A(IMAX,NBW),B(N)
C   ---- Reduction of the right hand side ----
      N1 = N - 1
      DO 2010 K = 1, N1
      NK = N - K + 1
      IF(NK .GT. NBW )NK = NBW
      DO 2010 I = 2, NK
      I1 = K + I - 1
      C1 = 1 / A(K, 1)
      B(I1) = B(I1) - C1 * A(K, I) * B(K)
 2010 CONTINUE
C   ---- Back Substitution ----
      B(N) = B(N) / A(N, 1)
      DO 2020 II = 1, N1
      I = N - II
      C1 = 1 / A(I, 1)
      NI = N - I + 1
      IF(NI .GT. NBW )NI = NBW
      B(I) = C1 * B(I)
      DO 2020 K = 2, NI
      B(I) = B(I) - C1 * A(I, K) * B(I + K - 1)
 2020 CONTINUE
      RETURN
      END
