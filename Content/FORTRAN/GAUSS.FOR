C   *****  Gauss Elimination Method (General Matrix) *****
C   *****           Row-wise Elimination             *****
      DIMENSION A(100, 100), B(100)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '======================================='
      PRINT *, '           PROGRAM GAUSS               '
      PRINT *, '  T.R.Chandrupatla and A.D.Belegundu   '
      PRINT *, '======================================='
      PRINT *,'Input Data File Name <DOS file name>'
      READ '(A)',FILE1
      LINP=10
      OPEN(UNIT=10,FILE=FILE1,STATUS='UNKNOWN')
      PRINT *,'Output Data File Name <DOS file name>'
      READ '(A)',FILE2
      LOUT=11
      OPEN(UNIT=11,FILE=FILE2,STATUS='UNKNOWN')
      READ(LINP,'(A)')TITLE
      READ(LINP,'(A)')DUMMY
C     ----- Number of Equations N
      READ(LINP,*) N
      READ(LINP,'(A)')DUMMY
C      ----- A in Ax = B -----
      READ(LINP,*) ((A(I,J),J=1,N),I=1,N)
      READ(LINP,'(A)')DUMMY
C     ----- Right Hand Side b
      READ(LINP,*) (B(I),I=1,N)
      N1 = N - 1
      DO 130 K = 1, N1
      K1 = K + 1
      DO 130 I = K1, N
      C = A(I, K) / A(K, K)
      DO 120 J = K1, N
  120 A(I, J) = A(I, J) - C * A(K, J)
  130 B(I) = B(I) - C * B(K)
C   *** Back Substitution ***
      B(N) = B(N) / A(N, N)
      DO 140 II = 1, N1
      I = N - II
      I1 = I + 1
      C = 1 / A(I, I)
      B(I) = C * B(I)
      DO 140 K = I1, N
  140 B(I) = B(I) - C * A(I, K) * B(K)
      PRINT *, 'Solution Vector'
      WRITE (*,'(1X,I4,E12.4)')(I, B(I),I = 1, N)
      WRITE(LOUT,*) 'Solution Vector'
      DO 150 I = 1, N
  150 WRITE(LOUT,'(1X,E15.5)') B(I)
      END
  
