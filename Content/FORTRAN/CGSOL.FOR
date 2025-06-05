      DIMENSION A(100, 100), B(100), X(100)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '*********   PROGRAM CGSOL   **********'
      PRINT *, '*     CONJUGATE GRADIENT METHOD      *'
      PRINT *, '*   FOR SOLVING AX=B, A Symmetric    *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu *'
      PRINT *, '**************************************'
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
      ITER = 0
      CALL CG(A, B, X, N, ITER)
      PRINT *, 'Iterations = ', ITER
      PRINT *, 'Solution Vector'
      WRITE (*,'(I4,E12.4)')(I, X(I),I = 1, N)
      END

      SUBROUTINE CG (A, B, X, N, ITER)
      DIMENSION A(100,100), B(100), X(100), G(100), D(100), AD(100)
      DO 100 I = 1, N
	X(I) = 0.
        G(I) = -B(I)
        D(I) = B(I)
  100 CONTINUE
      GG1 = 0.
      DO 110 I = 1, N
         GG1 = GG1 + G(I) * G(I)
  110 CONTINUE
  120 IF (GG1.LE.1E-6) GO TO 170
	ITER = ITER + 1
	DAD = 0.
        DO 140 I = 1, N
	   C = 0.
           DO 130 J = 1, N
              C = C + A(I, J) * D(J)
  130      CONTINUE
           AD(I) = C
           DAD = DAD + C * D(I)
  140   CONTINUE
        AL = GG1 / DAD
	GG2 = 0.
        DO 150 I = 1, N
           X(I) = X(I) + AL * D(I)
           G(I) = G(I) + AL * AD(I)
           GG2 = GG2 + G(I) * G(I)
  150   CONTINUE
        BT = GG2 / GG1
        DO 160 I = 1, N
           D(I) = -G(I) + BT * D(I)
  160	CONTINUE
        GG1 = GG2
      GO TO 120
  170 CONTINUE
      END
