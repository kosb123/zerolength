C   **********	 PROGRAM GENEIGEN   ********
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION S(100,100),GM(100,100),NORD(100),D(100),B(99)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY
      PRINT *, '********************************************'
      PRINT *, '*      GENERALIZED EIGENVALUE PROBLEM      *'
      PRINT *, '*              Kx = lambda Mx              *'
      PRINT *, '*    IMPLICIT SHIFT METHOD OF WILKINSON    *'
      PRINT *, '*    T.R.Chandrupatla and A.D.Belegundu    *'
      PRINT *, '********************************************'
      IMAX = 100
      PI = 3.14159
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
C   --- Banded Stiffness Matrix into S(NQ,NQ) ---
      READ(LINP,'(A)') DUMMY
C   --- Array D used for reading
      DO 110 I = 1, NQ
      READ(LINP,*) (D(JN),JN=1,NBW)
      DO 110 JN = 1, NBW
      J = I + JN - 1
      IF(J .GT. NQ) GO TO 110
      S(I, J) = D(JN)
      S(J, I) = D(JN)
  110 CONTINUE
C   --- Read in Mass Matrix into GM(NQ,NQ) ---
      READ(LINP,'(A)') DUMMY
      DO 120 I = 1, NQ
      READ(LINP,*) (D(JN),JN=1,NBW)
      DO 120 JN = 1, NBW
      J = I + JN - 1
      IF(J .GT. NQ) GO TO 120
      GM(I, J) = D(JN)
      GM(J, I) = D(JN)
  120 CONTINUE
      CLOSE(LINP)
C ----- Cholesky Factorization of Mass Matrix
      CALL CHLSKY(GM, NQ)
C ----- Update of Stiffness Matrix - Standard Form Ax=(lambda)x
      CALL UPDSTF(S, GM, NQ)
C ----- Tri-diagonalize D() Diagonal, B() Sub-diagonal
C ----- S() has the rotation matrix
      CALL TRDIAG(S, D, B, NQ)
C ----- Find Eigenvalues and Eigenvectors
      CALL EIGNTD(D, B, S, NQ, ITER)
C ----- Determine Eigenvectors
      CALL EIGNVC(S, GM, NQ)
C ----- Ascending Order of Eigenvalues
      CALL ASCEND(D, NORD, NQ)
C ----- Eigenvalues and Eigenvectors -----
      DO 320 I = 1, NQ
      II = NORD(I)
      C = D(II)
      OMEGA = SQRT(C)
      FREQ = .5 * OMEGA / PI
      PRINT '(1X,A,I4)' , ' Eigenvalue Number', I
      WRITE(LOUT,'(1X,A,I4)') 'Eigenvalue Number', I
      PRINT '(3(1X,A,E12.5),A)', ' Eigenvalue = ', C, '    Omega = ',
     1 OMEGA, '    Freq = ', FREQ, ' Hz'
      WRITE(LOUT,'(3(1X,A,E12.5),A)') 'Eigenvalue = ', C, '    Omega
     1= ', OMEGA, '    Freq = ', FREQ, ' Hz'
      PRINT *, 'Eigenvector '
      WRITE(LOUT,'(1X,A)') 'Eigenvector'
      PRINT '(1X,6E12.4)', (S(J, II),J=1,NQ)
      WRITE(LOUT,'(1X,6E12.4)') (S(J, II),J=1,NQ)
  320 CONTINUE
      CLOSE(LOUT)
      PRINT *, 'Results are in File ', FILE2
      END

      SUBROUTINE CHLSKY(A,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(100,100)
C ----- L into lower left triangle of A
      DO 1020 K = 1, N
      A(K, K) = SQRT(A(K, K))
      DO 1010 I = K + 1, N
 1010 A(I, K) = A(I, K) / A(K, K)
      DO 1020 J = K + 1, N
      DO 1020 I = J, N
 1020 A(I, J) = A(I, J) - A(I, K) * A(J, K)
      RETURN
      END

      SUBROUTINE UPDSTF(A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(100,100),B(100,100)
C --- Forward Substitution I  (invL)*A
      DO 2020 J = 1, N
        A(1, J) = A(1, J) / B(1, 1)
        DO 2020 I = 2, N
        DO 2010 K = 1, I - 1
        A(I, J) = A(I, J) - B(I, K) * A(K, J)
2010    CONTINUE 
      A(I, J) = A(I, J) / B(I, I)
2020  CONTINUE
C --- Forward Substitution II   (invL)*A*(invL)'
      DO 2040 J = 1, N
        A(J, 1) = A(J, 1) / B(1, 1)
      DO 2040 I = 2, N
          DO 2030 K = 1, I - 1
              A(J, I) = A(J, I) - B(I, K) * A(J, K)
 2030 CONTINUE
         A(J, I) = A(J, I) / B(I, I)
 2040 CONTINUE
      RETURN
      END

      SUBROUTINE TRDIAG(A, D, B, N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(100,100),D(N),B(N-1)
      DO 3005 I= 1, N - 1
      B(I) = 0.
 3005 D(I) = 0.
      D(N) = 0.
      DO 3050 I = 1, N - 2
       AA = 0.
       DO 3010 J = I + 1, N
          AA = AA + A(J, I) * A(J, I)
 3010  CONTINUE
       AA = SQRT(AA)
       WW = 2 * AA * (AA + ABS(A(I + 1, I)))
       WW = SQRT(WW)
C ----- Diagonal and Next to Diagonal Term
       D(I) = A(I, I)
       IA = 1
       IF (A(I+1, I).LT.0) IA = -1
       B(I) = -IA * AA
C ----- Unit Vector W() in Column I from Row I+1 to N
       DO 3020 J = I + 1, N
          A(J, I) = A(J, I) / WW
 3020  CONTINUE
       A(I + 1, I) = A(I + 1, I) +IA * AA / WW
C ----- W'A in Row I from Col I+1 to N
       BET = 0.
       DO 3040 J = I + 1, N
	  A(I, J) = 0.
          DO 3030 K = I + 1, N
              A(I, J) = A(I, J) + A(K, I) * A(K, J)
 3030     CONTINUE
          BET = BET + A(I, J) * A(J, I)
 3040  CONTINUE
C ----- Modified A()
       DO 3050 J = I + 1, N
          DO 3050 K = I + 1, N
          A(J,K) = A(J,K) - 2*A(J,I)*A(I,K) - 2*A(I,J)*A(K, I)
     1             + 4*BET*A(K, I)*A(J,I)
 3050  CONTINUE
      D(N - 1) = A(N - 1, N - 1)
      B(N - 1) = A(N - 1, N)
      D(N) = A(N, N)
      A(N - 1, N - 1) = 1.
      A(N, N) = 1
      A(N - 1, N) = 0.
      A(N, N - 1) = 0.
C ----- Now Create the Q matrix in A()
      DO 3080 I = 1, N - 2
       II = N - I - 1
       A(II, II) = 1
       DO 3070 J = 1, I + 1
          IJ = II + J
          A(II, IJ) = 0
          DO 3060 K = 1, I + 1
             IK = II + K
             A(II, IJ) = A(II, IJ) + A(IK, IJ) * A(IK, II)
 3060     CONTINUE
          DO 3070 K = 1, I + 1
             IK = II + K
             A(IK, IJ) = A(IK, IJ) - 2 * A(IK, II) * A(II, IJ)
 3070     CONTINUE
       DO 3080 J = 1, I + 1
          IJ = II + J
          A(II, IJ) = 0
          A(IJ, II) = 0
 3080  CONTINUE
      RETURN
      END


      SUBROUTINE EIGNTD(A, B, Q, N, ITER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(100),B(99),Q(100,100)
      ITER = 0
      M = N
 4010 ITER = ITER + 1
      DD = 0.5 * (A(M - 1) - A(M))
      BB = B(M - 1) * B(M - 1)
      BOT = DD + SIGN(SQRT(DD * DD + BB),DD)
      P = A(1) - A(M) + BB / BOT
      X = B(1)
      DO 4030 I = 1, M - 1
        PP = SQRT(P * P + X * X)
        CS = -P / PP
        SN = X / PP
        IF (I .GT. 1) B(I - 1) = CS * B(I - 1) - SN * X
        A1 = A(I)
        A2 = A(I + 1)
        B1 = B(I)
        A(I) = A1 * CS * CS - 2 * B1 * CS * SN + A2 * SN * SN
        B(I) = (A1 - A2) * CS * SN + B1 * (CS * CS - SN * SN)
        A(I + 1) = A1 * SN * SN + 2 * B1 * CS * SN + A2 * CS * CS
C ----- Update Q()
        DO 4020 K = 1, N
	   A1 = Q(K, I)
	   A2 = Q(K, I + 1)
           Q(K, I) = CS * A1 - SN * A2
           Q(K, I + 1) = SN * A1 + CS * A2
 4020  CONTINUE
        IF (I .EQ. M - 1) GO TO 4040
        X = -B(I + 1) * SN
        B(I + 1) = B(I + 1) * CS
        P = B(I)
 4030 CONTINUE
 4040 CONTINUE
      IF (M .LE. 1 .OR. ABS(B(M - 1)) .GT. 0.000001) GO TO 4050
	 M = M - 1
      GO TO 4040
 4050 IF (M .GT. 1) GO TO 4010
      RETURN
      END

      SUBROUTINE EIGNVC(A, B, N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(100,100),B(100,100)
C --- Backsubstitution ---
      DO 5020 J = 1, N
	A(N, J) = A(N, J) / B(N, N)
	N1 = N - 1
	DO 5020 I = N1, 1, -1
	   I1 = I + 1
	   DO 5010 K = N, I1, -1
              A(I, J) = A(I, J) - B(K, I) * A(K, J)
 5010   CONTINUE
        A(I, J) = A(I, J) / B(I, I)
 5020 CONTINUE
      RETURN
      END

      SUBROUTINE ASCEND(D, NORD, NQ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION D(100), NORD(100)
      DO 6010 I = 1, NQ
 6010 NORD(I) = I
C ----- Ascending Order of Eigenvalues
      N1 = NQ - 1
      DO 6030 I = 1, N1
        II = NORD(I)
        I1 = II
        C1 = D(II)
	J1 = I
	II1 = I + 1
	DO 6020 J = II1, NQ
           IJ = NORD(J)
	   IF (C1 .GT. D(IJ)) THEN
              C1 = D(IJ)
              I1 = IJ
              J1 = J
           END IF
 6020  CONTINUE
       NORD(I) = I1
       NORD(J1) = II
 6030 CONTINUE
      RETURN
      END
