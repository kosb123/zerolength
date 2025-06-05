C   ****  Sky Line Method  ****
C   **** Gauss Elimination ****
      DIMENSION ID(200), A(5000), B(200)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY,TITLE
      NMAX = 5000
      PRINT *,' ======================================='
      PRINT *,'            PROGRAM SKYLINE             '
      PRINT *,'   T.R.Chandrupatla and A.D.Belegundu   '
      PRINT *,' ======================================='
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
      NSUM = 0
C   *** Read Column Heights then convert, Diagonal Pointers ***
      READ(LINP,*) (ID(I),I=1,N)
      DO 100 I = 2, N
  100 ID(I) = ID(I) + ID(I - 1)
      READ(LINP,'(A)')DUMMY
      NSUM = ID(N)      
      NT = 0
      NI = 1
      DO 110 I = 1, N
      IF( I .EQ. 1 ) GO TO 105
      NT = ID(I - 1)
      NI = ID(I) - ID(I - 1)
  105 IF(I.EQ.1) NI = 1
  110 READ(LINP,*) (A(NT + J), J = 1, NI)
      READ(LINP,'(A)')DUMMY
C     ----- Right Hand Side b
      READ(LINP,*) (B(I),I=1,N)
C   --- Forward Elimination
      DO 180 J = 2, N
      NJ = ID(J) - ID(J - 1)
      IF(NJ .EQ. 1) GO TO 180
      K1 = 0
      NJ = J - NJ + 1
      DO 170 K = NJ, J - 1
      K1 = K1 + 1
      KJ = ID(J - 1) + K1
      KK = ID(K)
      C = A(KJ) / A(KK)
      DO 160 I = K + 1, J
      NI = ID(I) - ID(I - 1)
      IF((I - K + 1) .GT. NI ) GO TO 160
      IJ = ID(J) - J + I
      KI = ID(I) - I + K
      A(IJ) = A(IJ) - C * A(KI)
  160 CONTINUE
  170 CONTINUE
  180 CONTINUE
      N1 = N - 1
      DO 190 K = 1, N1
      KK = ID(K)
      K1 = K + 1
      C = B(K) / A(KK)
      DO 190 I = K1, N
      NI = ID(I) - ID(I - 1)
      IF (I - K + 1.GT. NI) GO TO 190
      KI = ID(I) - I + K
      B(I) = B(I) - C * A(KI)
  190 CONTINUE
C    **** Backsubstitution ****
      NS = ID(N)
      B(N) = B(N) / A(NS)
      DO 220 I1 = 1, N1
      I = N - I1
      II = ID(I)
      C = 1 / A(II)
      B(I) = C * B(I)
      II1=I + 1
      DO 210 J = II1, N
      J1 = J - I + 1
      NJ = ID(J) - ID(J - 1)
      IF( J1 .GT. NJ ) GO TO 210
      IJ = ID(J) - J + I
      B(I) = B(I) - C * A(IJ) * B(J)
  210 CONTINUE
  220 CONTINUE
      PRINT *, 'Solution Vector'
      WRITE (*,'(1X,I4,E12.4)')(I, B(I),I = 1, N)
      WRITE(LOUT,*) 'Solution Vector'
      DO 230 I = 1, N
  230 WRITE(LOUT,'(1X,E15.5)') B(I)
      END
