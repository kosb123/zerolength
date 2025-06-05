C     ----------------------  HEAT1D  -----------------------
      DIMENSION X(50),S(50,20),TC(40),F(50),V(15),H(15),NB(15)
      CHARACTER*5, BC(15)
      CHARACTER*16 FILE1,FILE2
      CHARACTER*81 DUMMY,TITLE
      PRINT *, '***************************************'
      PRINT *, '*          PROGRAM  HEAT1D           *'
      PRINT *, '* T.R.Chandrupatla and A.D.Belegundu  *'
      PRINT *, '***************************************'      
C     IMAX = FIRST DIMENSION OF THE S-MATRIX
      IMAX = 50
      PRINT *,'Input Data File Name <DOS file name>'
      READ '(A)',FILE1
      LINP=10
      OPEN(UNIT=10,FILE = FILE1,STATUS='UNKNOWN')
      PRINT *,'Output Data File Name <DOS file name>'
      READ '(A)',FILE2
      LOUT=11
      OPEN(UNIT=11,FILE=FILE2,STATUS='UNKNOWN')
      READ(LINP,'(A)')DUMMY
      READ(LINP,'(A)')TITLE
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) NE, NBC, NQ
      NN = NE + 1
      NBW = 2    
C     NBW IS THE HALF-BAND-WIDTH
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) (N,TC(N),I=1,NE)
      READ(LINP,'(A)')DUMMY
      READ(LINP,*) (N,X(N),I=1,NN)
      READ(LINP,'(A)')DUMMY
      DO 9 I = 1, NBC
	READ(LINP,'(I5,A)') NB(I), BC(I)
	IF (BC(I) .EQ.'TEMP' .OR. BC(I).EQ.'temp') THEN
          READ(LINP,*) V(I)
	ELSEIF (BC(I) .EQ.'FLUX' .OR. BC(I) .EQ.'flux') THEN
          READ(LINP,*) V(I)
	ELSEIF (BC(I).EQ. 'CONV' .OR. BC(I) .EQ.'conv') THEN
          READ(LINP,*) H(I), V(I)
        ENDIF
9     CONTINUE
C     --- CALCULATE AND INPUT NODAL HEAT SOURCE VECTOR ---
      DO 11 I=1,NN
      F(I)=0.
      DO 11 J=1,NBW
      S(I,J)=0.
11    CONTINUE
      READ(LINP,'(A)')DUMMY
      IF (NQ .GT. 0) THEN
        DO 13 I = 1,NQ
          READ(LINP,*) N, F(N)
13      CONTINUE
      END IF
      CLOSE (LINP)
C     --- STIFFNESS MATRIX ---
      DO 20 I = 1,NE
      I1 = I
      I2 = I + 1
      ELL = ABS(X(I2) - X(I1))
      EKL = TC(I) / ELL
      S(I1, 1) = S(I1, 1) + EKL
      S(I2, 1) = S(I2, 1) + EKL
      S(I1, 2) = S(I1, 2) - EKL
20    CONTINUE
C     --- ACCOUNT FOR B.C.'S ---
      AMAX = 0
      DO 30 I = 1, NN
        IF (S(I, 1) .GT. AMAX) AMAX = S(I, 1)
30    CONTINUE
      CNST = AMAX * 10000.
      DO 40 I = 1,NBC
      N = NB(I)
      IF (BC(I) .EQ. 'CONV' .OR. BC(I) .EQ. 'conv') THEN
        S(N, 1) = S(N, 1) + H(I)
        F(N) = F(N) + H(I) * V(I)
      ELSEIF (BC(I) .EQ. 'FLUX' .OR. BC(I) .EQ. 'flux') THEN
        F(N) = F(N) - V(I)
      ELSE
        S(N, 1) = S(N, 1) + CNST
        F(N) = F(N) + CNST * V(I)
      END IF
40    CONTINUE
      CALL BANSOL(NN,NBW,IMAX,S,F)
C     --- F CONTAINS THE SOLUTION. 'RHS' IS OVER-WRITTEN
      WRITE(LOUT,'(1X,A)') TITLE
      WRITE(6,'(1X,A)') TITLE
      WRITE(6,*)'        NODE#  TEMPERATURE'
      WRITE(LOUT,*)'        NODE#  TEMPERATURE'
      DO 45 I = 1,NN
        PRINT *, I, F(I)
        WRITE(LOUT,*)I,F(I)
45    CONTINUE
      CLOSE (LOUT)
      WRITE(6,*)'Output is in file  ',FILE2
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
