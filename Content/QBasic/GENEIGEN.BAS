DECLARE SUB EIGENVEC (A#(), B#(), N%)
DECLARE SUB TRIDIAG (A#(), D#(), B#(), N%)
DECLARE SUB EIGENTD (A#(), B#(), Q#(), N%, ITER%)
DECLARE SUB UPDTSTIFF (A#(), B#(), N%)
DECLARE SUB CHOLESKY (A#(), N%)
     '****** GENERALIZED EIGENVALUE PROBLEM ******
     '*              Kx = lambda Mx              *
     '*             T.R.Chandrupatla             *
     '********************************************
     DEFINT I-N
     DEFDBL A-H, O-Z
     CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "TRIDIAGONALIZATION & WILKINSON SHIFT"; SPACE$(15);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Name of Input File  ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     INPUT "Name of Output File  ", FILE2$
     OPEN FILE2$ FOR OUTPUT AS #2
     '--- Read in Number of Equations
     LINE INPUT #1, D$
     LINE INPUT #1, D$: INPUT #1, NQ, NBW
     DIM S(NQ, NQ), GM(NQ, NQ), D(NQ), B(NQ - 1), NORD(NQ)
     '----- Banded Stiffness Matrix into S(NQ,NQ) -----
     LINE INPUT #1, D$
     FOR I = 1 TO NQ: FOR JN = 1 TO NBW
        INPUT #1, STIFF: J = I + JN - 1
        IF J <= NQ THEN
        S(I, J) = STIFF: S(J, I) = STIFF
        END IF
     NEXT JN: NEXT I
     '----- Banded Mass Matrix into GM(NQ,NQ) -----
     LINE INPUT #1, D$
     FOR I = 1 TO NQ: FOR JN = 1 TO NBW
        INPUT #1, EMASS: J = I + JN - 1
        IF J <= NQ THEN
           GM(I, J) = EMASS: GM(J, I) = EMASS
        END IF
     NEXT JN: NEXT I
     CLOSE #1
     '----- Cholesky Factorization of Mass Matrix
     CALL CHOLESKY(GM(), NQ)
     '----- Update of Stiffness Matrix - Standard Form Ax=(lambda)x
     CALL UPDTSTIFF(S(), GM(), NQ)
     '----- Tri-diagonalize D() Diagonal, B() Sub-diagonal
     '----- S() has the rotation matrix
     CALL TRIDIAG(S(), D(), B(), NQ)
     '----- Find Eigenvalues and Eigenvectors
     CALL EIGENTD(D(), B(), S(), NQ, ITER)
     '----- Determine Eigenvectors
     CALL EIGENVEC(S(), GM(), NQ)
     '----- Print Eigenvalues and Eigenvectors
     FOR I = 1 TO NQ
        PRINT "Eigenvalue = "; D(I)
        PRINT "Eigenvector"
        FOR J = 1 TO NQ: PRINT S(J, I); : NEXT J: PRINT
     NEXT I
     '-----   RESULTS   -----
     '--- Ascending Order of Eigenvalues
     FOR I = 1 TO NQ: NORD(I) = I: NEXT I
     FOR I = 1 TO NQ - 1
        II = NORD(I): I1 = II
        C1 = D(II): J1 = I
        FOR J = I + 1 TO NQ
           IJ = NORD(J)
           IF C1 > D(IJ) THEN
              C1 = D(IJ): I1 = IJ: J1 = J
           END IF
        NEXT J
        NORD(I) = I1: NORD(J1) = II
     NEXT I
     PRINT "Eigenvalues and Eigenvectors for Data in File "; FILE1$
     PRINT #2, "Eigenvalues and Eigenvectors for Data in File "; FILE1$
     PI = 3.14159
     FOR I = 1 TO NQ
        II = NORD(I)
        PRINT "Eigenvalue Number  "; I
        PRINT #2, "Eigenvalue Number  "; I
        C = D(II): OMEGA = SQR(C): FREQ = .5 * OMEGA / PI
        PRINT USING "Eigenvalue = ##.###^^^^"; C;
        PRINT #2, USING "Eigenvalue = ##.####^^^^"; C;
        PRINT USING "   Omega = ##.###^^^^"; OMEGA;
        PRINT USING "   Freq = ##.###^^^^ Hz"; FREQ
        PRINT #2, USING "   Omega = ##.###^^^^"; OMEGA;
        PRINT #2, USING "   Freq = ##.###^^^^ Hz"; FREQ
        PRINT "Eigenvector ": PRINT #2, "Eigenvector "
        FOR J = 1 TO NQ
           PRINT USING "##.###^^^^ "; S(J, II);
           PRINT #2, USING "##.###^^^^ "; S(J, II);
           IFL = J - 7 * INT(J / 7)
           IF IFL = 0 THEN PRINT : PRINT #2,
        NEXT J
        PRINT : PRINT #2,
     NEXT I
     CLOSE #2
     END

     SUB CHOLESKY (A(), N)
     '----- L into lower left triangle of A
     FOR K = 1 TO N
        A(K, K) = SQR(A(K, K))
        FOR I = K + 1 TO N
           A(I, K) = A(I, K) / A(K, K)
        NEXT I
        FOR J = K + 1 TO N
           FOR I = J TO N
              A(I, J) = A(I, J) - A(I, K) * A(J, K)
           NEXT I
        NEXT J
     NEXT K
     '----- Inverse placed in upper triangle
     '---   Use 1/diag for inverse diagonal
     FOR K = 1 TO N - 1
        FOR I = 1 TO N - K
           IK = I + K
           A(I, IK) = 0
           FOR J = I TO IK - 1
              C = A(I, J): IF I = J THEN C = 1 / C
              A(I, IK) = A(I, IK) - A(IK, J) * C
           NEXT J
           A(I, IK) = A(I, IK) / A(IK, IK)
        NEXT I
     NEXT K
     END SUB

SUB EIGENTD (A(), B(), Q(), N, ITER)
     ITER = 0
     M = N
  DO
     ITER = ITER + 1
     D = .5 * (A(M - 1) - A(M))
     BB = B(M - 1) * B(M - 1)
     BOT = D + SGN(D) * SQR(D * D + BB)
     P = A(1) - A(M) + BB / BOT
     X = B(1)
     FOR I = 1 TO M - 1
        PP = SQR(P * P + X * X)
        C = -P / PP
        S = X / PP
        IF I > 1 THEN B(I - 1) = C * B(I - 1) - S * X
        A1 = A(I): A2 = A(I + 1): B1 = B(I)
        A(I) = A1 * C * C - 2 * B1 * C * S + A2 * S * S
        B(I) = (A1 - A2) * C * S + B1 * (C * C - S * S)
        A(I + 1) = A1 * S * S + 2 * B1 * C * S + A2 * C * C
        '----- Update Q()
        FOR K = 1 TO N
           A1 = Q(K, I): A2 = Q(K, I + 1)
           Q(K, I) = C * A1 - S * A2
           Q(K, I + 1) = S * A1 + C * A2
        NEXT K
        IF I = M - 1 THEN EXIT FOR
        X = -B(I + 1) * S
        B(I + 1) = B(I + 1) * C
        P = B(I)
     NEXT I
     DO WHILE M > 1 AND ABS(B(M - 1)) < .000001
           M = M - 1
     LOOP
  LOOP WHILE M > 1
  M = M - 1
END SUB

SUB EIGENVEC (A(), B(), N)
     FOR I = 1 TO N
        FOR J = 1 TO N
           T = 0
           FOR K = I TO N
              C = B(I, K)
              IF I = K THEN C = 1 / C
              T = T + C * A(K, J)
           NEXT K
           A(I, J) = T
        NEXT J
     NEXT I
END SUB

SUB TRIDIAG (A(), D(), B(), N)
    FOR I = 1 TO N - 2
       AA = 0
       FOR J = I + 1 TO N
          AA = AA + A(J, I) * A(J, I)
       NEXT J
       AA = SQR(AA)
       WW = 2 * AA * (AA + ABS(A(I + 1, I)))
       WW = SQR(WW)
       IA = SGN(A(I + 1, I))
       '----- Diagonal and Next to Diagonal Term
       D(I) = A(I, I)
       B(I) = -IA * AA
       '----- Unit Vector W() in Column I from Row I+1 to N
       FOR J = I + 1 TO N
          A(J, I) = A(J, I) / WW
       NEXT J
       A(I + 1, I) = A(I + 1, I) + IA * AA / WW
       '----- W'A in Row I from Col I+1 to N
       BET = 0
       FOR J = I + 1 TO N
          A(I, J) = 0
          FOR K = I + 1 TO N
              A(I, J) = A(I, J) + A(K, I) * A(K, J)
          NEXT K
          BET = BET + A(I, J) * A(J, I)
       NEXT J
       '----- Modified A()                                              
       FOR J = I + 1 TO N
          FOR K = I + 1 TO N
             A(J, K) = A(J, K) - 2 * A(J, I) * A(I, K) - 2 * A(I, J) * A(K, I)
             A(J, K) = A(J, K) + 4 * BET * A(K, I) * A(J, I)
          NEXT K
       NEXT J
    NEXT I
    D(N - 1) = A(N - 1, N - 1)
    B(N - 1) = A(N - 1, N)
    D(N) = A(N, N)
    A(N - 1, N - 1) = 1
    A(N, N) = 1
    A(N - 1, N) = 0
    A(N, N - 1) = 0
    '----- Now Create the Q matrix in A()
    FOR I = 1 TO N - 2
       II = N - I - 1
       A(II, II) = 1
       FOR J = 1 TO I + 1
          IJ = II + J
          A(II, IJ) = 0
          FOR K = 1 TO I + 1
             IK = II + K
             A(II, IJ) = A(II, IJ) + A(IK, IJ) * A(IK, II)
          NEXT K
          FOR K = 1 TO I + 1
             IK = II + K
             A(IK, IJ) = A(IK, IJ) - 2 * A(IK, II) * A(II, IJ)
          NEXT K
       NEXT J
       FOR J = 1 TO I + 1
          IJ = II + J
          A(II, IJ) = 0
          A(IJ, II) = 0
       NEXT J
    NEXT I
END SUB

SUB UPDTSTIFF (A(), B(), N)
     FOR J = N TO 1 STEP -1
        FOR I = 1 TO N
           T = 0
           FOR K = 1 TO J
              C = B(K, J)
              IF K = J THEN C = 1 / C
              T = T + C * A(I, K)
           NEXT K
           A(I, J) = T
        NEXT I
     NEXT J
     FOR I = N TO 1 STEP -1
        FOR J = 1 TO N
           T = 0
           FOR K = 1 TO I
              C = B(K, I)
              IF K = I THEN C = 1 / C
              T = T + C * A(K, J)
           NEXT K
           A(I, J) = T
        NEXT J
     NEXT I
END SUB

