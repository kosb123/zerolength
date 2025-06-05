     '********    PROGRAM GAUSS   **********
     '*      GAUSS ELIMINATION METHOD      *
     '*          GENERAL MATRIX            *
     '* T.R.Chandrupatla and A.D.Belegundu *
     '**************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "GAUSS ELIMINATION METHOD to solve Ax=B";
     PRINT SPACE$(13); "(C)  Chandrupatla & Belegundu": COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Input Data File Name <DOS file name> ", FILE1$
     INPUT "Output Data File Name ", FILE2$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, N
     DIM A(N, N), B(N)
     LINE INPUT #1, D$
     '--- Matrix A() in  Ax = B ---
     FOR I = 1 TO N
        FOR J = 1 TO N
           INPUT #1, A(I, J)
        NEXT J
     NEXT I
     LINE INPUT #1, D$
'--- Right hand side B() in Ax = B ---
     FOR I = 1 TO N
        INPUT #1, B(I)
     NEXT I
     '----- Forward Elimination -----
     FOR K = 1 TO N - 1
        FOR I = K + 1 TO N
           C = A(I, K) / A(K, K)
           FOR J = K + 1 TO N
              A(I, J) = A(I, J) - C * A(K, J)
           NEXT J
           B(I) = B(I) - C * B(K)
        NEXT I
     NEXT K
     '----- Back-substitution -----
     B(N) = B(N) / A(N, N)
     FOR II = 1 TO N - 1
        I = N - II
        C = 1 / A(I, I): B(I) = C * B(I)
        FOR K = I + 1 TO N
           B(I) = B(I) - C * A(I, K) * B(K)
        NEXT K
     NEXT II
     OPEN FILE2$ FOR OUTPUT AS #2
     PRINT TITLE$
     PRINT #2, TITLE$
     PRINT "Solution Vector"
     PRINT #2, "Solution Vector"
     FOR I = 1 TO N
        PRINT B(I)
        PRINT #2, B(I)
     NEXT I
     CLOSE #2
     END

