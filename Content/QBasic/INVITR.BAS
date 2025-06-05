     '*****       PROGRAM INVITR        *****
     '*      Inverse Iteration Method       *
     '*  for eigenvalues and eigenvectors   *
     '*        Searching in Subspace        *
     '*         for Banded Matrices         *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "INVERSE ITERATION METHOD"; SPACE$(27);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Name of Input File  ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     INPUT "Name of Output File  ", FILE2$
     OPEN FILE2$ FOR OUTPUT AS #2
     '--- Read in Number of Equations
     LINE INPUT #1, D$
     LINE INPUT #1, D$: INPUT #1, NQ, NBW
     DIM S(NQ, NBW), GM(NQ, NBW), EV1(NQ), EV2(NQ)
     DIM EVT(NQ), EVS(NQ), ST(NQ)
     TOL = .000001: ITMAX = 50: SH = 0: NEV = 0: PI = 3.14159
     INPUT "Tolerance  < default = 1E-6 > ", A$
     IF VAL(A$) > 0 THEN TOL = VAL(A$)
     INPUT "Number of Eigenvalues Desired ", NEV
     DIM EVC(NQ, NEV), EVL(NEV)
     '----- Read in Stiffness Matrix
     LINE INPUT #1, D$
     FOR I = 1 TO NQ: FOR J = 1 TO NBW
        INPUT #1, S(I, J): NEXT J: NEXT I
     '----- Read in Mass Matrix
     LINE INPUT #1, D$
     FOR I = 1 TO NQ: FOR J = 1 TO NBW
        INPUT #1, GM(I, J): NEXT J: NEXT I
     '----- Starting Vector for Inverse Iteration
     LINE INPUT #1, D$
     For I = 1 To NQ: INPUT #1, ST(I): Next I
     CLOSE #1
     INPUT "Shift Value for Eigenvalues  < default = 0 > ", A$
     IF VAL(A$) > 0 THEN SH = VAL(A$)
     If SH <> 0 THEN
         FOR I = 1 TO NQ: FOR J = 1 TO NBW
             S(I, J) = S(I, J) - SH * GM(I, J)
         NEXT J: NEXT I
     END IF
     PRINT "Eigenvalues and Eigenvectors for Data in File "; FILE1$
     PRINT #2, "Eigenvalues and Eigenvectors for Data in File "; FILE1$
     GOSUB BANSOL2        '<----Stiffness to Upper Triangle
     FOR NV = 1 TO NEV
        '--- Starting Value for Eigenvector
        FOR I = 1 TO NQ
           EV1(I) = ST(I)
        NEXT I
        EL2 = 0: ITER = 0
        DO
           EL1 = EL2
           ITER = ITER + 1
           IF ITER > ITMAX THEN
              PRINT "No Convergence for "; ITER; " Iterations"
              END
           END IF
           IF NV > 1 THEN
              '----  Starting Vector Orthogonal to
              '----       Evaluated Vectors
              FOR I = 1 TO NV - 1
                 CV = 0
                 FOR K = 1 TO NQ
                    KA = K - NBW + 1: KZ = K + NBW - 1
                    IF KA < 1 THEN KA = 1
                    IF KZ > NQ THEN KZ = NQ
                    FOR L = KA TO KZ
                       IF L < K THEN
                          K1 = L: L1 = K - L + 1
                       ELSE
                          K1 = K: L1 = L - K + 1
                       END IF
                       CV = CV + EVS(K) * GM(K1, L1) * EVC(L, I)
                    NEXT L
                 NEXT K
                 FOR K = 1 TO NQ
                    EV1(K) = EV1(K) - CV * EVC(K, I)
                 NEXT K
              NEXT I
           END IF
           FOR I = 1 TO NQ
              IA = I - NBW + 1: IZ = I + NBW - 1: EVT(I) = 0
              IF IA < 1 THEN IA = 1
              IF IZ > NQ THEN IZ = NQ
              FOR K = IA TO IZ
                 IF K < I THEN
                    I1 = K: K1 = I - K + 1
                 ELSE
                    I1 = I: K1 = K - I + 1
                 END IF
                 EVT(I) = EVT(I) + GM(I1, K1) * EV1(K)
              NEXT K
              EV2(I) = EVT(I)
           NEXT I
           GOSUB RHSOLVE       '<--- Reduce Right Side and Solve
           C1 = 0: C2 = 0
           FOR I = 1 TO NQ
              C1 = C1 + EV2(I) * EVT(I)
           NEXT I
           FOR I = 1 TO NQ
              IA = I - NBW + 1: IZ = I + NBW - 1: EVT(I) = 0
              IF IA < 1 THEN IA = 1
              IF IZ > NQ THEN IZ = NQ
              FOR K = IA TO IZ
                 IF K < I THEN
                    I1 = K: K1 = I - K + 1
                 ELSE
                    I1 = I: K1 = K - I + 1
                 END IF
                 EVT(I) = EVT(I) + GM(I1, K1) * EV2(K)
              NEXT K
           NEXT I
           FOR I = 1 TO NQ
              C2 = C2 + EV2(I) * EVT(I)
           NEXT I
           EL2 = C1 / C2
           C2 = SQR(C2)
           FOR I = 1 TO NQ
              EV1(I) = EV2(I) / C2
              EVS(I) = EV1(I)
           NEXT I
        LOOP WHILE ABS(EL2 - EL1) / ABS(EL2) > TOL
        FOR I = 1 TO NQ
           EVC(I, NV) = EV1(I)
        NEXT I
        PRINT "Eigenvalue Number "; NV;
        PRINT #2, "Eigenvalue Number "; NV;
        PRINT "     Iteration Number "; ITER
        PRINT #2, "     Iteration Number "; ITER
        EL2 = EL2 + SH: EVL(NV) = EL2
        PRINT USING "Eigenvalue = ##.###^^^^"; EL2;
        PRINT #2, USING "Eigenvalue = ##.####^^^^"; EL2;
        OMEGA = SQR(EL2): FREQ = .5 * OMEGA / PI
        PRINT USING "   Omega = ##.###^^^^"; OMEGA;
        PRINT USING "   Freq = ##.###^^^^ Hz"; FREQ
        PRINT #2, USING "   Omega = ##.###^^^^"; OMEGA;
        PRINT #2, USING "   Freq = ##.###^^^^ Hz"; FREQ
        PRINT "Eigenvector ": PRINT #2, "Eigenvector "
        FOR I = 1 TO NQ
           PRINT USING "##.###^^^^ "; EVC(I, NV);
           PRINT #2, USING "##.###^^^^ "; EVC(I, NV);
           IFL = I - 7 * INT(I / 7)
           IF IFL = 0 THEN PRINT : PRINT #2,
        NEXT I
        PRINT : PRINT #2,
     NEXT NV
     CLOSE #2
     END
BANSOL2:
'-----  Gauss Elimination LDU Approach (for Symmetric Banded Matrices)
     '----- Reduction to Upper Triangular Form
     FOR K = 1 TO NQ - 1
        NK = NQ - K + 1
        IF NK > NBW THEN NK = NBW
        FOR I = 2 TO NK
           C1 = S(K, I) / S(K, 1)
           I1 = K + I - 1
           FOR J = I TO NK
              J1 = J - I + 1
              S(I1, J1) = S(I1, J1) - C1 * S(K, J)
           NEXT J
        NEXT I
     NEXT K
     RETURN
RHSOLVE:
     '----- Reduction of the right hand side
     FOR K = 1 TO NQ - 1
        NK = NQ - K + 1
        IF NK > NBW THEN NK = NBW
        FOR I = 2 TO NK: I1 = K + I - 1
           C1 = 1 / S(K, 1)
           EV2(I1) = EV2(I1) - C1 * S(K, I) * EV2(K)
        NEXT I
     NEXT K
     '----- Back Substitution
     EV2(NQ) = EV2(NQ) / S(NQ, 1)
     FOR II = 1 TO NQ - 1
        I = NQ - II: C1 = 1 / S(I, 1)
        NI = NQ - I + 1
        IF NI > NBW THEN NI = NBW
        EV2(I) = C1 * EV2(I)
        FOR K = 2 TO NI
           EV2(I) = EV2(I) - C1 * S(I, K) * EV2(I + K - 1)
        NEXT K
     NEXT II
     RETURN

