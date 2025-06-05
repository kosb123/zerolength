     '*****       PROGRAM JACOBI            *    
     '*     Generalized Jacobi's Method     *
     '*       for symmetric matrices        *
     '*  T.R.Chandrupatla and A.D.Belegundu *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "GENERALIZED JACOBI'S METHOD"; SPACE$(24);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Name of Input File  ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     INPUT "Name of Output File  ", FILE2$
     OPEN FILE2$ FOR OUTPUT AS #2
     '--- Read in Number of Equations
     LINE INPUT #1, D$
     LINE INPUT #1, D$: INPUT #1, NQ, NBW
     DIM S(NQ, NQ), GM(NQ, NQ), EVL(NQ), EVC(NQ, NQ), NORD(NQ)
     REM NORD( ) is for ascending order of eigenvalues
     FOR I = 1 TO NQ: NORD(I) = I: NEXT I
     TOL = .000001: PI = 3.14159
     PRINT "Default Tolerance is  1E-6"
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
     '----- Initialize Eigenvector Matrix -----
     FOR I = 1 TO NQ: FOR J = 1 TO NQ: EVC(I, J) = 0
     NEXT J: EVC(I, I) = 1: NEXT I
     C1 = S(1, 1): C2 = GM(1, 1)
     FOR I = 1 TO NQ
        IF C1 > S(I, I) THEN C1 = S(I, I)
        IF C2 < GM(I, I) THEN C2 = GM(I, I)
     NEXT I
     TOLS = TOL * C1: TOLM = TOL * C2: NSWMAX = 50
     INPUT "Maximum Number of Sweeps < default = 50 > ", A$
     IF VAL(A$) > 0 THEN NSWMAX = VAL(A$)
     CLS
     K1 = 1: I1 = 1: NSW = 0
   DO
     NSW = NSW + 1
     IF NSW > NSWMAX THEN
        PRINT "No Convergence after "; NSWMAX; " sweeps"
        PRINT #2, "No Convergence after "; NSWMAX; " sweeps"
        CLOSE #2: END
     END IF
     PRINT "------  SWEEP NUMBER    "; NSW
     FOR K = K1 TO NQ - 1
     FOR I = I1 TO K
     J = NQ - K + I
     IF ABS(S(I, J)) > TOLS OR ABS(GM(I, J)) > TOLM THEN
     AA = S(I, I) * GM(I, J) - GM(I, I) * S(I, J)
     BB = S(J, J) * GM(I, J) - GM(J, J) * S(I, J)
     CC = S(I, I) * GM(J, J) - GM(I, I) * S(J, J)
     CAB = .25 * CC * CC + AA * BB
     IF CAB < 0 THEN
         PRINT "Square Root of Negative Term -- Check Matrices"
         PRINT #2, "Square Root of Negative Term -- Check Matrices"
         CLOSE #2: END
     END IF
     IF AA = 0 THEN
        BET = 0: ALP = -S(I, J) / S(I, I)
     ELSEIF BB = 0 THEN
        ALP = 0: BET = -S(I, J) / S(J, J)
     ELSE
        SQC = SQR(CAB): IF CC < 0 THEN SQC = -SQC
        ALP = (-.5 * CC + SQC) / AA
        BET = -AA * ALP / BB
     END IF
     '----- Only Upper Triangular Part is used in Diagonalization
     IF I > 1 THEN
     FOR N = 1 TO I - 1
     SI = S(N, I): SJ = S(N, J): EMI = GM(N, I): EMJ = GM(N, J)
     S(N, I) = SI + BET * SJ: S(N, J) = SJ + ALP * SI
     GM(N, I) = EMI + BET * EMJ: GM(N, J) = EMJ + ALP * EMI
     NEXT N
     END IF
     IF J < NQ THEN
     FOR N = J + 1 TO NQ
     SI = S(I, N): SJ = S(J, N): EMI = GM(I, N): EMJ = GM(J, N)
     S(I, N) = SI + BET * SJ: S(J, N) = SJ + ALP * SI
     GM(I, N) = EMI + BET * EMJ: GM(J, N) = EMJ + ALP * EMI
     NEXT N
     END IF
     IF I < J THEN
     FOR N = I + 1 TO J - 1
     SI = S(I, N): SJ = S(N, J): EMI = GM(I, N): EMJ = GM(N, J)
     S(I, N) = SI + BET * SJ: S(N, J) = SJ + ALP * SI
     GM(I, N) = EMI + BET * EMJ: GM(N, J) = EMJ + ALP * EMI
     NEXT N
     END IF
     SII = S(I, I): SIJ = S(I, J): SJJ = S(J, J)
     S(I, J) = 0: S(I, I) = SII + 2 * BET * SIJ + BET * BET * SJJ
     S(J, J) = SJJ + 2 * ALP * SIJ + ALP * ALP * SII
     EII = GM(I, I): EIJ = GM(I, J): EJJ = GM(J, J)
     GM(I, J) = 0: GM(I, I) = EII + 2 * BET * EIJ + BET * BET * EJJ
     GM(J, J) = EJJ + 2 * ALP * EIJ + ALP * ALP * EII
     REM *** EIGENVECTORS ***
     FOR N = 1 TO NQ
     EVI = EVC(N, I): EVJ = EVC(N, J)
     EVC(N, I) = EVI + BET * EVJ: EVC(N, J) = EVJ + ALP * EVI
     NEXT N
     END IF
     NEXT I: NEXT K
     FOR K = 1 TO NQ - 1
     FOR I = 1 TO K
     J = NQ - K + I
     IFL = 0
     IF ABS(S(I, J)) > TOLS OR ABS(GM(I, J)) > TOLM THEN
        K1 = K: I1 = I: IFL = 1
     END IF
     IF IFL = 1 THEN EXIT FOR
     NEXT I
     IF IFL = 1 THEN EXIT FOR
     NEXT K
   LOOP WHILE IFL = 1
     '-----  Calculation of Eigenvalues -----
     FOR I = 1 TO NQ
        IF ABS(GM(I, I)) < TOLM THEN GM(I, I) = TOLM
        EVL(I) = S(I, I) / GM(I, I)
     NEXT I
     '----- Scaling of Eigenvectors
     FOR I = 1 TO NQ
        GM2 = SQR(ABS(GM(I, I)))
        FOR J = 1 TO NQ
           EVC(J, I) = EVC(J, I) / GM2
        NEXT J
     NEXT I
     '-----   RESULTS   -----
     '--- Ascending Order of Eigenvalues
     FOR I = 1 TO NQ - 1
        II = NORD(I): I1 = II
        C1 = EVL(II): J1 = I
        FOR J = I + 1 TO NQ
           IJ = NORD(J)
           IF C1 > EVL(IJ) THEN
              C1 = EVL(IJ): I1 = IJ: J1 = J
           END IF
        NEXT J
        NORD(I) = I1: NORD(J1) = II
     NEXT I
     PRINT "Eigenvalues and Eigenvectors for Data in File "; FILE1$
     PRINT #2, "Eigenvalues and Eigenvectors for Data in File "; FILE1$
     FOR I = 1 TO NQ
        II = NORD(I)
        PRINT "Eigenvalue Number  "; I
        PRINT #2, "Eigenvalue Number  "; I
        C = EVL(II): OMEGA = SQR(C): FREQ = .5 * OMEGA / PI
        PRINT USING "Eigenvalue = ##.###^^^^"; C;
        PRINT #2, USING "Eigenvalue = ##.####^^^^"; C;
        PRINT USING "   Omega = ##.###^^^^"; OMEGA;
        PRINT USING "   Freq = ##.###^^^^ Hz"; FREQ
        PRINT #2, USING "   Omega = ##.###^^^^"; OMEGA;
        PRINT #2, USING "   Freq = ##.###^^^^ Hz"; FREQ
        PRINT "Eigenvector ": PRINT #2, "Eigenvector "
        FOR J = 1 TO NQ
           PRINT USING "##.###^^^^ "; EVC(J, II);
           PRINT #2, USING "##.###^^^^ "; EVC(J, II);
           IFL = J - 7 * INT(J / 7)
           IF IFL = 0 THEN PRINT : PRINT #2,
        NEXT J
        PRINT : PRINT #2,
     NEXT I
     CLOSE #2
     END

