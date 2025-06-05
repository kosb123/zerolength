     '***************************************
     '*          PROGRAM  HEAT1D           *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "   HEAT1    PROGRAM"; SPACE$(23);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     PRINT : INPUT "File Name for Input Data ", FILE1$
     INPUT "File Name for Output Data ", FILE2$
     OPEN FILE1$ FOR INPUT AS #1
     LINE INPUT #1, D$: LINE INPUT #1, TITLE$
     LINE INPUT #1, D$
     INPUT #1, NE, NBC, NQ
     NN = NE + 1
     NBW = 2    'NBW IS THE HALF-BAND-WIDTH
     DIM X(NN), S(NN, NBW), BC$(NBC), TC(NE), F(NN), V(NBC), H(NBC)
     DIM NB(NBC)
     LINE INPUT #1, D$
     FOR I = 1 TO NE
     INPUT #1, N, TC(N): NEXT I
     LINE INPUT #1, D$
     FOR I = 1 TO NN
     INPUT #1, N, X(N): NEXT I
     LINE INPUT #1, D$
     FOR I = 1 TO NBC
     INPUT #1, NB(I), BC$(I)
     IF BC$(I) = "TEMP" OR BC$(I) = "temp" THEN INPUT #1, V(I)
     IF BC$(I) = "FLUX" OR BC$(I) = "flux" THEN INPUT #1, V(I)
     IF BC$(I) = "CONV" OR BC$(I) = "conv" THEN INPUT #1, H(I), V(I)
     NEXT I
     '--- CALCULATE AND INPUT NODAL HEAT SOURCE VECTOR ---
     FOR I = 1 TO NN: F(I) = 0: NEXT I
     LINE INPUT #1, D$
     IF NQ > 0 THEN
       FOR I = 1 TO NQ
       INPUT #1, N, F(N): NEXT I
     END IF
     CLOSE #1
     '--- STIFFNESS MATRIX ---
     FOR J = 1 TO NBW
     FOR I = 1 TO NN: S(I, J) = 0: NEXT I: NEXT J
     FOR I = 1 TO NE
     I1 = I: I2 = I + 1
     ELL = ABS(X(I2) - X(I1))
     EKL = TC(I) / ELL
     S(I1, 1) = S(I1, 1) + EKL
     S(I2, 1) = S(I2, 1) + EKL
     S(I1, 2) = S(I1, 2) - EKL: NEXT I
     '--- ACCOUNT FOR B.C.'S ---
     AMAX = 0
     FOR I = 1 TO NN
     IF S(I, 1) > AMAX THEN AMAX = S(I, 1)
     NEXT I
     CNST = AMAX * 10000
     FOR I = 1 TO NBC
     N = NB(I)
     IF BC$(I) = "CONV" OR BC$(I) = "conv" THEN
       S(N, 1) = S(N, 1) + H(I)
       F(N) = F(N) + H(I) * V(I)
     ELSEIF BC$(I) = "FLUX" OR BC$(I) = "flux" THEN
       F(N) = F(N) - V(I)
     ELSE
       S(N, 1) = S(N, 1) + CNST
       F(N) = F(N) + CNST * V(I)
     END IF
     NEXT I
     GOSUB BANSOL
     '--- F CONTAINS THE SOLUTION. 'RHS' IS OVER-WRITTEN
     OPEN FILE2$ FOR OUTPUT AS #1
     PRINT TITLE$
     PRINT #1, TITLE$
     PRINT "NODE#  TEMPERATURE"
     PRINT #1, "NODE#  TEMPERATURE"
     FOR I = 1 TO NN
       PRINT I; "  "; F(I)
       PRINT #1, I; "  "; F(I)
     NEXT I
     CLOSE #1
     PRINT "Output is in file "; FILE2$
     END
BANSOL:
     '--- EQUATION SOLVING ... F CONTAINS THE TEMPERATURES ON EXIT
     N = NN
     FOR K = 1 TO N - 1
     NBK = N - K + 1
     IF N - K + 1 > NBW THEN NBK = NBW
     FOR I = K + 1 TO NBK + K - 1
     I1 = I - K + 1
     C = S(K, I1) / S(K, 1)
     FOR J = I TO NBK + K - 1
     J1 = J - I + 1
     J2 = J - K + 1
     S(I, J1) = S(I, J1) - C * S(K, J2): NEXT J
     F(I) = F(I) - C * F(K)
     NEXT I: NEXT K
     '--- BACK SUBSTITUTION ---
     F(N) = F(N) / S(N, 1)
     FOR II = 1 TO N - 1
     I = N - II
     NBI = N - I + 1
     IF N - I + 1 > NBW THEN NBI = NBW
     SUM = 0
     FOR J = 2 TO NBI
     SUM = SUM + S(I, J) * F(I + J - 1): NEXT J
     F(I) = (F(I) - SUM) / S(I, 1): NEXT II
     RETURN

