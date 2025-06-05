     '*****      PROGRAM BESTFIT        *****
     '*          BEST FIT PROGRAM           *
     '*        FOR 3-NODED TRIANGLES        *
     '* T.R.Chandrupatla and A.D.Belegundu  *
     '***************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "BESTFIT FOR TRIANGULAR ELEMENTS"; SPACE$(20);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Mesh Data File Name <path fn.ext>  ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     INPUT "File Name Containing Element Values  ", FILE3$
     OPEN FILE3$ FOR INPUT AS #3
     INPUT "Output Data File Name <path fn.ext>  ", FILE2$
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     IF NDIM <> 2 OR NEN <> 3 THEN
        PRINT "This program supports 3 Noded Triangles only"
        END
     END IF
     DIM X(NN, 2), NOC(NE, 3), FS(NE)
     '----- Coordinates
     LINE INPUT #1, D$
     FOR I = 1 TO NN: INPUT #1, N: FOR J = 1 TO 2
     INPUT #1, X(N, J): NEXT J: NEXT I
     '----- Connectivity
     LINE INPUT #1, D$
     FOR I = 1 TO NE: INPUT #1, N: FOR J = 1 TO NEN
     INPUT #1, NOC(N, J): NEXT J
     LINE INPUT #1, D$
     NEXT I
     '--- Element Function Values
     LINE INPUT #3, D$
     FOR I = 1 TO NE
     INPUT #3, FS(I): NEXT I
     CLOSE #1, #3
     '--- Bandwidth NBW from Connectivity NOC()
     NQ = NN: NBW = 0
     FOR I = 1 TO NE
        NMIN = NOC(I, 1): NMAX = NOC(I, 1)
        FOR J = 2 TO 3
           IF NMIN > NOC(I, J) THEN NMIN = NOC(I, J)
           IF NMAX < NOC(I, J) THEN NMAX = NOC(I, J)
        NEXT J
        NTMP = NMAX - NMIN + 1
        IF NBW < NTMP THEN NBW = NTMP
     NEXT I
     PRINT "Bandwidth is "; NBW
     DIM S(NQ, NBW), F(NQ)
     '---  Global Stiffness Matrix
     FOR N = 1 TO NE
        PRINT "Forming Stiffness Matrix of Element "; N
        GOSUB ELSTIF
        PRINT "... Placing in Global Locations"
        FOR II = 1 TO 3
           NR = NOC(N, II): F(NR) = F(NR) + FE(II)
           FOR JJ = 1 TO 3
              NC = NOC(N, JJ) - NR + 1
              IF NC > 0 THEN
                 S(NR, NC) = S(NR, NC) + SE(II, JJ)
              END IF
           NEXT JJ
        NEXT II
     NEXT N
     '--- Equation Solving
     GOSUB BANSOL
     OPEN FILE2$ FOR OUTPUT AS #2
     PRINT "Nodal Values "
     PRINT #2, "Nodal Values for Data in Files "; FILE1$; " & "; FILE3$
     FOR I = 1 TO NN
        PRINT #2, F(I)
     NEXT I
     CLOSE #2
     PRINT "Nodal Value Data is in File "; FILE2$
     PRINT "You may now run the CONTOURA or CONTOURB program"
     END
ELSTIF:
     '--- Element Stiffness Formation
     I1 = NOC(N, 1): I2 = NOC(N, 2): I3 = NOC(N, 3)
     X1 = X(I1, 1): Y1 = X(I1, 2)
     X2 = X(I2, 1): Y2 = X(I2, 2)
     X3 = X(I3, 1): Y3 = X(I3, 2)
     X21 = X2 - X1: X32 = X3 - X2: X13 = X1 - X3
     Y12 = Y1 - Y2: Y23 = Y2 - Y3: Y31 = Y3 - Y1
     DJ = X13 * Y23 - X32 * Y31      'DETERMINANT OF JACOBIAN
     AE = ABS(DJ) / 24
     SE(1, 1) = 2 * AE: SE(1, 2) = AE: SE(1, 3) = AE
     SE(2, 1) = AE: SE(2, 2) = 2 * AE: SE(2, 3) = AE
     SE(3, 1) = AE: SE(3, 2) = AE: SE(3, 3) = 2 * AE
     A1 = FS(N) * ABS(DJ) / 6
     FE(1) = A1: FE(2) = A1: FE(3) = A1
     RETURN
BANSOL:
     N1 = NQ - 1
     '----- Forward Elimination
     FOR K = 1 TO N1
     NK = NQ - K + 1
     IF NK > NBW THEN NK = NBW
     FOR I = 2 TO NK
     C1 = S(K, I) / S(K, 1)
     I1 = K + I - 1
     FOR J = I TO NK
     J1 = J - I + 1
     S(I1, J1) = S(I1, J1) - C1 * S(K, J): NEXT J
     F(I1) = F(I1) - C1 * F(K): NEXT I: NEXT K
     '----- Back Substitution
     F(NQ) = F(NQ) / S(NQ, 1)
     FOR KK = 1 TO N1
     K = NQ - KK
     C1 = 1 / S(K, 1)
     F(K) = C1 * F(K)
     NK = NQ - K + 1
     IF NK > NBW THEN NK = NBW
     FOR J = 2 TO NK
     F(K) = F(K) - C1 * S(K, J) * F(K + J - 1)
     NEXT J
     NEXT KK
     RETURN

