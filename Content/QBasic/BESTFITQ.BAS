        '**********   PROGRAM BESTFITQ  **********
        '*   BEST FIT PROGRAM FOR QUADRILATERAL  *
        '*    T.R.Chandrupatla & A.D.Belegundu   *
        '*****************************************
        DEFINT I-N: CLS : COLOR 1, 3
        LOCATE 1, 1: PRINT "BESTFIT FOR QUADRILATERAL"; SPACE$(26);
        PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
        VIEW PRINT 2 TO 25: PRINT
        INPUT "Mesh Data File Name <path fn.ext>  ", FILE1$
        OPEN FILE1$ FOR INPUT AS #1
        INPUT "File Name Containing Element Integ_Point Values  ", FILE2$
        INPUT "File Name for Nodal Value Data for CONTOUR ", FILE3$
        LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
        INPUT #1, NN, NE, NM, NDIM, NEN, NDN
        LINE INPUT #1, D$
        INPUT #1, ND, NL, NMPC
        IF NDIM <> 2 OR NEN <> 4 THEN
           PRINT "This program is for 4 noded quadrilaterals only"
           END
        END IF
        DIM X(NN, NDIM), NOC(NE, NEN), V(NE, NEN), F(NN)
'---------------  READ DATA  ----------------
        '----- Coordinates
        LINE INPUT #1, D$
        FOR I = 1 TO NN
           INPUT #1, N
              FOR J = 1 TO NDIM
                 INPUT #1, X(N, J)
              NEXT J
        NEXT I
        '----- Connectivity -----
        LINE INPUT #1, D$
        FOR I = 1 TO NE
          INPUT #1, N
          FOR J = 1 TO NEN
             INPUT #1, NOC(N, J)
          NEXT J
          LINE INPUT #1, D$
        NEXT I
        CLOSE #1
        OPEN FILE2$ FOR INPUT AS #2
        LINE INPUT #2, D$
        '----- Values at Integration Points
        FOR N = 1 TO NE
           FOR J = 1 TO NEN
              INPUT #2, V(N, J)
           NEXT J
        NEXT N
        CLOSE #2
'----- Bandwidth NBW from Connectivity NOC()
     NBW = 0
     FOR I = 1 TO NE
        NMIN = NOC(I, 1): NMAX = NOC(I, 1)
        FOR J = 2 TO NEN
           IF NMIN > NOC(I, J) THEN NMIN = NOC(I, J)
           IF NMAX < NOC(I, J) THEN NMAX = NOC(I, J)
        NEXT J
        NTMP = NMAX - NMIN + 1
        IF NBW < NTMP THEN NBW = NTMP
     NEXT I
     PRINT "The Bandwidth is"; NBW
     NQ = NN
     DIM S(NQ, NBW), SE(4, 4), SH(4, 4), FE(4)
     AL = .5773502692#
'----- Shape Function Values
     SH(1, 1) = .25 * (1 + AL) ^ 2
     SH(1, 2) = .25 * (1 - AL * AL)
     SH(1, 3) = .25 * (1 - AL) ^ 2
     SH(1, 4) = SH(1, 2)
     SH(2, 1) = SH(1, 2): SH(2, 2) = SH(1, 1)
     SH(2, 3) = SH(1, 2): SH(2, 4) = SH(1, 3)
     SH(3, 1) = SH(1, 3): SH(3, 2) = SH(1, 2)
     SH(3, 3) = SH(1, 1): SH(3, 4) = SH(1, 2)
     SH(4, 1) = SH(1, 2): SH(4, 2) = SH(1, 3)
     SH(4, 3) = SH(1, 2): SH(4, 4) = SH(1, 1)
'----- Element Stiffness
     FOR I = 1 TO 4
        FOR J = 1 TO 4
           C = 0
           FOR K = 1 TO 4
              C = C + SH(I, K) * SH(K, J)
           NEXT K
           SE(I, J) = C
        NEXT J
     NEXT I
'----- Stiffness and Loads
     FOR N = 1 TO NE
        FOR I = 1 TO 4
           C = 0
           FOR J = 1 TO 4
              C = C + SH(I, J) * V(N, J)
           NEXT J
           FE(I) = C
        NEXT I
        PRINT ".... Placing in Banded Locations"
        FOR I = 1 TO 4
           NR = NOC(N, I)
              FOR J = 1 TO 4
                 NC = NOC(N, J) - NR + 1
                    IF NC > 0 THEN
                       S(NR, NC) = S(NR, NC) + SE(I, J)
                    END IF
              NEXT J
              F(NR) = F(NR) + FE(I)
        NEXT I
     NEXT N
'----- Equation Solving
     PRINT ".... Solving Equations"
     GOSUB BANSOL
     OPEN FILE3$ FOR OUTPUT AS #3
     PRINT "Nodal Values "
     PRINT #3, "Nodal Values for Data in Files "; FILE1$; " & "; FILE2$
     FOR I = 1 TO NN
        PRINT F(I)
        PRINT #3, F(I)
     NEXT I
     CLOSE #3
     PRINT "Nodal Values are in File "; FILE3$
     END
BANSOL:
'-----  Band Solver  -----
     N1 = NQ - 1
     PRINT "--- Forward Elimination in Progress ---"
     FOR K = 1 TO N1
        NK = NQ - K + 1
        IF NK > NBW THEN NK = NBW
        FOR I = 2 TO NK
           C1 = S(K, I) / S(K, 1)
           I1 = K + I - 1
           FOR J = I TO NK
              J1 = J - I + 1
              S(I1, J1) = S(I1, J1) - C1 * S(K, J)
           NEXT J
           F(I1) = F(I1) - C1 * F(K)
        NEXT I
     NEXT K
     F(NQ) = F(NQ) / S(NQ, 1)
     PRINT "--- Back Substitution in Progress ---"
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

