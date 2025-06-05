DECLARE SUB GlobalNode ()
DECLARE SUB CoordConnect ()
DECLARE SUB SideDiv (I1%, I2%, IDIV%)
DECLARE SUB InputData ()
DECLARE SUB BlockXY (KW%, KSW%)
DECLARE SUB Shape (XI!, ETA!)
DECLARE SUB FileOut ()
'************************************************
'*               PROGRAM MESHGEN                *
'*  MESH GENERATOR FOR TWO DIMENSIONAL REGIONS  *
'*     (c) T.R.CHANDRUPATLA & A.D.BELEGUNDU     *
'************************************************
DEFINT I-N
DEFSNG A-H, O-Z
COMMON SHARED NS, NW, NSJ, NSR, NWR, NNS, NNW, NNT, NGN, NODE
COMMON SHARED NN, NE, NM, NEN
COMMON SHARED IDBLK(), NSD(), NWD(), NGCN(), SR(), WR(), SH()
COMMON SHARED X(), XB(), XP(), NOC(), MAT(), MERG(), NNAR()
COMMON SHARED Title AS STRING, File1 AS STRING, File2 AS STRING
COMMON SHARED Dummy AS STRING
     CLS : COLOR 1, 3: LOCATE 1, 1
     PRINT " 2D - MESH GENERATOR"; SPACE$(31);
     PRINT "(C)  Chandrupatla & Belegundu": COLOR 7, 0
     CALL InputData
     CALL GlobalNode
     CALL CoordConnect
     CALL FileOut
END

SUB BlockXY (KW, KSW)
     '======  Coordinates of 8-Nodes of the Block  ======
     N1 = KSW + KW - 1
     XP(1, 1) = XB(N1, 1): XP(1, 2) = XB(N1, 2)
     XP(3, 1) = XB(N1 + 1, 1): XP(3, 2) = XB(N1 + 1, 2)
     XP(5, 1) = XB(N1 + NS + 2, 1): XP(5, 2) = XB(N1 + NS + 2, 2)
     XP(7, 1) = XB(N1 + NS + 1, 1): XP(7, 2) = XB(N1 + NS + 1, 2)
     XP(2, 1) = SR(KSW, 1): XP(2, 2) = SR(KSW, 2)
     XP(6, 1) = SR(KSW + NS, 1): XP(6, 2) = SR(KSW + NS, 2)
     XP(8, 1) = WR(N1, 1): XP(8, 2) = WR(N1, 2)
     XP(4, 1) = WR(N1 + 1, 1): XP(4, 2) = WR(N1 + 1, 2)
END SUB

SUB CoordConnect
     '------------  Nodal Coordinates  ---------------
     NN = NODE: NELM = 0
     DIM X(NN, 2), XP(8, 2), NOC(2 * NNT, NEN), MAT(2 * NNT)
     FOR KW = 1 TO NW
        FOR KS = 1 TO NS
        KSW = NS * (KW - 1) + KS
        IF IDBLK(KSW) <> 0 THEN
           '---------  Extraction of Block Data  ----------
           NODW = NGCN(KSW + KW - 1) - NNS - 1
           FOR JW = 1 TO NWD(KW) + 1
              ETA = -1 + 2 * (JW - 1) / NWD(KW)
              NODW = NODW + NNS: NODS = NODW
              FOR JS = 1 TO NSD(KS) + 1
                 XI = -1 + 2 * (JS - 1) / NSD(KS)
                 NODS = NODS + 1: NODE = NNAR(NODS)
                 CALL BlockXY(KW, KSW)
                 CALL Shape(XI, ETA)
                 FOR J = 1 TO 2
                    C1 = 0
                    FOR I = 1 TO 8
                       C1 = C1 + SH(I) * XP(I, J)
                    NEXT I
                    X(NODE, J) = C1
                 NEXT J
                 '-----------------  Connectivity  ----------------
                 IF JS <> NSD(KS) + 1 AND JW <> NWD(KW) + 1 THEN
                    N1 = NODE: N2 = NNAR(NODS + 1)
                    N4 = NNAR(NODS + NNS): N3 = NNAR(NODS + NNS + 1)
                    NELM = NELM + 1
                    IF NEN = 3 THEN
                       '------------- Triangular Elements ------------
                       NOC(NELM, 1) = N1: NOC(NELM, 2) = N2
                       NOC(NELM, 3) = N3: MAT(NELM) = IDBLK(KSW)
                       NELM = NELM + 1: NOC(NELM, 1) = N3: NOC(NELM, 2) = N4
                       NOC(NELM, 3) = N1: MAT(NELM) = IDBLK(KSW)
                    ELSE
                       '------------- Quadrilateral Elements ----------
                       NOC(NELM, 1) = N1: NOC(NELM, 2) = N2: MAT(NELM) = IDBLK(KSW)
                       NOC(NELM, 3) = N3: NOC(NELM, 4) = N4
                    END IF
                 END IF
              NEXT JS
           NEXT JW
        END IF
     NEXT KS
   NEXT KW
     NE = NELM
     IF NEN = 3 THEN
     '--------- Readjustment for Triangle Connectivity ----------
        NE2 = NE / 2
        FOR I = 1 TO NE2
           I1 = 2 * I - 1: N1 = NOC(I1, 1): N2 = NOC(I1, 2)
           N3 = NOC(I1, 3): N4 = NOC(2 * I, 2)
           X13 = X(N1, 1) - X(N3, 1): Y13 = X(N1, 2) - X(N3, 2)
           X24 = X(N2, 1) - X(N4, 1): Y24 = X(N2, 2) - X(N4, 2)
           IF (X13 * X13 + Y13 * Y13) > 1.1 * (X24 * X24 + Y24 * Y24) THEN
              NOC(I1, 3) = N4: NOC(2 * I, 3) = N2
           END IF
        NEXT I
     END IF
END SUB

SUB FileOut
     '===== Print Displacements, Stresses, and Reactions
     INPUT "Output File d:\dir\fileName.ext ", File2
     OPEN File2 FOR OUTPUT AS #2
     PRINT #2, "Program MESHGEN - CHANDRUPATLA & BELEGUNDU"
     PRINT #2, Title
     NDIM = 2: NDN = 2
     PRINT #2, "NN  NE  NM  NDIM  NEN  NDN"
     PRINT #2, NN; NE; NM; NDIM; NEN; NDN
     PRINT #2, "ND   NL   NMPC"
     PRINT #2, ND; NL; NMPC
     PRINT #2, "Node#   X    Y"
     FOR I = 1 TO NN
        PRINT #2, I;
        FOR J = 1 TO NDIM
           PRINT #2, X(I, J);
        NEXT J
        PRINT #2,
     NEXT I
     PRINT #2, "Elem#  Node1  Node2  Node3";
     IF NEN = 3 THEN PRINT #2, "  Material#"
     IF NEN = 4 THEN PRINT #2, "  Node4  Material#"
     FOR I = 1 TO NE
        PRINT #2, I;
        FOR J = 1 TO NEN
           PRINT #2, NOC(I, J);
        NEXT J
        PRINT #2, MAT(I)
     NEXT I
     CLOSE #2
     PRINT "Data has been stored in the file "; File2
END SUB

SUB GlobalNode
     '------- Global Node Locations of Corner Nodes ---------
     NTMPI = 1
     FOR I = 1 TO NW + 1
        IF I = 1 THEN IINC = 0 ELSE IINC = NNS * NWD(I - 1)
        NTMPI = NTMPI + IINC: NTMPJ = 0
        FOR J = 1 TO NS + 1
           IJ = (NS + 1) * (I - 1) + J
           IF J = 1 THEN JINC = 0 ELSE JINC = NSD(J - 1)
           NTMPJ = NTMPJ + JINC
           NGCN(IJ) = NTMPI + NTMPJ
        NEXT J
     NEXT I
     '---------------- Node Point Array --------------------
     NNT = NNS * NNW
     DIM NNAR(NNT)
     FOR I = 1 TO NNT
        NNAR(I) = -1
     NEXT I
     '--------- Zero Non-Existing Node Locations ---------
     FOR KW = 1 TO NW
        FOR KS = 1 TO NS
           KSW = NS * (KW - 1) + KS
           IF IDBLK(KSW) <= 0 THEN
              '-------- Operation within an Empty Block --------
              K1 = (KW - 1) * (NS + 1) + KS: N1 = NGCN(K1)
              NS1 = 2: IF KS = 1 THEN NS1 = 1
              NW1 = 2: IF KW = 1 THEN NW1 = 1
              NS2 = NSD(KS) + 1
              IF KS < NS THEN
                 IF IDBLK(KSW + 1) > 0 THEN NS2 = NSD(KS)
              END IF
              NW2 = NWD(KW) + 1
              IF KW < NW THEN
                 IF IDBLK(KSW + NS) > 0 THEN NW2 = NWD(KW)
              END IF
              FOR I = NW1 TO NW2
                 IN1 = N1 + (I - 1) * NNS
                 FOR J = NS1 TO NS2
                    IJ = IN1 + J - 1
                    NNAR(IJ) = 0
                 NEXT J
              NEXT I
              ICT = 0
              IF NS2 = NSD(KS) OR NW2 = NWD(KW) THEN ICT = 1
              IF KS = NS OR KW = NW THEN ICT = 1
              IF ICT = 0 THEN
                 IF IDBLK(KSW + NS + 1) > 0 THEN NNAR(IJ) = -1
              END IF
           END IF
        NEXT KS
     NEXT KW
     '--------  Node Identification for Side Merging ------
     IF NSJ > 0 THEN
        FOR I = 1 TO NSJ
           I1 = MERG(I, 1): I2 = MERG(I, 2)
           CALL SideDiv(I1, I2, IDIV)
           IA1 = NGCN(I1): IA2 = NGCN(I2)
           IASTP = (IA2 - IA1) / IDIV
           I1 = MERG(I, 3): I2 = MERG(I, 4)
           CALL SideDiv(I1, I2, IDIV)
           IB1 = NGCN(I1): IB2 = NGCN(I2)
           IBSTP = (IB2 - IB1) / IDIV
           IAA = IA1 - IASTP
           FOR IBB = IB1 TO IB2 STEP IBSTP
              IAA = IAA + IASTP
              IF IBB = IAA THEN NNAR(IAA) = -1 ELSE NNAR(IBB) = IAA
           NEXT IBB
        NEXT I
     END IF
     '----------  Final Node Numbers in the Array  --------
     NODE = 0
     FOR I = 1 TO NNT
        IF NNAR(I) > 0 THEN
           II = NNAR(I): NNAR(I) = NNAR(II)
        ELSEIF NNAR(I) < 0 THEN
           NODE = NODE + 1: NNAR(I) = NODE
        END IF
     NEXT I
END SUB

SUB InputData
     INPUT "Input File d:\dir\fileName.ext ", File1
     OPEN File1 FOR INPUT AS #1
     '=============  READ DATA  ===============
     LINE INPUT #1, Dummy: LINE INPUT #1, Title
     LINE INPUT #1, Dummy
     INPUT #1, NEN          ' NEN = 3 for Triangle  4 for Quad
     IF NEN < 3 THEN NEN = 3
     IF NEN > 4 THEN NEN = 4
     'Hints:  A region is divided into 4-cornered blocks viewed as a
     '        mapping from a Checkerboard pattern of S- and W- Sides
     '        S- Side is one with lower number of final divisions
     '        Blocks, Corners, S- and W- Sides are labeled as shown in Fig. 12.2
     '        Make a sketch and identify void blocks and merging sides
     '----- Block Data -----
     '#S-Spans(NS)  #W-Spans(NW)  #PairsOfEdgesMerged(NSJ)
     LINE INPUT #1, Dummy: LINE INPUT #1, Dummy
     INPUT #1, NS, NW, NSJ
     NSW = NS * NW: NGN = (NS + 1) * (NW + 1): NM = 1
     DIM IDBLK(NSW), NSD(NS), NWD(NW), NGCN(NGN), SH(8)
     '------------- Span Divisions ---------------
     LINE INPUT #1, Dummy
     NNS = 1: NNW = 1
     '--- Number of divisions for each S-Span
     LINE INPUT #1, Dummy
     FOR KS = 1 TO NS
        INPUT #1, N
        INPUT #1, NSD(N)
        NNS = NNS + NSD(N)
     NEXT KS
     '--- Number of divisions for each W-Span
     LINE INPUT #1, Dummy
     FOR KW = 1 TO NW
        INPUT #1, N
        INPUT #1, NWD(N)
        NNW = NNW + NWD(N)
     NEXT KW
     '--- Block Material Data
     INPUT #1, Dummy: INPUT #1, Dummy
     '-------- Block Identifier / Material# (Default# is 1) --------
     FOR I = 1 TO NSW: IDBLK(I) = 1: NEXT I
     DO
        INPUT #1, NTMP
        IF NTMP = 0 THEN EXIT DO
        INPUT #1, IDBLK(NTMP)
        IF NM < IDBLK(NTMP) THEN NM = IDBLK(NTMP)
     LOOP
     '-----------------   Block Corner Data  ---------------
     NSR = NS * (NW + 1): NWR = NW * (NS + 1)
     DIM XB(NGN, 2), SR(NSR, 2), WR(NWR, 2)
     INPUT #1, Dummy: INPUT #1, Dummy
     DO
        INPUT #1, NTMP
        IF NTMP = 0 THEN EXIT DO
        INPUT #1, XB(NTMP, 1)
        INPUT #1, XB(NTMP, 2)
     LOOP
     '---------- Evaluate Mid-points of S-Sides -------------
     FOR I = 1 TO NW + 1
        FOR J = 1 TO NS
           IJ = (I - 1) * NS + J
           SR(IJ, 1) = .5 * (XB(IJ + I - 1, 1) + XB(IJ + I, 1))
           SR(IJ, 2) = .5 * (XB(IJ + I - 1, 2) + XB(IJ + I, 2))
        NEXT J
     NEXT I
     '---------- Evaluate Mid-points of W-Sides -------------
     FOR I = 1 TO NW
        FOR J = 1 TO NS + 1
           IJ = (I - 1) * (NS + 1) + J
           WR(IJ, 1) = .5 * (XB(IJ, 1) + XB(IJ + NS + 1, 1))
           WR(IJ, 2) = .5 * (XB(IJ, 2) + XB(IJ + NS + 1, 2))
        NEXT J
     NEXT I
     '------ Mid Points for Sides that are curved or graded ------
     LINE INPUT #1, Dummy: LINE INPUT #1, Dummy
     '--- S-Sides
     DO
        INPUT #1, NTMP
        IF NTMP = 0 THEN EXIT DO
        INPUT #1, SR(NTMP, 1)
        INPUT #1, SR(NTMP, 2)
     LOOP
     LINE INPUT #1, Dummy
     '--- W-Sides
     DO
        INPUT #1, NTMP
        IF NTMP = 0 THEN EXIT DO
        INPUT #1, WR(NTMP, 1)
        INPUT #1, WR(NTMP, 2)
     LOOP
     '--------- Merging Sides ----------
     IF NSJ > 0 THEN
        INPUT #1, Dummy: INPUT #1, Dummy
        DIM MERG(NSJ, 4)
        FOR I = 1 TO NSJ
           INPUT #1, N
           INPUT #1, L1
           INPUT #1, L2
           CALL SideDiv(L1, L2, IDIV1)
           INPUT #1, L3
           INPUT #1, L4
           CALL SideDiv(L3, L4, IDIV2)
           IF IDIV1 <> IDIV2 THEN
              PRINT "#Div don't match. Check merge data."
              END
           END IF
           MERG(I, 1) = L1: MERG(I, 2) = L2
           MERG(I, 3) = L3: MERG(I, 4) = L4
        NEXT I
     END IF
     CLOSE #1
END SUB

SUB Shape (XI, ETA)
     '==============  Shape Functions  ================
     SH(1) = -(1 - XI) * (1 - ETA) * (1 + XI + ETA) / 4
     SH(2) = (1 - XI * XI) * (1 - ETA) / 2
     SH(3) = -(1 + XI) * (1 - ETA) * (1 - XI + ETA) / 4
     SH(4) = (1 - ETA * ETA) * (1 + XI) / 2
     SH(5) = -(1 + XI) * (1 + ETA) * (1 - XI - ETA) / 4
     SH(6) = (1 - XI * XI) * (1 + ETA) / 2
     SH(7) = -(1 - XI) * (1 + ETA) * (1 + XI - ETA) / 4
     SH(8) = (1 - ETA * ETA) * (1 - XI) / 2
END SUB

SUB SideDiv (I1, I2, IDIV)
     '===========  Number of Divisions  for Side I1,I2  ===========
     IMIN = I1: IMAX = I2
     IF IMIN > I2 THEN
        IMIN = I2
        IMAX = I1
     END IF
     IF (IMAX - IMIN) = 1 THEN
        IDIV = NGCN(IMAX) - NGCN(IMIN)
     ELSE
        IDIV = (NGCN(IMAX) - NGCN(IMIN)) / NNS
     END IF
END SUB

