     '*****       PROGRAM HEXAFRON          *****
     '*   3-D STRESS ANALYSIS USING  8-NODE     *
     '*    ISOPARAMETRIC HEXAHEDRAL ELEMENT     *
     '*          USING FRONTAL SOLVER           *
     '*    T.R.Chandrupatla and A.DBelegundu    *
     '*******************************************
     DEFDBL A-H, O-Z: DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "3-D HEXAHEDRAL ELEMENT"; SPACE$(29);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25
     INPUT "Name of Data File  < d:fn.ext > ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     INPUT "Name of File for Output < d:fn.ext > ", FILE2$
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     NPR = 3   'Material properties E, Nu, Alpha
     '----- Total dof is  NQ
     NQ = NDN * NN
     DIM X(NN, NDIM), NOC(NE, NEN), NU(ND), U(ND), MAT(NE), F(NQ)
     DIM SE(24, 24), PM(NM, NPR), MPC(NMPC, 2), BT(NMPC, 3)
     DIM XI(3, 8), XNI(3, 8), D(6, 6), GN(3, 8), H(9, 24), TJ(3, 3)
     DIM AJ(3, 3), G(6, 9), B(6, 24), DB(6, 24), QT(24), STR(6), DT(NE)
     GOSUB GETDATA
     NEDF = NEN * NDN
     GOSUB PREFRONT
     OPEN "R", 3, "SCRATCH.DAT", 8      'Scratch file for writing
     FIELD 3, 4 AS VAR$, 4 AS COEF$     'Field definition
     ICOUNT = 0
     '=====  FRONTAL ASSEMBLY & ELIMINATON ETC.  =====
     '----- Corner Nodes and Integration Points
     GOSUB INTEG: MTN1 = 0
     FOR N = 1 TO NE
        PRINT "... Forming Stiffness Matrix of Element "; N
        MTN = MAT(N): IF MTN <> MTN1 THEN GOSUB DMAT
        GOSUB ELSTIF
        IF N = 1 THEN
           CNST = 0
           FOR I = 1 TO NEDF: CNST = CNST + SE(I, I): NEXT I
           CNST = 1E+11 * CNST
           GOSUB MPCFRON
        END IF
        '----- Account for temperature loads QT()
        FOR I = 1 TO NEN
          IL = 3 * (I - 1): IG = 3 * (ABS(NOC(N, I)) - 1)
          FOR J = 1 TO 3
             IL = IL + 1: IG = IG + 1
             F(IG) = F(IG) + QT(IL)
          NEXT J
        NEXT I
        GOSUB FRONT          'Frontal assembly  and Forward Elimination
     NEXT N
     '----- Assembly and reduction are complete
     '----- Now Backsubstitute
     GOSUB BACKSUB
     CLOSE #3
     KILL "SCRATCH.DAT"
     OPEN FILE2$ FOR OUTPUT AS #2
     PRINT "Node#   X-Displ.    Y-Displ.    Z-Displ."
     PRINT #2, "Node#   X-Displ.    Y-Displ.    Z-Displ."
     FOR I = 1 TO NN
        PRINT USING " ###"; I;
        PRINT #2, USING " ###"; I;
        II = 3 * (I - 1)
        FOR J = 1 TO 3
           PRINT USING " ##.####^^^^"; F(II + J);
           PRINT #2, USING " ##.####^^^^"; F(II + J);
        NEXT J
        PRINT : PRINT #2,
     NEXT I
'----- Reaction Calculation
     PRINT "DOF#        Reaction"
     PRINT #2, "DOF#        Reaction"
     FFF1$ = "  ###     ##.####^^^^"
     FOR I = 1 TO ND
        N = NU(I)
        R = CNST * (U(I) - F(N))
        PRINT USING FFF1$; N; R
        PRINT #2, USING FFF1$; N; R
     NEXT I
'-----  Stress Calculations -----
     MTN1 = 0
     FOR N = 1 TO NE
        PRINT #2, "Von Mises Stress at 8 Integation Pts. in Elem#  "; N
        MTN = MAT(N): IF MTN <> MTN1 THEN GOSUB DMAT
        CAL = AL * DT(N)
        FOR IP = 1 TO 8
           '--- Von Mises Stress at Integration Points
           GOSUB DBMAT
                 '--- Element Nodal Displacements stored in QT()
           FOR I = 1 TO 8
              IN = 3 * (ABS(NOC(N, I)) - 1)
              II = 3 * (I - 1)
              FOR J = 1 TO 3
                 QT(II + J) = F(IN + J)
              NEXT J
           NEXT I
           '--- Stress Calculation STR = DB * Q
           FOR I = 1 TO 6
              STR(I) = 0
              FOR J = 1 TO 24
                 STR(I) = STR(I) + DB(I, J) * QT(J)
              NEXT J
              STR(I) = STR(I) - CAL * (D(I, 1) + D(I, 2) + D(I, 3))
           NEXT I
           '--- Calculation of Von Mises Stress at IP
           SIV1 = STR(1) + STR(2) + STR(3)
           SIV2 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1)
           SIV2 = SIV2 - STR(4) ^ 2 - STR(5) ^ 2 - STR(6) ^ 2
           VM = SQR(SIV1 * SIV1 - 3 * SIV2)
           IF IP = 5 THEN PRINT #2,
           PRINT #2, USING "   ##.####^^^^"; VM;
        NEXT IP
        PRINT #2,
     NEXT N
     PRINT "The Results are saved in the file  "; FILE2$
     CLOSE #2
     END
GETDATA:
     '===============  READ DATA  ====================
     '----- Coordinates
     LINE INPUT #1, D$
     FOR I = 1 TO NN
        INPUT #1, N
        FOR J = 1 TO NDIM
           INPUT #1, X(N, J)
        NEXT J
     NEXT I
     '----- Connectivity, Material, Temp-change
     LINE INPUT #1, D$
     FOR I = 1 TO NE
        INPUT #1, N
        FOR J = 1 TO NEN
           INPUT #1, NOC(N, J)
        NEXT J
        INPUT #1, MAT(N), DT(N)
     NEXT I
     '----- Displacement BC
     LINE INPUT #1, D$
     FOR I = 1 TO ND: INPUT #1, NU(I), U(I): NEXT I
     '----- Component Loads
     LINE INPUT #1, D$
     FOR I = 1 TO NL: INPUT #1, N, F(N): NEXT I
     '----- Material Properties
     LINE INPUT #1, D$
     FOR I = 1 TO NM
        INPUT #1, N
        FOR J = 1 TO NPR
           INPUT #1, PM(N, J)
        NEXT J
     NEXT I
     IF NMPC > 0 THEN
        '-----  Multi-point Constraints
        LINE INPUT #1, D$
        FOR I = 1 TO NMPC
           INPUT #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        NEXT I
     END IF
     CLOSE #1
     RETURN
INTEG:
'------- Integration Points XNI() --------
     C = .57735026919#
     XI(1, 1) = -1: XI(2, 1) = -1: XI(3, 1) = -1
     XI(1, 2) = 1: XI(2, 2) = -1: XI(3, 2) = -1
     XI(1, 3) = 1: XI(2, 3) = 1: XI(3, 3) = -1
     XI(1, 4) = -1: XI(2, 4) = 1: XI(3, 4) = -1
     XI(1, 5) = -1: XI(2, 5) = -1: XI(3, 5) = 1
     XI(1, 6) = 1: XI(2, 6) = -1: XI(3, 6) = 1
     XI(1, 7) = 1: XI(2, 7) = 1: XI(3, 7) = 1
     XI(1, 8) = -1: XI(2, 8) = 1: XI(3, 8) = 1
     FOR I = 1 TO 8
        XNI(1, I) = C * XI(1, I)
        XNI(2, I) = C * XI(2, I)
        XNI(3, I) = C * XI(3, I)
     NEXT I
     RETURN
DMAT:
     '--- D() Matrix relating Stresses to Strains
     E = PM(MTN, 1): PNU = PM(MTN, 2): AL = PM(MTN, 3)
     C1 = E / ((1 + PNU) * (1 - 2 * PNU))
     C2 = .5 * E / (1 + PNU)
     FOR I = 1 TO 6: FOR J = 1 TO 6: D(I, J) = 0: NEXT J: NEXT I
     D(1, 1) = C1 * (1 - PNU): D(1, 2) = C1 * PNU: D(1, 3) = D(1, 2)
     D(2, 1) = D(1, 2): D(2, 2) = D(1, 1): D(2, 3) = D(1, 2)
     D(3, 1) = D(1, 3): D(3, 2) = D(2, 3): D(3, 3) = D(1, 1)
     D(4, 4) = C2: D(5, 5) = C2: D(6, 6) = C2
     MTN1 = MTN
     RETURN
ELSTIF:
'--------  Element Stiffness  -----
     FOR I = 1 TO 24: FOR J = 1 TO 24
     SE(I, J) = 0: NEXT J: QT(I) = 0: NEXT I
     DTE = DT(N)
     '--- Weight Factor is ONE
     '--- Loop on Integration Points
     FOR IP = 1 TO 8
        PRINT "Integration Point = "; IP
        GOSUB DBMAT
        '--- Element Stiffness Matrix  SE
        FOR I = 1 TO 24
           FOR J = 1 TO 24
              FOR K = 1 TO 6
                 SE(I, J) = SE(I, J) + B(K, I) * DB(K, J) * DJ
              NEXT K
           NEXT J
        NEXT I
        '--- Determine Temperature Load QT()
        C = AL * DTE
        FOR I = 1 TO 24
           DSUM = DB(1, I) + DB(2, I) + DB(3, I)
           QT(I) = QT(I) + C * ABS(DJ) * DSUM / 6
        NEXT I
     NEXT IP
     RETURN
DBMAT:
'-------  DB()  MATRIX  ------
     '--- Gradient of Shape Functions - The GN() Matrix
     FOR I = 1 TO 3
        FOR J = 1 TO 8
           C = 1
           FOR K = 1 TO 3
              IF K <> I THEN
                 C = C * (1 + XI(K, J) * XNI(K, IP))
              END IF
           NEXT K
           GN(I, J) = .125 * XI(I, J) * C
        NEXT J
     NEXT I
     '--- Formation of Jacobian  TJ
     FOR I = 1 TO 3
        FOR J = 1 TO 3
           TJ(I, J) = 0
           FOR K = 1 TO 8
              KN = ABS(NOC(N, K))
              TJ(I, J) = TJ(I, J) + GN(I, K) * X(KN, J)
           NEXT K
        NEXT J
     NEXT I
     '--- Determinant of the JACOBIAN
     DJ1 = TJ(1, 1) * (TJ(2, 2) * TJ(3, 3) - TJ(3, 2) * TJ(2, 3))
     DJ2 = TJ(1, 2) * (TJ(2, 3) * TJ(3, 1) - TJ(3, 3) * TJ(2, 1))
     DJ3 = TJ(1, 3) * (TJ(2, 1) * TJ(3, 2) - TJ(3, 1) * TJ(2, 2))
     DJ = DJ1 + DJ2 + DJ3
     '--- Inverse of the Jacobian AJ()
     AJ(1, 1) = (TJ(2, 2) * TJ(3, 3) - TJ(2, 3) * TJ(3, 2)) / DJ
     AJ(1, 2) = (TJ(3, 2) * TJ(1, 3) - TJ(3, 3) * TJ(1, 2)) / DJ
     AJ(1, 3) = (TJ(1, 2) * TJ(2, 3) - TJ(1, 3) * TJ(2, 2)) / DJ
     AJ(2, 1) = (TJ(2, 3) * TJ(3, 1) - TJ(2, 1) * TJ(3, 3)) / DJ
     AJ(2, 2) = (TJ(1, 1) * TJ(3, 3) - TJ(1, 3) * TJ(3, 1)) / DJ
     AJ(2, 3) = (TJ(1, 3) * TJ(2, 1) - TJ(1, 1) * TJ(2, 3)) / DJ
     AJ(3, 1) = (TJ(2, 1) * TJ(3, 2) - TJ(2, 2) * TJ(3, 1)) / DJ
     AJ(3, 2) = (TJ(1, 2) * TJ(3, 1) - TJ(1, 1) * TJ(3, 2)) / DJ
     AJ(3, 3) = (TJ(1, 1) * TJ(2, 2) - TJ(1, 2) * TJ(2, 1)) / DJ
     '--- H() Matrix relates local derivatives of  u  to local
     '    displacements  q
     FOR I = 1 TO 9
        FOR J = 1 TO 24
           H(I, J) = 0
        NEXT J
     NEXT I
     FOR I = 1 TO 3
        FOR J = 1 TO 3
           IR = 3 * (I - 1) + J
           FOR K = 1 TO 8
              IC = 3 * (K - 1) + I
              H(IR, IC) = GN(J, K)
           NEXT K
        NEXT J
     NEXT I
     '--- G() Matrix relates strains to local derivatives of  u
     FOR I = 1 TO 6
        FOR J = 1 TO 9
           G(I, J) = 0
        NEXT J
     NEXT I
     G(1, 1) = AJ(1, 1): G(1, 2) = AJ(1, 2): G(1, 3) = AJ(1, 3)
     G(2, 4) = AJ(2, 1): G(2, 5) = AJ(2, 2): G(2, 6) = AJ(2, 3)
     G(3, 7) = AJ(3, 1): G(3, 8) = AJ(3, 2): G(3, 9) = AJ(3, 3)
     G(4, 4) = AJ(3, 1): G(4, 5) = AJ(3, 2): G(4, 6) = AJ(3, 3)
          G(4, 7) = AJ(2, 1): G(4, 8) = AJ(2, 2): G(4, 9) = AJ(2, 3)
     G(5, 1) = AJ(3, 1): G(5, 2) = AJ(3, 2): G(5, 3) = AJ(3, 3)
          G(5, 7) = AJ(1, 1): G(5, 8) = AJ(1, 2): G(5, 9) = AJ(1, 3)
     G(6, 1) = AJ(2, 1): G(6, 2) = AJ(2, 2): G(6, 3) = AJ(2, 3)
          G(6, 4) = AJ(1, 1): G(6, 5) = AJ(1, 2): G(6, 6) = AJ(1, 3)
     '--- B() Matrix relates strains to  q
     FOR I = 1 TO 6
        FOR J = 1 TO 24
           B(I, J) = 0
           FOR K = 1 TO 9
              B(I, J) = B(I, J) + G(I, K) * H(K, J)
           NEXT K
        NEXT J
     NEXT I
     '--- DB() Matrix relates stresses to  q
     FOR I = 1 TO 6
        FOR J = 1 TO 24
           DB(I, J) = 0
           FOR K = 1 TO 6
              DB(I, J) = DB(I, J) + D(I, K) * B(K, J)
           NEXT K
        NEXT J
     NEXT I
     RETURN
PREFRONT:
        '----- Mark Last Appearance of Node / Make it negative in NOC()
        ' Last appearance is first appearance for reverse element order
        FOR I = 1 TO NN
           FOR J = NE TO 1 STEP -1
              FOR K = 1 TO NEN
                 IF I = NOC(J, K) GOTO LABEL1
              NEXT K
           NEXT J
LABEL1:
           NOC(J, K) = -I
        NEXT I
        '===== Block Size Determination
        NQ = NN * NDN
        DIM IDE(NQ)
        FOR I = 1 TO NQ: IDE(I) = 0: NEXT I
        FOR I = 1 TO NMPC: FOR J = 1 TO 2: IDE(MPC(I, J)) = 1: NEXT J: NEXT I
        IFRON = 0: FOR I = 1 TO NQ: IFRON = IFRON + IDE(I): NEXT I
        IBL = IFRON
        FOR N = 1 TO NE
           INEG = 0
           FOR I = 1 TO NEN
              I1 = NOC(N, I): IA = NDN * (ABS(I1) - 1)
              FOR J = 1 TO NDN
                 IA = IA + 1
                 IF IDE(IA) = 0 THEN
                    IFRON = IFRON + 1: IDE(IA) = 1
                 END IF
              NEXT J
              IF I1 < 0 THEN INEG = INEG + 1
           NEXT I
           IF IBL < IFRON THEN IBL = IFRON
           IFRON = IFRON - NDN * INEG
        NEXT N
        ERASE IDE
        PRINT "Block size = "; IBL
        DIM ISBL(IBL), S(IBL, IBL), IEBL(NEDF), INDX(IBL)
        NFRON = 0: NTOGO = 0: NDCNT = 0
        FOR I = 1 TO IBL: INDX(I) = I: NEXT I
        RETURN
MPCFRON:
        '----- Modifications for Multipoint Constraints by Penalty Method
        FOR I = 1 TO NMPC
           I1 = MPC(I, 1)
           IFL = 0
           FOR J = 1 TO NFRON
              J1 = INDX(J)
              IF I1 = ISBL(J1) THEN
                 IFL = 1: EXIT FOR
              END IF
           NEXT J
           IF IFL = 0 THEN
              NFRON = NFRON + 1: J1 = INDX(NFRON): ISBL(J1) = I1
           END IF
           I2 = MPC(I, 2)
           IFL = 0
           FOR K = 1 TO NFRON
              K1 = INDX(K)
              IF K1 = ISBL(K1) THEN
                 IFL = 1: EXIT FOR
              END IF
           NEXT K
           IF IFL = 0 THEN
              NFRON = NFRON + 1: K1 = INDX(NFRON): ISBL(K1) = I2
           END IF
           '----- Stiffness Modification
           S(J1, J1) = S(J1, J1) + CNST * BT(I, 1) ^ 2
           S(K1, K1) = S(K1, K1) + CNST * BT(I, 2) ^ 2
           S(J1, K1) = S(J1, K1) + CNST * BT(I, 1) * BT(I, 2)
           S(K1, J1) = S(J1, K1)
           '----- Force Modification
           F(I1) = F(I1) + CNST * BT(I, 3) * BT(I, 1)
           F(I2) = F(I2) + CNST * BT(I, 3) * BT(I, 2)
        NEXT I
        RETURN
FRONT:
        '----- Frontal Method Assembly and Elimination
        FIELD 3, 4 AS VAR$, 4 AS COEF$     'Field definition
'----------------  Assembly of Element N  --------------------
        FOR I = 1 TO NEN
           I1 = NOC(N, I): IA = ABS(I1): IS1 = SGN(I1)
           IDF = NDN * (IA - 1): IE1 = NDN * (I - 1)
           FOR J = 1 TO NDN
              IDF = IDF + 1: IE1 = IE1 + 1: IFL = 0
              IF NFRON > NTOGO THEN
                 FOR II = NTOGO + 1 TO NFRON
                    IX = INDX(II)
                    IF IDF = ISBL(IX) THEN
                       IFL = 1: EXIT FOR
                    END IF
                 NEXT II
              END IF
              IF IFL = 0 THEN
                 NFRON = NFRON + 1: II = NFRON: IX = INDX(II)
              END IF
              ISBL(IX) = IDF: IEBL(IE1) = IX
              IF IS1 = -1 THEN
                 NTOGO = NTOGO + 1
                 ITEMP = INDX(NTOGO)
                 INDX(NTOGO) = INDX(II)
                 INDX(II) = ITEMP
              END IF
           NEXT J
        NEXT I
        FOR I = 1 TO NEDF
           I1 = IEBL(I)
           FOR J = 1 TO NEDF
              J1 = IEBL(J)
              S(I1, J1) = S(I1, J1) + SE(I, J)
           NEXT J
        NEXT I
'------------------------------------------------------------------
     IF NDCNT < ND THEN
'-----  Modification for displacement BCs / Penalty Approach  -----
        FOR I = 1 TO NTOGO
           I1 = INDX(I)
           IG = ISBL(I1)
              FOR J = 1 TO ND
                 IF IG = NU(J) THEN
                    S(I1, I1) = S(I1, I1) + CNST
                    F(IG) = F(IG) + CNST * U(J)
                    NDCNT = NDCNT + 1       'Counter for check
                    EXIT FOR
                 END IF
              NEXT J
        NEXT I
     END IF
'------------   Elimination of completed variables   ---------------
        NTG1 = NTOGO
        FOR II = 1 TO NTG1
           IPV = INDX(1): IPG = ISBL(IPV)
           PIVOT = S(IPV, IPV)
        '-----  Write separator "0" and PIVOT value to disk  -----
           LSET COEF$ = MKS$(PIVOT): LSET VAR$ = MKL$(0)
           ICOUNT = ICOUNT + 1
           PUT 3, ICOUNT: S(IPV, IPV) = 0
           FOR I = 2 TO NFRON
              I1 = INDX(I): IG = ISBL(I1)
              IF S(I1, IPV) <> 0 THEN
                  C = S(I1, IPV) / PIVOT: S(I1, IPV) = 0
                  FOR J = 2 TO NFRON
                     J1 = INDX(J)
                     IF S(IPV, J1) <> 0 THEN
                        S(I1, J1) = S(I1, J1) - C * S(IPV, J1)
                     END IF
                  NEXT J
                  F(IG) = F(IG) - C * F(IPG)
              END IF
           NEXT I
           FOR J = 2 TO NFRON
        '-----  Write Variable# and Reduced Coeff/PIVOT to disk  -----
              J1 = INDX(J)
              IF S(IPV, J1) <> 0 THEN
                 ICOUNT = ICOUNT + 1: IBA = ISBL(J1)
                 LSET COEF$ = MKS$(S(IPV, J1) / PIVOT)
                 LSET VAR$ = MKL$(IBA)
                 PUT 3, ICOUNT: S(IPV, J1) = 0
              END IF
           NEXT J
           ICOUNT = ICOUNT + 1
        '-----  Write Eliminated Variable# and RHS/PIVOT to disk  -----
           LSET COEF$ = MKS$(F(IPG) / PIVOT): LSET VAR$ = MKL$(IPG)
           F(IPG) = 0
           PUT 3, ICOUNT
        '----- (NTOGO) into (1); (NFRON) into (NTOGO)
        '----- IPV into (NFRON) and reduce front & NTOGO sizes by 1
           IF NTOGO > 1 THEN
              INDX(1) = INDX(NTOGO)
           END IF
           INDX(NTOGO) = INDX(NFRON): INDX(NFRON) = IPV
           NFRON = NFRON - 1: NTOGO = NTOGO - 1
        NEXT II
        RETURN
BACKSUB:
        '===== Backsubstitution
        FIELD 3, 4 AS VAR$, 4 AS COEF$     'Field definition
STEP1:
        IF ICOUNT <= 0 THEN RETURN
        GET 3, ICOUNT: ICOUNT = ICOUNT - 1
        N1 = CVL(VAR$): F(N1) = CVS(COEF$)
STEP2:
        GET 3, ICOUNT: ICOUNT = ICOUNT - 1
        N2 = CVL(VAR$)
        IF N2 = 0 GOTO STEP1
        F(N1) = F(N1) - CVS(COEF$) * F(N2)
        GOTO STEP2
        RETURN

