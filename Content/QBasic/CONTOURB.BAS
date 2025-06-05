     '**********   PROGRAM CONTOURB   *********
     '*             CONTOUR BANDS             *
     '*   T.R.Chandrupatla and A.D.Belegundu  *
     '*****************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "CONTOUR BAND PLOTTING"; SPACE$(30);
     PRINT "(C)  Chandrupatla & Belegundu": : COLOR 7, 0
     VIEW PRINT 2 TO 25: PRINT
     INPUT "Finite Element Input File ", FILE1$
     OPEN FILE1$ FOR INPUT AS #1
     INPUT "File Name for Contour Input Data ", FILE2$
     LINE INPUT #1, D$: INPUT #1, TITLE$: LINE INPUT #1, D$
     INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$
     INPUT #1, ND, NL, NMPC
     IF NDIM <> 2 OR NEN < 3 OR NEN > 4 THEN
        PRINT "This program supports triangular and quadrilateral"
        PRINT "Elements only."
        END
     END IF
     DIM X(NN, NDIM), NOC(NE, NEN), FF(NN)
     DIM XX(3), YY(3), U(3), IC(10)
     '=============  COLOR DATA  ===============
     FOR I = 1 TO 10: READ IC(I): NEXT I
     DATA 13,5,9,1,2,10,14,6,4,12
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
     FOR I = 1 TO NE: INPUT #1, N: FOR J = 1 TO NEN
     INPUT #1, NOC(N, J): NEXT J: LINE INPUT #1, D$
     NEXT I
     CLOSE #1
     OPEN FILE2$ FOR INPUT AS #2
     '----- Nodal Values
     LINE INPUT #2, D$
     FOR I = 1 TO NN
     INPUT #2, FF(I): NEXT I
     CLOSE #2
     XMAX = X(1, 1): YMAX = X(1, 2): FMAX = FF(1)
     XMIN = X(1, 1): YMIN = X(1, 2): FMIN = FF(1)
     FOR I = 2 TO NN
     IF XMAX < X(I, 1) THEN XMAX = X(I, 1)
     IF YMAX < X(I, 2) THEN YMAX = X(I, 2)
     IF XMIN > X(I, 1) THEN XMIN = X(I, 1)
     IF YMIN > X(I, 2) THEN YMIN = X(I, 2)
     IF FMAX < FF(I) THEN FMAX = FF(I)
     IF FMIN > FF(I) THEN FMIN = FF(I)
     NEXT I
     NCL = 10
     STP = (FMAX - FMIN) / NCL
     XL = (XMAX - XMIN): YL = (YMAX - YMIN)
     X0 = XMIN - XL / 10: Y0 = YMIN - YL / 10
     SCREEN 9: F$ = "####.##": CLS
     ASP = .65 '***  Change Aspect Ratio if figure is not proportional
     MMIN = 95: MMAX = 600: NMIN = 10: NMAX = 300
     VIEW (MMIN, NMIN)-(MMAX, NMAX), , 12
     ML = MMAX - MMIN: NL = NMAX - NMIN
     AA = ML * ASP / NL
     IF XL / YL > AA THEN YL = XL / AA
     IF XL / YL < AA THEN XL = YL * AA
     XL = 1.2 * XL: YL = 1.2 * YL
     XMAX = X0 + XL: YMAX = Y0 + YL
     WINDOW (X0, Y0)-(XMAX, YMAX)
     LOCATE 1, 5: PRINT USING F$; YMAX
     LOCATE 23, 73: PRINT USING F$; XMAX
     LOCATE 23, 10: PRINT USING F$; X0
     LOCATE 22, 5: PRINT USING F$; Y0
'-----  Color Legend  -----
     Y1 = .33 * YMAX + .67 * Y0
     X1 = .02 * XMAX + .98 * X0
     DY = .075 * (YMAX - Y1)
     IR = 15: FT = FMIN
     LOCATE IR, 3: PRINT USING "##.#^^^^"; FMIN
     FOR I = 1 TO NCL
        ICO = IC(I): IR = IR - 1
        LINE (X0, Y1)-(X1, Y1 + DY), ICO, BF
        Y1 = Y1 + DY: FT = FT + STP
        LOCATE IR, 3: PRINT USING "##.#^^^^"; FT
     NEXT I
     WINDOW
     VIEW SCREEN (MMIN, NMIN)-(MMAX, NMAX), , 12
'===========  Contour Plotting  ===========
     DU = STP
     FOR IE = 1 TO NE
        IF NEN = 3 THEN
           FOR IEN = 1 TO NEN
              IEE = NOC(IE, IEN)
              U(IEN) = FF(IEE)
              XX(IEN) = X(IEE, 1)
              YY(IEN) = X(IEE, 2)
           NEXT IEN
           GOSUB TRIPLOT
        ELSEIF NEN = 4 THEN
           X5 = 0: Y5 = 0: U5 = 0
           FOR IT = 1 TO NEN
              NIT = NOC(IE, IT)
              X5 = X5 + .25 * X(NIT, 1)
              Y5 = Y5 + .25 * X(NIT, 2)
              U5 = U5 + .25 * FF(NIT)
           NEXT IT
           FOR IT = 1 TO NEN
              IT1 = IT + 1: IF IT1 > 4 THEN IT1 = 1
              XX(1) = X5: YY(1) = Y5: U(1) = U5
              NIE = NOC(IE, IT)
              XX(2) = X(NIE, 1): YY(2) = X(NIE, 2): U(2) = FF(NIE)
              NIE = NOC(IE, IT1)
              XX(3) = X(NIE, 1): YY(3) = X(NIE, 2): U(3) = FF(NIE)
              GOSUB TRIPLOT
           NEXT IT
        ELSE
           PRINT "NUMBER OF ELEMENT NODES > 4 IS NOT SUPPORTED"
           END
        END IF
     NEXT IE
     END
TRIPLOT:
        FOR I = 1 TO 2
           C = YY(I): II = I
           FOR J = I + 1 TO 3
              IF C > YY(J) THEN
                 C = YY(J): II = J
              END IF
           NEXT J
           YY(II) = YY(I): YY(I) = C
           C1 = XX(II): XX(II) = XX(I): XX(I) = C1
           C1 = U(II): U(II) = U(I): U(I) = C1
        NEXT I
        AL = (YY(2) - YY(1)) / (YY(3) - YY(1))
        XT = XX(1) + AL * (XX(3) - XX(1))
        YT = YY(2)
        UT = U(1) + AL * (U(3) - U(1))
        '----- Interchange to make X2 < X3
        IF XT > XX(2) THEN
           X2 = XX(2): Y2 = YY(2): U2 = U(2)
           X3 = XT: Y3 = YT: U3 = UT
        ELSE
           X2 = XT: Y2 = YT: U2 = UT
           X3 = XX(2): Y3 = YY(2): U3 = U(2)
        END IF
        I23 = INT((X3 - X2) / XL * ML + .5)
        UINC = (U3 - U2) / I23
        '----- Lower Triangle -----
        X1 = XX(1): Y1 = YY(1): U1 = U(1): IYS = 1
        IY2 = NMIN + INT((YMAX - Y2) / YL * NL + .5)
        IY1 = NMIN + INT((YMAX - Y1) / YL * NL + .5)
        IF IY2 < IY1 THEN GOSUB TRIPIXEL
          '----- Upper Triangle -----
        X1 = XX(3): Y1 = YY(3): U1 = U(3): IYS = -1
        IY1 = NMIN + INT((YMAX - Y1) / YL * NL + .5)
        IF IY1 < IY2 THEN GOSUB TRIPIXEL
        RETURN
TRIPIXEL:
        FOR I = IY2 TO IY1 STEP IYS
            YI = Y0 + (NMAX - I) / NL * YL
            AL = (YI - Y1) / (Y2 - Y1)
            A = X1 + AL * (X2 - X1)
            B = X1 + AL * (X3 - X1)
            UA = U1 + AL * (U2 - U1)
            UB = U1 + AL * (U3 - U1)
            IA = MMIN + INT((A - X0) / XL * ML + .5)
            IB = MMIN + INT((B - X0) / XL * ML + .5)
            GOSUB COLINE
        NEXT I
        RETURN
COLINE:
        UU = UA - UINC
        FOR J = IA TO IB
            UU = UU + UINC
            IU = INT((UU - FMIN) / DU) + 1
            IF IU > 10 THEN IU = 10
            ICO = IC(IU)
            PSET (J, I), ICO
        NEXT J
        RETURN

