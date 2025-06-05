     '**********   PROGRAM CONTOURA   *********
     '*             CONTOUR LINES             *
     '*   T.R.Chandrupatla and A.D.Belegundu  *
     '*****************************************
     DEFINT I-N: CLS : COLOR 1, 3
     LOCATE 1, 1: PRINT "CONTOUR LINE PLOTTING"; SPACE$(30);
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
     DIM NCON(NE, NEN), XX(3), YY(3), U(3), IC(10)
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
     VIEW (95, 10)-(600, 300), , 12
     AA = 505 * ASP / 290
     IF XL / YL > AA THEN YL = XL / AA
     IF XL / YL < AA THEN XL = YL * AA
     XMAX = X0 + 1.2 * XL: YMAX = Y0 + 1.2 * YL
     WINDOW (X0, Y0)-(XMAX, YMAX)
     LOCATE 1, 5: PRINT USING F$; YMAX
     LOCATE 23, 73: PRINT USING F$; XMAX
     LOCATE 23, 10: PRINT USING F$; X0
     LOCATE 21, 5: PRINT USING F$; Y0
     '----- Print Color Legend -----
     Y1 = .38 * YMAX + .62 * Y0
     X1 = .03 * XMAX + .97 * X0
     DY = .077 * (YMAX - Y1)
     IR = 15: FT = FMIN
     LOCATE IR, 3: PRINT USING "##.#^^^^"; FMIN
     FOR I = 1 TO NCL
        ICO = IC(I): IR = IR - 1
        LINE (X0, Y1)-(X1, Y1), ICO, BF
        Y1 = Y1 + DY: FT = FT + STP
        LOCATE IR, 3: PRINT USING "##.#^^^^"; FT
     NEXT I
'=============  Find Boundary Lines  ===============
     'Edges defined by nodes in NOC to nodes in NCON
     FOR IE = 1 TO NE
     FOR I = 1 TO NEN
     I1 = I + 1: IF I1 > NEN THEN I1 = 1
     NCON(IE, I) = NOC(IE, I1): NEXT I: NEXT IE
     FOR IE = 1 TO NE
     FOR I = 1 TO NEN
     I1 = NCON(IE, I): I2 = NOC(IE, I)
     INDX = 0
     FOR JE = IE + 1 TO NE
     FOR J = 1 TO NEN
     IF NCON(JE, J) = 0 GOTO STEP1
     IF I1 <> NCON(JE, J) AND I1 <> NOC(JE, J) GOTO STEP1
     IF I2 <> NCON(JE, J) AND I2 <> NOC(JE, J) GOTO STEP1
     NCON(JE, J) = 0: INDX = INDX + 1
STEP1:
     NEXT J: NEXT JE
     IF INDX > 0 THEN NCON(IE, I) = 0
     NEXT I: NEXT IE
     '============  Draw Boundary  ==============
     FOR IE = 1 TO NE
     FOR I = 1 TO NEN
     IF NCON(IE, I) > 0 THEN
        I1 = NCON(IE, I): I2 = NOC(IE, I)
        LINE (X(I1, 1), X(I1, 2))-(X(I2, 1), X(I2, 2))
     END IF
     NEXT I: NEXT IE
     '===========  Contour Plotting  ===========
     FOR IE = 1 TO NE
        IF NEN = 3 THEN
           FOR IEN = 1 TO NEN
              IEE = NOC(IE, IEN)
              U(IEN) = FF(IEE)
              XX(IEN) = X(IEE, 1)
              YY(IEN) = X(IEE, 2)
           NEXT IEN
           GOSUB LPLOT
        ELSEIF NEN = 4 THEN
           XB = 0: YB = 0: UB = 0
           FOR IT = 1 TO NEN
              NIT = NOC(IE, IT)
              XB = XB + .25 * X(NIT, 1)
              YB = YB + .25 * X(NIT, 2)
              UB = UB + .25 * FF(NIT)
           NEXT IT
           FOR IT = 1 TO NEN
              IT1 = IT + 1: IF IT1 > 4 THEN IT1 = 1
              XX(1) = XB: YY(1) = YB: U(1) = UB
              NIE = NOC(IE, IT)
              XX(2) = X(NIE, 1): YY(2) = X(NIE, 2): U(2) = FF(NIE)
              NIE = NOC(IE, IT1)
              XX(3) = X(NIE, 1): YY(3) = X(NIE, 2): U(3) = FF(NIE)
              GOSUB LPLOT
           NEXT IT
        ELSE
           PRINT "NUMBER OF ELEMENT NODES > 4 IS NOT SUPPORTED"
           END
        END IF
      NEXT IE
     END
LPLOT:
   'THREE POINTS IN ASCENDING ORDER
        FOR I = 1 TO 2
           C = U(I): II = I
           FOR J = I + 1 TO 3
              IF C > U(J) THEN
                 C = U(J): II = J
              END IF
           NEXT J
           U(II) = U(I): U(I) = C
           C1 = XX(II): XX(II) = XX(I): XX(I) = C1
           C1 = YY(II): YY(II) = YY(I): YY(I) = C1
        NEXT I
           II = INT((U(1) - FMIN) / STP): UT = FMIN + II * STP
        DO
           II = II + 1: IF II > 10 THEN II = 10
           ICO = IC(II)
           UT = UT + STP
        IF UT >= U(3) THEN EXIT DO
           X1 = ((U(3) - UT) * XX(1) + (UT - U(1)) * XX(3)) / (U(3) - U(1))
           Y1 = ((U(3) - UT) * YY(1) + (UT - U(1)) * YY(3)) / (U(3) - U(1))
           L = 1: IF UT > U(2) THEN L = 3
           X2 = ((U(L) - UT) * XX(2) + (UT - U(2)) * XX(L)) / (U(L) - U(2))
           Y2 = ((U(L) - UT) * YY(2) + (UT - U(2)) * YY(L)) / (U(L) - U(2))
           LINE (X1, Y1)-(X2, Y2), ICO
        LOOP
   RETURN

