     CLS
     PRINT "*************   PLOT2D  ***************"
     PRINT "*       TWO DIMENSIONAL PLOT          *"
     PRINT "* T.R.Chandrupatla and A.D.Belegundu  *"
     PRINT "***************************************"
     '========  Screen for Graphics  =========
     DEFINT I-N: F$ = "####.##"
     LOCATE 7, 1
     INPUT "Plot File Name <path fn.ext>  ", FILE$
     OPEN FILE$ FOR INPUT AS #1
     LINE INPUT #1, D$: INPUT #1, TITLE$
     LINE INPUT #1, D$: INPUT #1, NN, NE, NM, NDIM, NEN, NDN
     LINE INPUT #1, D$: INPUT #1, ND, NL, NMPC
     IF NDIM <> 2 THEN
        PRINT "THE PROGRAM SUPPORTS TWO DIMENSIONAL PLOTS ONLY"
        PRINT "THE DIMENSION OF THE DATA IS  "; NDIM
        END
     END IF
     DIM X(NN, NDIM), NOC(NE, NEN)
     '=============  READ DATA  ===============
     LINE INPUT #1, D$
     FOR I = 1 TO NN: INPUT #1, N: FOR J = 1 TO NDIM
     INPUT #1, X(N, J): NEXT J: NEXT I
     LINE INPUT #1, D$
     FOR I = 1 TO NE: INPUT #1, N: FOR J = 1 TO NEN
     INPUT #1, NOC(N, J): NEXT J: INPUT #1, D$
     NEXT I
     CLOSE #1
     XMAX = X(1, 1): YMAX = X(1, 2): XMIN = X(1, 1): YMIN = X(1, 2)
     FOR I = 2 TO NN
        IF XMAX < X(I, 1) THEN XMAX = X(I, 1)
        IF YMAX < X(I, 2) THEN YMAX = X(I, 2)
        IF XMIN > X(I, 1) THEN XMIN = X(I, 1)
        IF YMIN > X(I, 2) THEN YMIN = X(I, 2)
     NEXT I
     SCREEN 9
     CLS : XL = (XMAX - XMIN): YL = (YMAX - YMIN)
     X0 = XMIN - XL / 10: Y0 = YMIN - YL / 10
     VIEW (101, 1)-(638, 300), , 2
     ASP = 1.3  'Ratio of larger side to smaller side
                'of border rectangle  <change this ratio if needed>
     IF ASP * YL > XL THEN XL = ASP * YL ELSE YL = XL / ASP
     XMAX = X0 + 1.2 * XL: YMAX = Y0 + 1.2 * YL
     WINDOW (X0, Y0)-(XMAX, YMAX)
     LOCATE 1, 5: PRINT USING F$; YMAX
     LOCATE 23, 73: PRINT USING F$; XMAX
     LOCATE 23, 10: PRINT USING F$; X0
     LOCATE 21, 5: PRINT USING F$; Y0
     IFLIP = 1: GOSUB DRAWEL
     DO
       A$ = INKEY$
       IF A$ = CHR$(0) + CHR$(59) THEN
          IFLIP = -IFLIP
          IF IFLIP = 1 THEN GOSUB DRAWEL
          IF IFLIP = -1 THEN GOSUB NODES
       END IF
       IF A$ = CHR$(0) + CHR$(60) THEN EXIT DO
     LOOP
     END
DRAWEL:
     '===========  Draw Elements  ================
     CLS : FOR IE = 1 TO NE
     FOR II = 1 TO NEN
     X1 = X(NOC(IE, II), 1): Y1 = X(NOC(IE, II), 2)
     IF II = NEN THEN
        X2 = X(NOC(IE, 1), 1): Y2 = X(NOC(IE, 1), 2)
        LINE (X1, Y1)-(X2, Y2)
     ELSE
        X2 = X(NOC(IE, II + 1), 1): Y2 = X(NOC(IE, II + 1), 2)
        LINE (X1, Y1)-(X2, Y2)
     END IF
     NEXT II: NEXT IE
     LOCATE 24, 8: PRINT "F1 - to Show or Remove Node Numbers";
     LOCATE 24, 50: PRINT "F2 - to Exit Program";
     RETURN
NODES:
     '============  Show Node Numbers  ==============
     FOR I = 1 TO NN
        ICOL = (100 + (537 * (X(I, 1) - X0) / (1.2 * XL))) / 8
        IROW = (300 - 299 * (X(I, 2) - Y0) / (1.2 * YL)) / 14 + .5
        LOCATE IROW, ICOL: PRINT I;
     NEXT I
     RETURN

