VERSION 5.00
Begin VB.Form frmContourB 
   Caption         =   "ContourB"
   ClientHeight    =   4920
   ClientLeft      =   60
   ClientTop       =   348
   ClientWidth     =   6912
   LinkTopic       =   "Form1"
   ScaleHeight     =   4920
   ScaleWidth      =   6912
   StartUpPosition =   3  'Windows Default
   WindowState     =   2  'Maximized
   Begin MSComDlg.CommonDialog cdlDialog 
      Left            =   8040
      Top             =   5640
      _ExtentX        =   847
      _ExtentY        =   847
      _Version        =   393216
      FilterIndex     =   2
   End
   Begin VB.PictureBox picBox2 
      Height          =   3495
      Left            =   6480
      ScaleHeight     =   3444
      ScaleWidth      =   1764
      TabIndex        =   9
      Top             =   120
      Width           =   1815
   End
   Begin VB.CommandButton cmdLLeft 
      Caption         =   "ZoomLL"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6480
      TabIndex        =   8
      Top             =   4320
      Width           =   855
   End
   Begin VB.CommandButton cmdLRight 
      Caption         =   "ZoomLR"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7440
      TabIndex        =   7
      Top             =   4320
      Width           =   855
   End
   Begin VB.CommandButton cmdOriginal 
      Caption         =   "Original"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7440
      TabIndex        =   6
      Top             =   5160
      Width           =   855
   End
   Begin VB.CommandButton cmdPrevious 
      Caption         =   "Pevious"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6480
      TabIndex        =   5
      Top             =   5160
      Width           =   855
   End
   Begin VB.CommandButton cmdURight 
      Caption         =   "ZoomUR"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7440
      TabIndex        =   4
      Top             =   3720
      Width           =   855
   End
   Begin VB.CommandButton cmdULeft 
      Caption         =   "ZoomUL"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6480
      TabIndex        =   3
      Top             =   3720
      Width           =   855
   End
   Begin VB.CommandButton cmdEnd 
      Caption         =   "END"
      Height          =   375
      Left            =   7440
      TabIndex        =   2
      Top             =   5880
      Width           =   855
   End
   Begin VB.CommandButton cmdPlot 
      Caption         =   "START"
      Height          =   375
      Left            =   6480
      TabIndex        =   1
      Top             =   5880
      Width           =   855
   End
   Begin VB.PictureBox picBox 
      Height          =   6135
      Left            =   120
      ScaleHeight     =   6084
      ScaleWidth      =   6204
      TabIndex        =   0
      Top             =   120
      Width           =   6255
   End
End
Attribute VB_Name = "frmContourB"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
'********         CONTOURB          **********
'*     CONTOUR PLOTTING - CONTOUR BANDS      *
'*     FOR 2D TRIANGLES AND QUADRILATERALS   *
'*     T.R.Chandrupatla and A.D.Belegundu    *
'*********************************************
DefInt I-N
Dim XMIN, XMAX, YMIN, YMAX, XL, XH, YL, YH, XOL, XOH, YOL, YOH
Dim NN, NE, NM, NDIM, NEN, NDN, ND, NL, NMPC
Dim X(), NOC(), FF(), NCON(), XX(), YY(), U(), IC()
Dim FMIN, FMAX, STP, MMIN, NMIN, MMAX, NMAX, DU
Dim X1, X2, X3, Y1, Y2, Y3, U1, U2, U3
Dim File1 As String, File2 As String, Dummy As String, Title As String
Const AL = 0.67, NCL = 10
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdPlot_Click()
     Call InputData
     Call DrawLimits(XMIN, YMIN, XMAX, YMAX)
     Call DrawLegend
     Call DrawContours
     cmdPlot.Enabled = False
     cmdULeft.Enabled = True
     cmdURight.Enabled = True
     cmdLLeft.Enabled = True
     cmdLRight.Enabled = True
End Sub
Private Sub InputData()
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     cdlDialog.Filter = "All Files (*.*)|*.*|<Nodal val file> FileName.n|*.n"
     cdlDialog.FileName = ""
     cdlDialog.ShowOpen
     File2 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     If NDIM <> 2 Or NEN < 3 Or NEN > 4 Then
        picBox.Print "This program supports triangular and quadrilateral"
        picBox.Print "Elements only."
        End
     End If
     ReDim X(NN, NDIM), NOC(NE, NEN), FF(NN), NCON(NE, NEN)
     ReDim XX(3), YY(3), U(3), IC(10)
     '=============  COLOR DATA  ===============
     IC(1) = 13: IC(2) = 5: IC(3) = 9: IC(4) = 1: IC(5) = 2
     IC(6) = 10: IC(7) = 14: IC(8) = 6: IC(9) = 4: IC(10) = 12
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN: Input #1, N: For J = 1 To NDIM
     Input #1, X(N, J): Next J: Next I
     '----- Connectivity
     Line Input #1, Dummy
     For I = 1 To NE: Input #1, N: For J = 1 To NEN
     Input #1, NOC(N, J): Next J: Line Input #1, Dummy
     Next I
     Close #1
     Open File2 For Input As #2
     '----- Nodal Values
     Line Input #2, Dummy
     For I = 1 To NN
     Input #2, FF(I): Next I
     Close #2
End Sub
Private Sub DrawLimits(XMIN, YMIN, XMAX, YMAX)
     XMAX = X(1, 1): YMAX = X(1, 2): FMAX = FF(1)
     XMIN = X(1, 1): YMIN = X(1, 2): FMIN = FF(1)
     For I = 2 To NN
        If XMAX < X(I, 1) Then XMAX = X(I, 1)
        If YMAX < X(I, 2) Then YMAX = X(I, 2)
        If XMIN > X(I, 1) Then XMIN = X(I, 1)
        If YMIN > X(I, 2) Then YMIN = X(I, 2)
        If FMAX < FF(I) Then FMAX = FF(I)
        If FMIN > FF(I) Then FMIN = FF(I)
     Next I
     STP = (FMAX - FMIN) / NCL
     XL = (XMAX - XMIN): YL = (YMAX - YMIN)
     A = XL: If A < YL Then A = YL
     XB = 0.5 * (XMIN + XMAX)
     YB = 0.5 * (YMIN + YMAX)
     XMIN = XB - 0.55 * A: XMAX = XB + 0.55 * A
     YMIN = YB - 0.55 * A: YMAX = YB + 0.55 * A
     XL = XMIN: YL = YMIN: XH = XMAX: YH = YMAX
     XOL = XL: YOL = YL: XOH = XH: YOH = YH
End Sub
Private Sub DrawContours()
     picBox.ScaleMode = 3
     picBox.Cls
     MMIN = picBox.ScaleLeft
     MMAX = MMIN + picBox.ScaleWidth
     NMIN = picBox.ScaleTop
     NMAX = NMIN + picBox.ScaleHeight
     '===========  Contour Plotting  ===========
     DU = STP
     For IE = 1 To NE
        If NEN = 3 Then
           For IEN = 1 To NEN
              IEE = NOC(IE, IEN)
              U(IEN) = FF(IEE)
              XX(IEN) = X(IEE, 1)
              YY(IEN) = X(IEE, 2)
           Next IEN
           Call ElementPlot
        ElseIf NEN = 4 Then
           X5 = 0: Y5 = 0: U5 = 0
           For IT = 1 To NEN
              NIT = NOC(IE, IT)
              X5 = X5 + 0.25 * X(NIT, 1)
              Y5 = Y5 + 0.25 * X(NIT, 2)
              U5 = U5 + 0.25 * FF(NIT)
           Next IT
           For IT = 1 To NEN
              IT1 = IT + 1: If IT1 > 4 Then IT1 = 1
              XX(1) = X5: YY(1) = Y5: U(1) = U5
              NIE = NOC(IE, IT)
              XX(2) = X(NIE, 1): YY(2) = X(NIE, 2): U(2) = FF(NIE)
              NIE = NOC(IE, IT1)
              XX(3) = X(NIE, 1): YY(3) = X(NIE, 2): U(3) = FF(NIE)
              Call ElementPlot
           Next IT
        Else
           picBox.Print "NUMBER OF ELEMENT NODES > 4 IS NOT SUPPORTED"
           End
        End If
     Next IE
End Sub
Private Sub DrawLegend()
    picBox2.Scale (0, 1)-(1, 0)
    Y = 0
    F = FMIN + 0.5 * STP
    For I = 1 To NCL
       picBox2.Line (0, Y)-(0.1, Y + 0.1), QBColor(IC(I)), BF
       picBox2.CurrentX = 0.15
       picBox2.CurrentY = Y + 0.08
       picBox2.Print Format(F, "0.000E+00  ")
       Y = Y + 0.1
       F = F + STP
    Next I
End Sub
Private Sub ElementPlot()
        For I = 1 To 2
           C = YY(I): II = I
           For J = I + 1 To 3
              If C > YY(J) Then
                 C = YY(J): II = J
              End If
           Next J
           YY(II) = YY(I): YY(I) = C
           C1 = XX(II): XX(II) = XX(I): XX(I) = C1
           C1 = U(II): U(II) = U(I): U(I) = C1
        Next I
        BET = (YY(2) - YY(1)) / (YY(3) - YY(1))
        XT = XX(1) + BET * (XX(3) - XX(1))
        YT = YY(2)
        UT = U(1) + BET * (U(3) - U(1))
        '----- Interchange to make X2 < X3
        If XT > XX(2) Then
           X2 = XX(2): Y2 = YY(2): U2 = U(2)
           X3 = XT: Y3 = YT: U3 = UT
        Else
           X2 = XT: Y2 = YT: U2 = UT
           X3 = XX(2): Y3 = YY(2): U3 = U(2)
        End If
        A = XH - XL: B = YH - YL
        ML = MMAX - MMIN
        NL = NMAX - NMIN
        X2 = MMIN + (X2 - XL) / A * ML
        Y2 = NMIN + (YH - Y2) / B * NL
        X3 = MMIN + (X3 - XL) / A * ML
        Y3 = NMIN + (YH - Y3) / B * NL
        '----- Lower Triangle -----
        X1 = XX(1): Y1 = YY(1): U1 = U(1): ISTEP = -1
        X1 = MMIN + (X1 - XL) / A * ML
        Y1 = NMIN + (YH - Y1) / B * NL
        IY1 = Int(Y1)
        IY2 = NMAX - Int(NMAX - Y2)
        If IY2 <= IY1 Then Call TriPixel(IY1, IY2, ISTEP)
        '----- Upper Triangle -----
        X1 = XX(3): Y1 = YY(3): U1 = U(3): ISTEP = 1
        X1 = MMIN + (X1 - XL) / A * ML
        Y1 = NMIN + (YH - Y1) / B * NL
        IY1 = NMAX - Int(NMAX - Y1)
        IY2 = Int(Y2)
        If IY1 <= IY2 Then Call TriPixel(IY1, IY2, ISTEP)
End Sub
Private Sub TriPixel(IY1, IY2, ISTEP)
        For I = IY1 To IY2 Step ISTEP
            A = X1 + (I - Y1) / (Y2 - Y1) * (X2 - X1)
            B = X1 + (I - Y1) / (Y3 - Y1) * (X3 - X1)
            UA = U1 + (I - Y1) / (Y2 - Y1) * (U2 - U1)
            UB = U1 + (I - Y1) / (Y3 - Y1) * (U3 - U1)
            IA = MMAX - Int(MMAX - A)
            IB = Int(B)
            If IA <= IB Then Call CoLine(UA, UB, A, B, IA, IB, I)
        Next I
End Sub
Private Sub CoLine(UA, UB, A, B, IA, IB, I)
        For J = IA To IB
            UU = UA + (J - IA) / (B - A) * (UB - UA)
            IU = Int((UU - FMIN) / DU) + 1
            If IU < 0 Then IU = 1
            If IU > 10 Then IU = 10
            ICO = IC(IU)
            picBox.PSet (J, I), QBColor(ICO)
        Next J
End Sub
Private Sub cmdPrevious_Click()
     T = XL: XL = XOL: XOL = T
     T = YL: YL = YOL: YOL = T
     T = XH: XH = XOH: XOH = T
     T = YH: YH = YOH: YOH = T
     Call DrawContours
     cmdOriginal.Enabled = True
End Sub
Private Sub cmdOriginal_Click()
     XL = XMIN: YL = YMIN: XH = XMAX: YH = YMAX
     XOL = XL: YOL = YL: XOH = XH: YOH = YH
     Call DrawContours
     cmdOriginal.Enabled = False
     cmdPrevious.Enabled = False
End Sub
Private Sub cmdULeft_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XH = (1 - AL) * XL + AL * XH
     YL = AL * YL + (1 - AL) * YH
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdURight_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XL = AL * XL + (1 - AL) * XH
     YL = AL * YL + (1 - AL) * YH
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdLLeft_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XH = (1 - AL) * XL + AL * XH
     YH = (1 - AL) * YL + AL * YH
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdLRight_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XL = AL * XL + (1 - AL) * XH
     YH = (1 - AL) * YL + AL * YH
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
