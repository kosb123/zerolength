VERSION 5.00
Begin VB.Form frmContourA 
   Caption         =   "ContourA"
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
   Begin VB.CommandButton cmdLLeft 
      Caption         =   "ZoomLL"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6600
      TabIndex        =   8
      Top             =   2400
      Width           =   855
   End
   Begin VB.CommandButton cmdLRight 
      Caption         =   "ZoomLR"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7560
      TabIndex        =   7
      Top             =   2400
      Width           =   855
   End
   Begin VB.CommandButton cmdOriginal 
      Caption         =   "Original"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7560
      TabIndex        =   6
      Top             =   3360
      Width           =   855
   End
   Begin VB.CommandButton cmdPrevious 
      Caption         =   "Pevious"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6600
      TabIndex        =   5
      Top             =   3360
      Width           =   855
   End
   Begin VB.CommandButton cmdURight 
      Caption         =   "ZoomUR"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7560
      TabIndex        =   4
      Top             =   1800
      Width           =   855
   End
   Begin VB.CommandButton cmdULeft 
      Caption         =   "ZoomUL"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6600
      TabIndex        =   3
      Top             =   1800
      Width           =   855
   End
   Begin VB.CommandButton cmdEnd 
      Caption         =   "END"
      Height          =   375
      Left            =   7560
      TabIndex        =   2
      Top             =   4560
      Width           =   855
   End
   Begin VB.CommandButton cmdPlot 
      Caption         =   "START"
      Height          =   375
      Left            =   6600
      TabIndex        =   1
      Top             =   4560
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
Attribute VB_Name = "frmContourA"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
'********         CONTOURA          **********
'*     CONTOUR PLOTTING - CONTOUR LINES      *
'*     FOR 2D TRIANGLES AND QUADRILATERALS   *
'*     T.R.Chandrupatla and A.D.Belegundu    *
'*********************************************
DefInt I-N
Dim XMIN, XMAX, YMIN, YMAX, XL, XH, YL, YH, XOL, XOH, YOL, YOH
Dim NN, NE, NM, NDIM, NEN, NDN, ND, NL, NMPC
Dim X(), NOC(), FF(), NCON(), XX(), YY(), U(), IC(), ID()
Dim FMIN, FMAX, STP
Dim File1 As String, File2 As String, Dummy As String, Title As String
Const AL = 0.67, NCL = 10
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdPlot_Click()
     Call InputData
     Call FindBoundary
     Call DrawLimits(XMIN, YMIN, XMAX, YMAX)
     Call DrawBoundary
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
     cdlDialog.Filter = "All Files (*.*)|*.*|<Nodal cal file> FileName.n|*.n"
     cdlDialog.FileName = ""
     cdlDialog.ShowOpen
     File2 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, D$: Input #1, Title$
     Line Input #1, D$: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, D$: Input #1, ND, NL, NMPC
     If NDIM <> 2 Or NEN < 3 Or NEN > 4 Then
        picBox.Print "This program supports triangular and quadrilateral"
        picBox.Print "Elements only."
        End
     End If
     ReDim X(NN, NDIM), NOC(NE, NEN), FF(NN), NCON(NE, NEN)
     ReDim XX(3), YY(3), U(3), IC(10), ID(10)
     '=============  COLOR DATA  ===============
     IC(1) = 13: IC(2) = 5: IC(3) = 9: IC(4) = 1: IC(5) = 2
     IC(6) = 10: IC(7) = 14: IC(8) = 6: IC(9) = 4: IC(10) = 12
     For I = 1 To 10: ID(I) = 0: Next I
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, D$
     For I = 1 To NN: Input #1, N: For J = 1 To NDIM
     Input #1, X(N, J): Next J: Next I
     '----- Connectivity
     Line Input #1, D$
     For I = 1 To NE: Input #1, N: For J = 1 To NEN
     Input #1, NOC(N, J): Next J: Line Input #1, D$
     Next I
     Close #1
     Open File2 For Input As #2
     '----- Nodal Values
     Line Input #2, D$
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
Private Sub FindBoundary()
'=============  Find Boundary Lines  ===============
     'Edges defined by nodes in NOC to nodes in NCON
     For IE = 1 To NE
       For I = 1 To NEN
         I1 = I + 1: If I1 > NEN Then I1 = 1
         NCON(IE, I) = NOC(IE, I1)
       Next I
     Next IE
     For IE = 1 To NE
       For I = 1 To NEN
         I1 = NCON(IE, I): I2 = NOC(IE, I)
         INDX = 0
         For JE = IE + 1 To NE
           For J = 1 To NEN
             If NCON(JE, J) <> 0 Then
               If I1 = NCON(JE, J) Or I1 = NOC(JE, J) Then
                 If I2 = NCON(JE, J) Or I2 = NOC(JE, J) Then
                   NCON(JE, J) = 0: INDX = INDX + 1
                 End If
               End If
             End If
           Next J
         Next JE
         If INDX > 0 Then NCON(IE, I) = 0
       Next I
     Next IE
End Sub
Private Sub DrawBoundary()
     picBox.Scale (XL, YH)-(XH, YL)
     picBox.Cls
     '============  Draw Boundary  ==============
     For IE = 1 To NE
       For I = 1 To NEN
         If NCON(IE, I) > 0 Then
           I1 = NCON(IE, I): I2 = NOC(IE, I)
           picBox.Line (X(I1, 1), X(I1, 2))-(X(I2, 1), X(I2, 2))
         End If
       Next I
     Next IE
End Sub
Private Sub DrawContours()
     '===========  Contour Plotting  ===========
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
           XB = 0: YB = 0: UB = 0
           For IT = 1 To NEN
              NIT = NOC(IE, IT)
              XB = XB + 0.25 * X(NIT, 1)
              YB = YB + 0.25 * X(NIT, 2)
              UB = UB + 0.25 * FF(NIT)
           Next IT
           For IT = 1 To NEN
              IT1 = IT + 1: If IT1 > 4 Then IT1 = 1
              XX(1) = XB: YY(1) = YB: U(1) = UB
              NIE = NOC(IE, IT)
              XX(2) = X(NIE, 1): YY(2) = X(NIE, 2): U(2) = FF(NIE)
              NIE = NOC(IE, IT1)
              XX(3) = X(NIE, 1): YY(3) = X(NIE, 2): U(3) = FF(NIE)
              Call ElementPlot
           Next IT
        Else
           Print "NUMBER OF ELEMENT NODES > 4 IS NOT SUPPORTED"
           End
        End If
      Next IE
     For I = 1 To 10: ID(I) = 0: Next I
End Sub
Private Sub ElementPlot()
   'THREE POINTS IN ASCENDING ORDER
        For I = 1 To 2
           C = U(I): II = I
           For J = I + 1 To 3
              If C > U(J) Then
                 C = U(J): II = J
              End If
           Next J
           U(II) = U(I): U(I) = C
           C1 = XX(II): XX(II) = XX(I): XX(I) = C1
           C1 = YY(II): YY(II) = YY(I): YY(I) = C1
        Next I
        SU = (U(1) - FMIN) / STP
        II = Int(SU)
        If II <= SU Then II = II + 1
        UT = FMIN + II * STP
        Do While UT <= U(3)
           ICO = IC(II)
           X1 = ((U(3) - UT) * XX(1) + (UT - U(1)) * XX(3)) / (U(3) - U(1))
           Y1 = ((U(3) - UT) * YY(1) + (UT - U(1)) * YY(3)) / (U(3) - U(1))
           L = 1: If UT > U(2) Then L = 3
           X2 = ((U(L) - UT) * XX(2) + (UT - U(2)) * XX(L)) / (U(L) - U(2))
           Y2 = ((U(L) - UT) * YY(2) + (UT - U(2)) * YY(L)) / (U(L) - U(2))
           picBox.Line (X1, Y1)-(X2, Y2), QBColor(ICO)
           If ID(II) = 0 Then
              picBox.CurrentX = X1
              picBox.CurrentY = Y1
              If (XL < X1 And X1 < XH) And (YL < Y1 And Y1 < YH) Then
                 picBox.Print II
                 ID(II) = 1
              End If
           End If
           UT = UT + STP
           II = II + 1
        Loop
End Sub
Private Sub cmdPrevious_Click()
     T = XL: XL = XOL: XOL = T
     T = YL: YL = YOL: YOL = T
     T = XH: XH = XOH: XOH = T
     T = YH: YH = YOH: YOH = T
     Call DrawBoundary
     Call DrawContours
     cmdOriginal.Enabled = True
End Sub
Private Sub cmdOriginal_Click()
     XL = XMIN: YL = YMIN: XH = XMAX: YH = YMAX
     XOL = XL: YOL = YL: XOH = XH: YOH = YH
     Call DrawBoundary
     Call DrawContours
     cmdOriginal.Enabled = False
     cmdPrevious.Enabled = False
End Sub
Private Sub cmdULeft_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XH = (1 - AL) * XL + AL * XH
     YL = AL * YL + (1 - AL) * YH
     Call DrawBoundary
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdURight_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XL = AL * XL + (1 - AL) * XH
     YL = AL * YL + (1 - AL) * YH
     Call DrawBoundary
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdLLeft_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XH = (1 - AL) * XL + AL * XH
     YH = (1 - AL) * YL + AL * YH
     Call DrawBoundary
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdLRight_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XL = AL * XL + (1 - AL) * XH
     YH = (1 - AL) * YL + AL * YH
     Call DrawBoundary
     Call DrawContours
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
