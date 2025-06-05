VERSION 5.00
Begin VB.Form frmPlot2D 
   Caption         =   "Plot2D"
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
      TabIndex        =   10
      Top             =   2400
      Width           =   855
   End
   Begin VB.CommandButton cmdLRight 
      Caption         =   "ZoomLR"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7560
      TabIndex        =   9
      Top             =   2400
      Width           =   855
   End
   Begin VB.CommandButton cmdOriginal 
      Caption         =   "Original"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7560
      TabIndex        =   8
      Top             =   3360
      Width           =   855
   End
   Begin VB.CommandButton cmdPrevious 
      Caption         =   "Pevious"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6600
      TabIndex        =   7
      Top             =   3360
      Width           =   855
   End
   Begin VB.CommandButton cmdURight 
      Caption         =   "ZoomUR"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7560
      TabIndex        =   6
      Top             =   1800
      Width           =   855
   End
   Begin VB.CommandButton cmdULeft 
      Caption         =   "ZoomUL"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6600
      TabIndex        =   5
      Top             =   1800
      Width           =   855
   End
   Begin VB.CommandButton cmdEnd 
      Caption         =   "END"
      Height          =   375
      Left            =   7680
      TabIndex        =   4
      Top             =   5880
      Width           =   975
   End
   Begin VB.CommandButton cmdNonode 
      Caption         =   "Erase  Node#"
      Enabled         =   0   'False
      Height          =   495
      Left            =   7680
      TabIndex        =   3
      Top             =   5280
      Width           =   975
   End
   Begin VB.CommandButton cmdNode 
      Caption         =   "Show   Node#"
      Enabled         =   0   'False
      Height          =   495
      Left            =   6600
      TabIndex        =   2
      Top             =   5280
      Width           =   975
   End
   Begin VB.CommandButton cmdPlot 
      Caption         =   "START"
      Height          =   375
      Left            =   6600
      TabIndex        =   1
      Top             =   5880
      Width           =   975
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
Attribute VB_Name = "frmPlot2D"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
'************************************************
'*               PROGRAM PLOT2D                 *
'*     PLOTS 2D MESHES - TRIANGLES AND QUADS    *
'*     (c) T.R.CHANDRUPATLA & A.D.BELEGUNDU     *
'************************************************
DefInt I-N
Dim NN, NE, NM, NDIM, NEN, NDN, ND, NL, NMPC
Dim X(), NOC()
Dim XMIN, XMAX, YMIN, YMAX, XL, XH, YL, YH, XOL, XOH, YOL, YOH
Dim File1 As String, Dummy As String, Title As String
Const AL = 0.67
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdNode_Click()
    '============  Show Node Numbers  ==============
     For I = 1 To NN
        picBox.CurrentX = X(I, 1)
        picBox.CurrentY = X(I, 2)
        picBox.Print I
     Next I
     cmdNode.Enabled = False
     cmdNonode.Enabled = True
End Sub
Private Sub cmdNonode_Click()
     Call DrawElements
     cmdNonode.Enabled = False
End Sub
Private Sub cmdPlot_Click()
     Call InputData
     Call DrawLimits(XMIN, YMIN, XMAX, YMAX)
     Call DrawElements
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
     Open File1 For Input As #1
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     If NDIM <> 2 Then
        picBox.Print "THE PROGRAM SUPPORTS TWO DIMENSIONAL PLOTS ONLY"
        picBox.Print "THE DIMENSION OF THE DATA IS  "; NDIM
        End
     End If
     ReDim X(NN, NDIM), NOC(NE, NEN)
     '=============  READ DATA  ===============
     Line Input #1, Dummy
     For I = 1 To NN: Input #1, N: For J = 1 To NDIM
     Input #1, X(N, J): Next J: Next I
     Line Input #1, Dummy
     For I = 1 To NE: Input #1, N: For J = 1 To NEN
     Input #1, NOC(N, J): Next J: Line Input #1, Dummy
     Next I
     Close #1
End Sub
Private Sub DrawLimits(XMIN, YMIN, XMAX, YMAX)
     XMAX = X(1, 1): YMAX = X(1, 2): XMIN = X(1, 1): YMIN = X(1, 2)
     For I = 2 To NN
        If XMAX < X(I, 1) Then XMAX = X(I, 1)
        If YMAX < X(I, 2) Then YMAX = X(I, 2)
        If XMIN > X(I, 1) Then XMIN = X(I, 1)
        If YMIN > X(I, 2) Then YMIN = X(I, 2)
     Next I
     XL = (XMAX - XMIN): YL = (YMAX - YMIN)
     A = XL: If A < YL Then A = YL
     XB = 0.5 * (XMIN + XMAX)
     YB = 0.5 * (YMIN + YMAX)
     XMIN = XB - 0.55 * A: XMAX = XB + 0.55 * A
     YMIN = YB - 0.55 * A: YMAX = YB + 0.55 * A
     XL = XMIN: YL = YMIN: XH = XMAX: YH = YMAX
     XOL = XL: YOL = YL: XOH = XH: YOH = YH
End Sub
Private Sub DrawElements()
     '===========  Draw Elements  ================
     picBox.Scale (XL, YH)-(XH, YL)
     picBox.Cls
     For IE = 1 To NE
       For II = 1 To NEN
        I2 = II + 1
        If II = NEN Then I2 = 1
        X1 = X(NOC(IE, II), 1): Y1 = X(NOC(IE, II), 2)
        X2 = X(NOC(IE, I2), 1): Y2 = X(NOC(IE, I2), 2)
        picBox.Line (X1, Y1)-(X2, Y2), QBColor(1)
        If NEN = 2 Then Exit For
       Next II
     Next IE
     cmdNode.Enabled = True
End Sub
Private Sub cmdPrevious_Click()
     T = XL: XL = XOL: XOL = T
     T = YL: YL = YOL: YOL = T
     T = XH: XH = XOH: XOH = T
     T = YH: YH = YOH: YOH = T
     Call DrawElements
     cmdNonode.Enabled = False
     cmdOriginal.Enabled = True
End Sub
Private Sub cmdOriginal_Click()
     XL = XMIN: YL = YMIN: XH = XMAX: YH = YMAX
     XOL = XL: YOL = YL: XOH = XH: YOH = YH
     Call DrawElements
     cmdNonode.Enabled = False
     cmdOriginal.Enabled = False
     cmdPrevious.Enabled = False
End Sub
Private Sub cmdULeft_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XH = (1 - AL) * XL + AL * XH
     YL = AL * YL + (1 - AL) * YH
     Call DrawElements
     cmdNonode.Enabled = False
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdURight_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XL = AL * XL + (1 - AL) * XH
     YL = AL * YL + (1 - AL) * YH
     Call DrawElements
     cmdNonode.Enabled = False
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdLLeft_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XH = (1 - AL) * XL + AL * XH
     YH = (1 - AL) * YL + AL * YH
     Call DrawElements
     cmdNonode.Enabled = False
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
Private Sub cmdLRight_Click()
     XOL = XL: XOH = XH: YOL = YL: YOH = YH
     XL = AL * XL + (1 - AL) * XH
     YH = (1 - AL) * YL + AL * YH
     Call DrawElements
     cmdNonode.Enabled = False
     cmdOriginal.Enabled = True
     cmdPrevious.Enabled = True
End Sub
