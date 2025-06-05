VERSION 5.00
Begin VB.Form frmFemCB 
   Caption         =   "FemCB"
   ClientHeight    =   6168
   ClientLeft      =   60
   ClientTop       =   348
   ClientWidth     =   8712
   LinkTopic       =   "Form1"
   ScaleHeight     =   6168
   ScaleWidth      =   8712
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
   Begin VB.CommandButton cmdView 
      Caption         =   "View Results"
      Enabled         =   0   'False
      Height          =   375
      Left            =   3720
      TabIndex        =   4
      Top             =   5640
      Width           =   1212
   End
   Begin VB.TextBox txtView 
      Height          =   4932
      Left            =   600
      MultiLine       =   -1  'True
      ScrollBars      =   3  'Both
      TabIndex        =   3
      Top             =   480
      Visible         =   0   'False
      Width           =   7572
   End
   Begin VB.CommandButton cmdEnd 
      Caption         =   "END"
      Height          =   375
      Left            =   6120
      TabIndex        =   2
      Top             =   5640
      Width           =   972
   End
   Begin VB.CommandButton cmdStart 
      Caption         =   "START"
      Height          =   375
      Left            =   1680
      TabIndex        =   1
      Top             =   5640
      Width           =   1095
   End
   Begin VB.PictureBox picBox 
      Height          =   5292
      Left            =   360
      ScaleHeight     =   5244
      ScaleWidth      =   7884
      TabIndex        =   0
      Top             =   240
      Width           =   7932
   End
End
Attribute VB_Name = "frmFemCB"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
'*****      PROGRAM BESTFIT        *****
'*          BEST FIT PROGRAM           *
'*        FOR 3-NODED TRIANGLES        *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NMPC, NQ, NBW
Dim X(), NOC(), FS(), S(), F(), SE(), FE()
Dim File1 As String, File2 As String, File3 As String
Dim Title As String, Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Bandwidth
     Call Stiffness
     Call BandSolver
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     Dim msg As String
     cdlDialog.Filter = "All Files (*.*)|*.*|<Mesh file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     cdlDialog.Filter = "All Files (*.*)|*.*|<Element values> FileName.e|*.e"
     cdlDialog.ShowOpen
     File2 = cdlDialog.FileName
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     If NDIM <> 2 Or NEN <> 3 Then
        picBox.Print "This program supports 3 Noded Triangles only"
        Close #1
        End
     End If
     Open File2 For Input As #2
     ReDim X(NN, 2), NOC(NE, 3), FS(NE), SE(3, 3), FE(3)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        For J = 1 To NDIM
           Input #1, X(N, J)
        Next J
     Next I
     '----- Connectivity
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
        For J = 1 To NEN
           Input #1, NOC(N, J)
        Next J
        Line Input #1, Dummy
     Next I
     Close #1
     '--- Element Function Values
     Line Input #2, Dummy
     For I = 1 To NE
         Input #2, FS(I)
     Next I
     Close #2
End Sub
Private Sub Bandwidth()
     '----- Bandwidth Evaluation -----
     NQ = NN: NBW = 0
     For I = 1 To NE
        NMIN = NOC(I, 1): NMAX = NOC(I, 1)
        For J = 2 To 3
           If NMIN > NOC(I, J) Then NMIN = NOC(I, J)
           If NMAX < NOC(I, J) Then NMAX = NOC(I, J)
        Next J
        NTMP = NMAX - NMIN + 1
        If NBW < NTMP Then NBW = NTMP
     Next I
     picBox.Print "The Bandwidth is"; NBW
End Sub
Private Sub Stiffness()
     ReDim S(NQ, NBW), F(NQ)
     '---  Global Stiffness Matrix
     For N = 1 To NE
        Call ElemStiff(N)
        For II = 1 To 3
           NR = NOC(N, II): F(NR) = F(NR) + FE(II)
           For JJ = 1 To 3
              NC = NOC(N, JJ) - NR + 1
              If NC > 0 Then
                 S(NR, NC) = S(NR, NC) + SE(II, JJ)
              End If
           Next JJ
        Next II
     Next N
     picBox.Print "Stiffness Formation completed..."
End Sub
Private Sub ElemStiff(N)
     '--- Element Stiffness Formation
     I1 = NOC(N, 1): I2 = NOC(N, 2): I3 = NOC(N, 3)
     X1 = X(I1, 1): Y1 = X(I1, 2)
     X2 = X(I2, 1): Y2 = X(I2, 2)
     X3 = X(I3, 1): Y3 = X(I3, 2)
     X21 = X2 - X1: X32 = X3 - X2: X13 = X1 - X3
     Y12 = Y1 - Y2: Y23 = Y2 - Y3: Y31 = Y3 - Y1
     DJ = X13 * Y23 - X32 * Y31      'DETERMINANT OF JACOBIAN
     AE = Abs(DJ) / 24
     SE(1, 1) = 2 * AE: SE(1, 2) = AE: SE(1, 3) = AE
     SE(2, 1) = AE: SE(2, 2) = 2 * AE: SE(2, 3) = AE
     SE(3, 1) = AE: SE(3, 2) = AE: SE(3, 3) = 2 * AE
     A1 = FS(N) * Abs(DJ) / 6
     FE(1) = A1: FE(2) = A1: FE(3) = A1
End Sub
Private Sub BandSolver()
     '----- Band Solver -----
     N1 = NQ - 1
     '--- Forward Elimination
     For K = 1 To N1
        NK = NN - K + 1
        If NK > NBW Then NK = NBW
           For I = 2 To NK
              C1 = S(K, I) / S(K, 1)
              I1 = K + I - 1
              For J = I To NK
                 J1 = J - I + 1
                 S(I1, J1) = S(I1, J1) - C1 * S(K, J)
              Next J
              F(I1) = F(I1) - C1 * F(K)
           Next I
        Next K
     '--- Back-substitution
     F(NN) = F(NN) / S(NN, 1)
     For KK = 1 To N1
        K = NN - KK
        C1 = 1 / S(K, 1)
        F(K) = C1 * F(K)
        NK = NN - K + 1
        If NK > NBW Then NK = NBW
        For J = 2 To NK
           F(K) = F(K) - C1 * S(K, J) * F(K + J - 1)
        Next J
     Next KK
End Sub
Private Sub Output()
     '===== Print Nodal Values of Stresses
     cdlDialog.Filter = "All Files (*.*)|*.*|<Nodal values> FileName.n|*.n"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File3 = cdlDialog.FileName
     Open File3 For Output As #3
     Print #3, "Nodal Values for Data in Files "; File1; " & "; File2
     For I = 1 To NN
        Print #3, F(I)
     Next I
     Close #3
     picBox.Print "Complete results are in file "; File3
     picBox.Print "Run CONTOURA or CONTOURB program"
End Sub
Private Sub cmdView_Click()
   Dim ALine As String, CRLF As String
   CRLF = Chr$(13) + Chr$(10)
   picBox.Visible = False
   txtView.Visible = True
   txtView.Text = ""
   Open File3 For Input As #1
   Do While Not EOF(1)
     Line Input #1, ALine
     txtView.Text = txtView.Text + ALine + CRLF
   Loop
   Close #1
End Sub

