VERSION 5.00
Begin VB.Form frmFemCB 
   Caption         =   "FemCB"
   ClientHeight    =   6228
   ClientLeft      =   60
   ClientTop       =   348
   ClientWidth     =   8616
   LinkTopic       =   "Form1"
   ScaleHeight     =   6228
   ScaleWidth      =   8616
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
'***************************************
'*          PROGRAM  HEAT1D            *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NE, NN, NBC, NQ, NBW
Dim X(), S(), TC(), F(), V(), H(), NB()
Dim BC() As String
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Stiffness
     Call ModifyForBC
     Call BandSolver
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy
     Input #1, NE, NBC, NQ
     NN = NE + 1
     NBW = 2    'NBW is half the bandwidth
     ReDim X(NN), S(NN, NBW), TC(NE), F(NN), V(NBC), H(NBC), NB(NBC)
     ReDim BC(NBC)
     '=============  READ DATA  ===============
     '----- Element Thermal Conductivity -----
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
        Input #1, TC(N)
     Next I
     '----- Coordinates -----
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        Input #1, X(N)
     Next I
     '----- Boundary Conditions -----
     Line Input #1, Dummy
     For I = 1 To NBC
        Input #1, NB(I), BC$(I)
        If BC(I) = "TEMP" Or BC(I) = "temp" Then Input #1, V(I)
        If BC(I) = "FLUX" Or BC(I) = "flux" Then Input #1, V(I)
        If BC(I) = "CONV" Or BC(I) = "conv" Then Input #1, H(I), V(I)
     Next I
     '----- Calculate and Input Nodal Heat Source Vector -----
     Line Input #1, Dummy
     For I = 1 To NN: F(I) = 0: Next I
     For I = 1 To NQ
        Input #1, N
        Input #1, F(N)
     Next I
     Close #1
End Sub
Private Sub Stiffness()
     ReDim S(NN, NBW)
     '----- Stiffness Matrix -----
     For J = 1 To NBW
     For I = 1 To NN: S(I, J) = 0: Next I: Next J
     For I = 1 To NE
     I1 = I: I2 = I + 1
     ELL = Abs(X(I2) - X(I1))
     EKL = TC(I) / ELL
     S(I1, 1) = S(I1, 1) + EKL
     S(I2, 1) = S(I2, 1) + EKL
     S(I1, 2) = S(I1, 2) - EKL: Next I
End Sub
Private Sub ModifyForBC()
     '----- Decide Penalty Parameter CNST -----
     AMAX = 0
     For I = 1 To NN
        If S(I, 1) > AMAX Then AMAX = S(I, 1)
     Next I
     CNST = AMAX * 10000
     For I = 1 To NBC
        N = NB(I)
        If BC(I) = "CONV" Or BC(I) = "conv" Then
           S(N, 1) = S(N, 1) + H(I)
           F(N) = F(N) + H(I) * V(I)
        ElseIf BC(I) = "FLUX" Or BC(I) = "flux" Then
           F(N) = F(N) - V(I)
        Else
           S(N, 1) = S(N, 1) + CNST
           F(N) = F(N) + CNST * V(I)
        End If
     Next I
End Sub
Private Sub BandSolver()
     '----- Equation Solving using Band Solver -----
          N = NN
     '----- Forward Elimination -----
     For K = 1 To N - 1
        NBK = N - K + 1
        If N - K + 1 > NBW Then NBK = NBW
        For I = K + 1 To NBK + K - 1
           I1 = I - K + 1
           C = S(K, I1) / S(K, 1)
           For J = I To NBK + K - 1
              J1 = J - I + 1
              J2 = J - K + 1
              S(I, J1) = S(I, J1) - C * S(K, J2)
           Next J
           F(I) = F(I) - C * F(K)
        Next I
     Next K
     '----- Back Substitution -----
     F(N) = F(N) / S(N, 1)
     For II = 1 To N - 1
        I = N - II
        NBI = N - I + 1
        If N - I + 1 > NBW Then NBI = NBW
        Sum = 0!
        For J = 2 To NBI
           Sum = Sum + S(I, J) * F(I + J - 1)
        Next J
        F(I) = (F(I) - Sum) / S(I, 1)
     Next II
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Heat1D - CHANDRUPATLA & BELEGUNDU"
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE#   TEMPERATURE"
     For I = 1 To NN
        Print #2, I, Format(F(I), "0.0000E+00")
     Next I
     Close #2
     picBox.Print "RESULTS ARE IN FILE "; File2
End Sub

Private Sub cmdView_Click()
   Dim ALine As String, CRLF As String, File1 As String
   CRLF = Chr$(13) + Chr$(10)
   picBox.Visible = False
   txtView.Visible = True
   txtView.Text = ""
   Open File2 For Input As #1
   Do While Not EOF(1)
     Line Input #1, ALine
     txtView.Text = txtView.Text + ALine + CRLF
   Loop
   Close #1
End Sub

