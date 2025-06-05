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
'*          PROGRAM FEM1D2             *
'*         1-D BAR ELEMENT             *
'*    WITH MULTI-POINT CONSTRAINTS     *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), F(), AREA(), MAT(), DT(), S()
Dim PM(), NU(), U(), MPC(), BT(), Stress(), React()
Dim CNST
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Bandwidth
     Call Stiffness
     Call ModifyForBC
     Call BandSolver
     Call StressCalc
     Call ReactionCalc
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
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     NPR = 2  ' Two properties E and Nu
     ReDim X(NN), NOC(NE, NEN), F(NN), AREA(NE), MAT(NE), DT(NE)
     ReDim PM(NM, NPR), NU(ND), U(ND), MPC(NMPC, 2), BT(NMPC, 3)
     '=============  READ DATA  ===============
     '----- Coordinates -----
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        Input #1, X(N)
     Next I
     '----- Connectivity -----
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
        Input #1, NOC(N, 1), NOC(N, 2)
        Input #1, MAT(N), AREA(N), DT(N)
     Next I
     '----- Specified Displacements -----
     Line Input #1, Dummy
     For I = 1 To ND
        Input #1, NU(I), U(I)
     Next I
     '----- Component Loads -----
     Line Input #1, Dummy
     For I = 1 To NL
        Input #1, N
        Input #1, F(N)
     Next I
     '----- Material Properties -----
     Line Input #1, Dummy
     For I = 1 To NM
        Input #1, N
        For J = 1 To NPR
           Input #1, PM(N, J)
        Next J
     Next I
     '----- Multi-point Constraints B1*Qi+B2*Qj=B0 -----
     If NMPC > 0 Then
        Line Input #1, Dummy
        For I = 1 To NMPC
           Input #1, BT(I, 1), MPC(I, 1), BT(I, 2), MPC(I, 2), BT(I, 3)
        Next I
     End If
     Close #1
End Sub
Private Sub Bandwidth()
     '----- Bandwidth Evaluation -----
     NBW = 0
     For N = 1 To NE
        NABS = Abs(NOC(N, 1) - NOC(N, 2)) + 1
        If NBW < NABS Then NBW = NABS
     Next N
     For I = 1 To NMPC
        NABS = Abs(MPC(I, 1) - MPC(I, 2)) + 1
        If NBW < NABS Then NBW = NABS
     Next I
End Sub
Private Sub Stiffness()
     ReDim S(NN, NBW)
     '----- Stiffness Matrix -----
     For N = 1 To NE
        N1 = NOC(N, 1): N2 = NOC(N, 2): N3 = MAT(N)
        X21 = X(N2) - X(N1): EL = Abs(X21)
        EAL = PM(N3, 1) * AREA(N) / EL
        TL = PM(N3, 1) * PM(N3, 2) * DT(N) * AREA(N) * EL / X21
        '----- Temperature Loads -----
        F(N1) = F(N1) - TL
        F(N2) = F(N2) + TL
        '----- Element Stiffness in Global Locations -----
        S(N1, 1) = S(N1, 1) + EAL
        S(N2, 1) = S(N2, 1) + EAL
        IR = N1: If IR > N2 Then IR = N2
        IC = Abs(N2 - N1) + 1
        S(IR, IC) = S(IR, IC) - EAL
     Next N
End Sub
Private Sub ModifyForBC()
     '----- Decide Penalty Parameter CNST -----
     CNST = 0
     For I = 1 To NN
        If CNST < S(I, 1) Then CNST = S(I, 1)
     Next I
     CNST = CNST * 10000
     '----- Modify for Boundary Conditions -----
        '--- Displacement BC ---
        For I = 1 To ND
           N = NU(I)
           S(N, 1) = S(N, 1) + CNST
           F(N) = F(N) + CNST * U(I)
        Next I
        '--- Multi-point Constraints ---
        For I = 1 To NMPC
           I1 = MPC(I, 1): I2 = MPC(I, 2)
           S(I1, 1) = S(I1, 1) + CNST * BT(I, 1) * BT(I, 1)
           S(I2, 1) = S(I2, 1) + CNST * BT(I, 2) * BT(I, 2)
           IR = I1: If IR > I2 Then IR = I2
           IC = Abs(I2 - I1) + 1
           S(IR, IC) = S(IR, IC) + CNST * BT(I, 1) * BT(I, 2)
           F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
           F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
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
Private Sub StressCalc()
     ReDim Stress(NE)
     '----- Stress Calculation -----
     For N = 1 To NE
        N1 = NOC(N, 1): N2 = NOC(N, 2): N3 = MAT(N)
        EPS = (F(N2) - F(N1)) / (X(N2) - X(N1))
        Stress(N) = PM(N3, 1) * (EPS - PM(N3, 2) * DT(N))
     Next N
End Sub
Private Sub ReactionCalc()
     ReDim React(ND)
     '----- Reaction Calculation -----
     For I = 1 To ND
        N = NU(I)
        React(I) = CNST * (U(I) - F(N))
     Next I
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program FEM1D - CHANDRUPATLA & BELEGUNDU"
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE NO.", "DISPLACEMENT"
     For I = 1 To NN
        Print #2, I, Format(F(I), "0.0000E+00")
     Next I
     '----- Stresses -----
     Print #2, "ELEM NO.", "STRESS"
     For N = 1 To NE
        Print #2, N, Format(Stress(N), "0.0000E+00")
     Next N
     '----- Reactions -----
     Print #2, "NODE NO.", "REACTION"
     For I = 1 To ND
        N = NU(I)
        Print #2, N, Format(React(I), "0.0000E+00")
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

