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
'********    PROGRAM FRAME2D      ********
'*      2-D   FRAME ANALYSIS BY FEM      *
'*  T.R.Chandrupatla and A.D.Belegundu   *
'*****************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), MAT(), PM(), ARIN(), React(), S()
Dim NU(), U(), F(), SEP(), SE(), MPC(), UDL()
Dim BT(), DCOS(), ALAMBDA(), ED(), EDP(), EF()
Dim CNST, NQ, ISTF
Dim Title As String, File1 As String, File2 As String, Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Bandwidth
     Call Stiffness
     Call AddLoads
     Call ModifyForBC
     Call BandSolver
     Call EndActions
     Call ReactionCalc
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     Dim msg As String, File1 As String
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     '----- Total dof is  NQ
     NQ = NDN * NN
     NPR = 1   'One Material Property (E)
     ReDim X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR), ARIN(NE, 2)
     ReDim NU(ND), U(ND), F(NQ), SEP(6, 6), SE(6, 6), MPC(NMPC, 2), UDL(NE)
     ReDim BT(NMPC, 3), DCOS(3, 3), ALAMBDA(6, 6), ED(6), EDP(6)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        For J = 1 To NDIM
           Input #1, X(N, J)
        Next J
     Next I
     '----- Connectivity, Material, Area, Mom_Inertia, Distributed Load
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
       For J = 1 To NEN
           Input #1, NOC(N, J)
        Next J
        Input #1, MAT(N), ARIN(N, 1), ARIN(N, 2), UDL(N)
     Next I
     '----- Displacement BC
     Line Input #1, Dummy
     For I = 1 To ND
        Input #1, NU(I), U(I)
     Next I
     '----- Component Loads
     Line Input #1, Dummy
     For I = 1 To NL
        Input #1, N
        Input #1, F(N)
     Next I
     '----- Material Properties
     Line Input #1, Dummy
     For I = 1 To NM
        Input #1, N
        For J = 1 To NPR
           Input #1, PM(N, J)
        Next J
     Next I
     If NMPC > 0 Then
        '-----  Multi-point Constraints
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
        NABS = NDN * (Abs(NOC(N, 1) - NOC(N, 2)) + 1)
        If NBW < NABS Then NBW = NABS
     Next N
     For I = 1 To NMPC
        NABS = Abs(MPC(I, 1) - MPC(I, 2)) + 1
        If NBW < NABS Then NBW = NABS
     Next I
     picBox.Print "The Bandwidth is"; NBW
End Sub
Private Sub Stiffness()
     ReDim S(NQ, NBW)
     '----- Global Stiffness Matrix -----
     For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
        ISTF = 2
        Call Elstif(N)
        picBox.Print ".... Placing in Global Locations"
        For II = 1 To NEN
           NRT = NDN * (NOC(N, II) - 1)
           For IT = 1 To NDN
              NR = NRT + IT
              I = NDN * (II - 1) + IT
              For JJ = 1 To NEN
                 NCT = NDN * (NOC(N, JJ) - 1)
                 For JT = 1 To NDN
                    J = NDN * (JJ - 1) + JT
                    NC = NCT + JT - NR + 1
                    If NC > 0 Then
                       S(NR, NC) = S(NR, NC) + SE(I, J)
                    End If
                 Next JT
              Next JJ
           Next IT
        Next II
     Next N
End Sub
Private Sub Elstif(N)
     '----- Element Stiffness Matrix -----
     I1 = NOC(N, 1): I2 = NOC(N, 2): M = MAT(N)
     X21 = X(I2, 1) - X(I1, 1)
     Y21 = X(I2, 2) - X(I1, 2)
     EL = Sqr(X21 * X21 + Y21 * Y21)
     EAL = PM(M, 1) * ARIN(N, 1) / EL
     EIZL = PM(M, 1) * ARIN(N, 2) / EL
     For I = 1 To 6
     For J = 1 To 6
       SEP(I, J) = 0!
     Next J: Next I
     SEP(1, 1) = EAL: SEP(1, 4) = -EAL: SEP(4, 4) = EAL
     SEP(2, 2) = 12 * EIZL / EL ^ 2: SEP(2, 3) = 6 * EIZL / EL
     SEP(2, 5) = -SEP(2, 2): SEP(2, 6) = SEP(2, 3)
     SEP(3, 3) = 4 * EIZL: SEP(3, 5) = -6 * EIZL / EL: SEP(3, 6) = 2 * EIZL
     SEP(5, 5) = 12 * EIZL / EL ^ 2: SEP(5, 6) = -6 * EIZL / EL
     SEP(6, 6) = 4 * EIZL
     For I = 1 To 6
     For J = I To 6
       SEP(J, I) = SEP(I, J)
     Next J: Next I
'CONVERT ELEMENT STIFFNESS MATRIX TO GLOBAL SYSTEM
     DCOS(1, 1) = X21 / EL: DCOS(1, 2) = Y21 / EL: DCOS(1, 3) = 0
     DCOS(2, 1) = -DCOS(1, 2): DCOS(2, 2) = DCOS(1, 1): DCOS(2, 3) = 0
     DCOS(3, 1) = 0: DCOS(3, 2) = 0: DCOS(3, 3) = 1
     For I = 1 To 6
     For J = 1 To 6
       ALAMBDA(I, J) = 0!
     Next J: Next I
     For K = 1 To 2
       IK = 3 * (K - 1)
       For I = 1 To 3
       For J = 1 To 3
         ALAMBDA(I + IK, J + IK) = DCOS(I, J)
       Next J: Next I
     Next K
     If ISTF = 1 Then Exit Sub
     For I = 1 To 6
     For J = 1 To 6
       SE(I, J) = 0
       For K = 1 To 6
         SE(I, J) = SE(I, J) + SEP(I, K) * ALAMBDA(K, J)
       Next K
     Next J: Next I
     For I = 1 To 6: For J = 1 To 6: SEP(I, J) = SE(I, J): Next J: Next I
     For I = 1 To 6
     For J = 1 To 6
       SE(I, J) = 0
       For K = 1 To 6
         SE(I, J) = SE(I, J) + ALAMBDA(K, I) * SEP(K, J)
       Next K
     Next J: Next I
End Sub
Private Sub AddLoads()
'----- Loads due to uniformly distributed load on element
     For N = 1 To NE
     If Abs(UDL(N)) > 0 Then
       ISTF = 1
       Call Elstif(N)
       I1 = NOC(N, 1): I2 = NOC(N, 2)
       X21 = X(I2, 1) - X(I1, 1)
       Y21 = X(I2, 2) - X(I1, 2)
       EL = Sqr(X21 * X21 + Y21 * Y21)
       ED(1) = 0: ED(4) = 0
       ED(2) = UDL(N) * EL / 2: ED(5) = ED(2)
       ED(3) = UDL(N) * EL ^ 2 / 12: ED(6) = -ED(3)
       For I = 1 To 6
         EDP(I) = 0
         For K = 1 To 6
           EDP(I) = EDP(I) + ALAMBDA(K, I) * ED(K)
         Next K
       Next I
       For I = 1 To 3
         F(3 * I1 - 3 + I) = F(3 * I1 - 3 + I) + EDP(I)
         F(3 * I2 - 3 + I) = F(3 * I2 - 3 + I) + EDP(I + 3)
       Next I
     End If
     Next N
End Sub
Private Sub ModifyForBC()
     '----- Decide Penalty Parameter CNST -----
     CNST = 0
     For I = 1 To NQ
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
     '----- Band Solver -----
     N1 = NQ - 1
     '--- Forward Elimination
     For K = 1 To N1
        NK = NQ - K + 1
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
     F(NQ) = F(NQ) / S(NQ, 1)
     For KK = 1 To N1
        K = NQ - KK
        C1 = 1 / S(K, 1)
        F(K) = C1 * F(K)
        NK = NQ - K + 1
        If NK > NBW Then NK = NBW
        For J = 2 To NK
           F(K) = F(K) - C1 * S(K, J) * F(K + J - 1)
        Next J
     Next KK
End Sub
Private Sub EndActions()
     ReDim EF(NE, 6)
     '----- Calculating Member End-Forces
     For N = 1 To NE
       ISTF = 1
       Call Elstif(N)
       I1 = NOC(N, 1): I2 = NOC(N, 2)
       X21 = X(I2, 1) - X(I1, 1)
       Y21 = X(I2, 2) - X(I1, 2)
       EL = Sqr(X21 * X21 + Y21 * Y21)
       For I = 1 To 3
         ED(I) = F(3 * I1 - 3 + I): ED(I + 3) = F(3 * I2 - 3 + I)
       Next I
       For I = 1 To 6
         EDP(I) = 0
         For K = 1 To 6
           EDP(I) = EDP(I) + ALAMBDA(I, K) * ED(K)
         Next K
       Next I
'----- END FORCES DUE TO DISTRIBUTED LOADS
       If Abs(UDL(N)) > 0 Then
         ED(1) = 0: ED(4) = 0
         ED(2) = -UDL(N) * EL / 2: ED(5) = ED(2)
         ED(3) = -UDL(N) * EL ^ 2 / 12: ED(6) = -ED(3)
       Else
         For K = 1 To 6: ED(K) = 0: Next K
       End If
       For I = 1 To 6
         EF(N, I) = ED(I)
         For K = 1 To 6
           EF(N, I) = EF(N, I) + SEP(I, K) * EDP(K)
         Next K
       Next I
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
     Print #2, "Program Frame2D - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE#   X-Displ     Y-Displ    Z-Rotation"
     For I = 1 To NN
        I1 = 3 * I - 2: I2 = I1 + 1: I3 = I1 + 2
        Print #2, Format(I, "@@@@@  "); Format(F(I1), "0.0000E+00  ");
        Print #2, Format(F(I2), "0.0000E+00  "); Format(F(I3), "0.0000E+00")
     Next I
     '----- Member End Actions -----
     Print #2, "Member End-Forces"
     For N = 1 To NE
        Print #2, "Member# "; N
        For I = 1 To 2
          II = (I - 1) * 3
          Print #2, Format(EF(N, II + 1), "0.0000E+00  ");
          Print #2, Format(EF(N, II + 2), "0.0000E+00  ");
          Print #2, Format(EF(N, II + 3), "0.0000E+00")
        Next I
     Next N
     '----- Reactions -----
     Print #2, "DOF#", "Reaction"
     For I = 1 To ND
        N = NU(I)
        Print #2, Format(N, "@@@@@  "); Format(React(I), "0.0000E+00")
     Next I
     Close #2
     picBox.Print "Complete results are in file "; File2
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

