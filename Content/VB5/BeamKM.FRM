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
'**************************************
'*           PROGRAM BEAMKM           *
'*   STIFFNESS AND MASS GENERATION    *
'* T.R.Chandrupatla and A.D.Belegundu *
'**************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), MAT(), PM()
Dim AREA(), SMI(), EM()
Dim NU(), U(), F(), SE(), MPC(), BT(), S(), GM()
Dim CNST, NQ
Dim Title As String, File1 As String, File2 As String, Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Bandwidth
     Call StiffnMass
     Call ModifyForBC
     Call AddSprMass
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
     NPR = 2  'Material Properties E, Rho (density)
     ReDim X(NN), NOC(NE, NEN), MAT(NE), PM(NM, NPR), SMI(NE)
     ReDim NU(ND), U(ND), F(NQ), SE(4, 4), MPC(NMPC, 2), BT(NMPC, 3)
     ReDim AREA(NE), EM(4, 4)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        Input #1, X(N)
     Next I
     '----- Connectivity, Material, MomentOfInertia, Area
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
       For J = 1 To NEN
           Input #1, NOC(N, J)
        Next J
        Input #1, MAT(N), SMI(N), AREA(N)
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
Private Sub StiffnMass()
     ReDim S(NQ, NBW), GM(NQ, NBW)
     '----- Global Stiffness and Mass Matrices -----
     For N = 1 To NE
        picBox.Print "Forming Stiffness and Mass Matrices of Element "; N
        Call ElemStiffMass(N)
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
                       GM(NR, NC) = GM(NR, NC) + EM(I, J)
                    End If
                 Next JT
              Next JJ
           Next IT
        Next II
     Next N
End Sub
Private Sub ElemStiffMass(N)
'--------  Element Stiffness and Mass Matrices  --------
      N1 = NOC(N, 1)
      N2 = NOC(N, 2)
      M = MAT(NE)
      EL = Abs(X(N1) - X(N2))
      '--- Element Stiffness
      EIL = PM(M, 1) * SMI(N) / EL ^ 3
         SE(1, 1) = 12 * EIL
         SE(1, 2) = EIL * 6 * EL
         SE(1, 3) = -12 * EIL
         SE(1, 4) = EIL * 6 * EL
            SE(2, 1) = SE(1, 2)
            SE(2, 2) = EIL * 4 * EL * EL
            SE(2, 3) = -EIL * 6 * EL
            SE(2, 4) = EIL * 2 * EL * EL
         SE(3, 1) = SE(1, 3)
         SE(3, 2) = SE(2, 3)
         SE(3, 3) = EIL * 12
         SE(3, 4) = -EIL * 6 * EL
            SE(4, 1) = SE(1, 4)
            SE(4, 2) = SE(2, 4)
            SE(4, 3) = SE(3, 4)
            SE(4, 4) = EIL * 4 * EL * EL
      '--- Element Mass
      RHO = PM(M, 2)
      C1 = RHO * AREA(N) * EL / 420
      EM(1, 1) = 156 * C1
      EM(1, 2) = 22 * EL * C1
      EM(1, 3) = 54 * C1
      EM(1, 4) = -13 * EL * C1
         EM(2, 1) = EM(1, 2)
         EM(2, 2) = 4 * EL * EL * C1
         EM(2, 3) = 13 * EL * C1
         EM(2, 4) = -3 * EL * EL * C1
      EM(3, 1) = EM(1, 3)
      EM(3, 2) = EM(2, 3)
      EM(3, 3) = 156 * C1
      EM(3, 4) = -22 * EL * C1
         EM(4, 1) = EM(1, 4)
         EM(4, 2) = EM(2, 4)
         EM(4, 3) = EM(3, 4)
         EM(4, 4) = 4 * EL * EL * C1
End Sub
Private Sub AddSprMass()
'-----  Additional Springs and Lumped Masses  -----
     Do
        N = InputBox("Dof# = 0 Exits this", "DOF# with Spring Support", 0)
        If N = 0 Then Exit Do
        C = InputBox(N, "Support Stiffness at Dof ", 0)
        S(N, 1) = S(N, 1) + C
     Loop
     Do
        N = InputBox("Dof# = 0 Exits this", "DOF# with Lumped Mass", 0)
        If N = 0 Then Exit Do
        C = InputBox(N, "Lumped Mass at", 0)
        GM(N, 1) = GM(N, 1) + C
     Loop
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
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input for Eigen> FileName.inp|*.inp"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Stiffness and Mass for Data in File "; File1
     Print #2, "Num. of DOF    Bandwidth"
     Print #2, NQ, NBW
     '----- Stiffness and Mass Data -----
     Print #2, "Banded Stiffness Matrix"
     For I = 1 To NQ
        For J = 1 To NBW
           Print #2, S(I, J);
        Next J
        Print #2,
     Next I
     Print #2, "Banded Mass Matrix"
     For I = 1 To NQ
        For J = 1 To NBW
           Print #2, GM(I, J);
        Next J
        Print #2,
     Next I
     Print #2, "Starting Vector for Inverse Iteration "
     For I = 1 To NQ
        Print #2, "1 ";
     Next I
     Print #2,
     Close #2
     picBox.Print "Global Stiffness and Mass Matrices are in file "; File2
     picBox.Print "Run INVITR or JACOBI or GENEIGEN program to get ";
     picBox.Print "Eigenvalues and Eigenvectors"
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

