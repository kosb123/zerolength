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
'***************************************
'*          PROGRAM  CSTKM             *
'*   2-D  CONSTANT STRAIN TRIANGLE     *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), MAT(), PM()
Dim TH(), DT(), NU(), U(), F(), MPC(), BT(), GM()
Dim D(), B(), DB(), SE(), EM(), Q(), STR(), TL(), S()
Dim Stress(), React(), PrinStress(), PltStress()
Dim CNST, LC, NQ, DJ
Dim File1 As String, File2 As String, File3 As String
Dim Title As String, Dummy As String
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
     msg = "1. Plane Stress Analysis" & Chr(13)
     msg = msg + "2. Plane Strain Analysis" & Chr(13)
     msg = msg + "   Your Choice <1 or 2> ?"
     LC = InputBox(msg, "Type of Problem", 1)
     Open File1 For Input As #1
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     NPR = 4   'Material properties E, Nu, Alpha, Density(Rho)
     '----- Total dof is  NQ
     NQ = NDN * NN
     ReDim X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     ReDim TH(NE), DT(NE), NU(ND), U(ND), F(NQ), MPC(NMPC, 2), BT(NMPC, 3)
     ReDim D(3, 3), B(3, 6), DB(3, 6), SE(6, 6), EM(6, 6), Q(6), STR(3), TL(6)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        For J = 1 To NDIM
           Input #1, X(N, J)
        Next J
     Next I
     '----- Connectivity, Material, Thickness, Temp-change
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
        For J = 1 To NEN
           Input #1, NOC(N, J)
        Next J
        Input #1, MAT(N), TH(N), DT(N)
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
      For I = 1 To NE
        NMIN = NOC(I, 1): NMAX = NOC(I, 1)
        For J = 2 To 3
           If NMIN > NOC(I, J) Then NMIN = NOC(I, J)
           If NMAX < NOC(I, J) Then NMAX = NOC(I, J)
        Next J
        NTMP = NDN * (NMAX - NMIN + 1)
        If NBW < NTMP Then NBW = NTMP
     Next I
     For I = 1 To NMPC
        NABS = Abs(MPC(I, 1) - MPC(I, 2)) + 1
        If NBW < NABS Then NBW = NABS
     Next I
     picBox.Print "The Bandwidth is"; NBW
End Sub
Private Sub StiffnMass()
     ReDim S(NQ, NBW), GM(NQ, NBW)
     '----- Global Stiffness Matrix -----
     For N = 1 To NE
     picBox.Print "Forming Stiffness Matrix of Element "; N
     Call DbMat(N, 1)
     Call ElemKM(N)
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
Private Sub ElemKM(N)
     '--- Element Stiffness SE()
     For I = 1 To 6
        For J = 1 To 6
           C = 0
           For K = 1 To 3
              C = C + 0.5 * Abs(DJ) * B(K, I) * DB(K, J) * TH(N)
           Next K
           SE(I, J) = C
        Next J
     Next I
     '--- Element Mass  EM()
     RHO = PM(MAT(N), 4)
     CM = RHO * TH(N) * 0.5 * Abs(DJ) / 12
     For I = 1 To 6: For J = 1 To 6: EM(I, J) = 0: Next J: Next I
     '--- Non-zero elements of mass matrix are defined
     EM(1, 1) = 2 * CM: EM(1, 3) = CM: EM(1, 5) = CM
     EM(2, 2) = 2 * CM: EM(2, 4) = CM: EM(2, 6) = CM
     EM(3, 1) = CM: EM(3, 3) = 2 * CM: EM(3, 5) = CM
     EM(4, 2) = CM: EM(4, 4) = 2 * CM: EM(4, 6) = CM
     EM(5, 1) = CM: EM(5, 3) = CM: EM(5, 5) = 2 * CM
     EM(6, 2) = CM: EM(6, 4) = CM: EM(6, 6) = 2 * CM
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
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     File2 = InputBox("Output File d:\dir\fileName.ext", "Name of File")
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file for Eigen> FileName.inp|*.inp"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
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
     picBox.Print "Run INVITR or JACOBI program to get ";
     picBox.Print "Eigenvalues and Eigenvectors"
End Sub
Private Sub DbMat(N, ISTR)
     '----- D(), B() and DB() matrices
     '--- First the D-Matrix
     M = MAT(N): E = PM(M, 1): PNU = PM(M, 2): AL = PM(M, 3)
     '--- D() Matrix
     If LC = 1 Then
        '--- Plane Stress
        C1 = E / (1 - PNU ^ 2): C2 = C1 * PNU
     Else
        '--- Plane Strain
        C = E / ((1 + PNU) * (1 - 2 * PNU))
        C1 = C * (1 - PNU): C2 = C * PNU
     End If
     C3 = 0.5 * E / (1 + PNU)
     D(1, 1) = C1: D(1, 2) = C2: D(1, 3) = 0
     D(2, 1) = C2: D(2, 2) = C1: D(2, 3) = 0
     D(3, 1) = 0: D(3, 2) = 0: D(3, 3) = C3
     '--- Strain-Displacement Matrix B()
     I1 = NOC(N, 1): I2 = NOC(N, 2): I3 = NOC(N, 3)
     X1 = X(I1, 1): Y1 = X(I1, 2)
     X2 = X(I2, 1): Y2 = X(I2, 2)
     X3 = X(I3, 1): Y3 = X(I3, 2)
     X21 = X2 - X1: X32 = X3 - X2: X13 = X1 - X3
     Y12 = Y1 - Y2: Y23 = Y2 - Y3: Y31 = Y3 - Y1
     DJ = X13 * Y23 - X32 * Y31   'DJ is determinant of Jacobian
     '--- Definition of B() Matrix
     B(1, 1) = Y23 / DJ: B(2, 1) = 0: B(3, 1) = X32 / DJ
     B(1, 2) = 0: B(2, 2) = X32 / DJ: B(3, 2) = Y23 / DJ
     B(1, 3) = Y31 / DJ: B(2, 3) = 0: B(3, 3) = X13 / DJ
     B(1, 4) = 0: B(2, 4) = X13 / DJ: B(3, 4) = Y31 / DJ
     B(1, 5) = Y12 / DJ: B(2, 5) = 0: B(3, 5) = X21 / DJ
     B(1, 6) = 0: B(2, 6) = X21 / DJ: B(3, 6) = Y12 / DJ
     '--- DB Matrix DB = D*B
     For I = 1 To 3
        For J = 1 To 6
           C = 0
           For K = 1 To 3
              C = C + D(I, K) * B(K, J)
           Next K
           DB(I, J) = C
        Next J
     Next I
     If ISTR = 2 Then
     '----- Stress Evaluation
     Q(1) = F(2 * I1 - 1): Q(2) = F(2 * I1)
     Q(3) = F(2 * I2 - 1): Q(4) = F(2 * I2)
     Q(5) = F(2 * I3 - 1): Q(6) = F(2 * I3)
     C1 = AL * DT(N): If LC = 2 Then C1 = C1 * (1 + PNU)
     For I = 1 To 3
        C = 0
        For K = 1 To 6
           C = C + DB(I, K) * Q(K)
        Next K
        STR(I) = C - C1 * (D(I, 1) + D(I, 2))
     Next I
     End If
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

