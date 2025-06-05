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
'****************************************
'*         PROGRAM TRUSSKY              *
'*      TWO DIMENSIONAL TRUSSES         *
'*  T.R.Chandrupatla and A.D.Belegundu  *
'****************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), F(), AREA(), MAT()
Dim PM(), NU(), U(), DT(), SE(), TL(), A()
Dim MPC(), BT(), Stress(), React(), ID()
Dim CNST, NQ, NSUM
Dim Title As String, File1 As String, File2 As String, Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call SkyStore
     Call Stiffness
     Call ModifyBCsky
     Call SkySolve
     Call StressCalc
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
     NPR = 2 ' Two properties E and Alpha
     ReDim X(NN, 2), NOC(NE, NEN), F(NQ), AREA(NE), MAT(NE)
     ReDim PM(NM, NPR), NU(ND), U(ND), DT(NE), SE(4, 4), TL(4)
     ReDim MPC(NMPC, 2), BT(NMPC, 3), ID(NQ)
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
        Input #1, MAT(N), AREA(N), DT(N)
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
Private Sub SkyStore()
     '----- Preparation for Skyline Storage -----
     For I = 1 To NE
        II = NOC(I, 1): If II > NOC(I, 2) Then II = NOC(I, 2)
        II1 = NDN * (II - 1) + 1
        For J = 1 To NEN
           NIJ = NDN * (NOC(I, J) - 1)
           For JJ = 1 To NDN
              NDG = NIJ + JJ: NHT = NDG - II1 + 1
              If NHT > ID(NDG) Then ID(NDG) = NHT
           Next JJ
        Next J
     Next I
     '----- Skyline Height adjustment for MPC
     For I = 1 To NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        NDG = I1: If NDG < I2 Then NDG = I2
        NHT = Abs(I2 - I1) + 1
        If NHT > ID(NDG) Then ID(NDG) = NHT
     Next I
     For I = 2 To NQ
        ID(I) = ID(I) + ID(I - 1)
     Next I
     NSUM = ID(NQ)
End Sub
Private Sub Stiffness()
     ReDim A(NSUM)
     '-----  Stiffness Matrix -----
     For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
        I1 = NOC(N, 1): I2 = NOC(N, 2)
        I3 = MAT(N)
        X21 = X(I2, 1) - X(I1, 1)
        Y21 = X(I2, 2) - X(I1, 2)
        EL = Sqr(X21 * X21 + Y21 * Y21)
        EAL = PM(I3, 1) * AREA(N) / EL
        CS = X21 / EL: SN = Y21 / EL
     '----------- Element Stiffness Matrix SE() -----------
        SE(1, 1) = CS * CS * EAL
        SE(1, 2) = CS * SN * EAL: SE(2, 1) = SE(1, 2)
        SE(1, 3) = -CS * CS * EAL: SE(3, 1) = SE(1, 3)
        SE(1, 4) = -CS * SN * EAL: SE(4, 1) = SE(1, 4)
        SE(2, 2) = SN * SN * EAL
        SE(2, 3) = -CS * SN * EAL: SE(3, 2) = SE(2, 3)
        SE(2, 4) = -SN * SN * EAL: SE(4, 2) = SE(2, 4)
        SE(3, 3) = CS * CS * EAL
        SE(3, 4) = CS * SN * EAL: SE(4, 3) = SE(3, 4)
        SE(4, 4) = SN * SN * EAL
     '-------------- Temperature Load TL() ---------------
        EE0 = PM(I3, 2) * DT(N) * PM(I3, 1) * AREA(N)
        TL(1) = -EE0 * CS: TL(2) = -EE0 * SN
        TL(3) = EE0 * CS: TL(4) = EE0 * SN
        picBox.Print ".... Placing in Global Locations"
        For II = 1 To NEN
           NCT = NDN * (NOC(N, II) - 1)
           For IT = 1 To NDN
              NC = NCT + IT: IID = ID(NC)
              I = NDN * (II - 1) + IT
              For JJ = 1 To NEN
                 NRT = NDN * (NOC(N, JJ) - 1)
                 For JT = 1 To 2
                    J = NDN * (JJ - 1) + JT
                    NR = NRT + JT
                    If NR <= NC Then
                       NLOC = IID - (NC - NR)
                       A(NLOC) = A(NLOC) + SE(I, J)
                    End If
                 Next JT
              Next JJ
              F(NC) = F(NC) + TL(I)
           Next IT
        Next II
     Next N
End Sub
Private Sub ModifyBCsky()
     '----- Decide Penalty Parameter CNST -----
     CNST = 0
     For I = 1 To NQ
        II = ID(I)
        If CNST < A(II) Then CNST = A(II)
     Next I
     CNST = CNST * 10000
     '----- Modify for Boundary Conditions -----
     '--- Displacement BC ---
     For I = 1 To ND
        N = NU(I): II = ID(N)
        A(II) = A(II) + CNST
        F(N) = F(N) + CNST * U(I)
     Next I
        '--- Multi-point Constraints ---
     For I = 1 To NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        A(ID(I1)) = A(ID(I1)) + CNST * BT(I, 1) * BT(I, 1)
        A(ID(I2)) = A(ID(I2)) + CNST * BT(I, 2) * BT(I, 2)
        II = I1: If II < I2 Then II = I2
        IL = ID(II) - Abs(I2 - I1)
        A(IL) = A(IL) + CNST * BT(I, 1) * BT(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
     Next I
End Sub
Private Sub SkySolve()
     '-----  Sky Line Solver -----
     '--- Forward Elimination ---
     For J = 2 To NQ
        NJ = ID(J) - ID(J - 1)
        If NJ <> 1 Then
           K1 = 0: NJ = J - NJ + 1
           For K = NJ To J - 1
              K1 = K1 + 1: KJ = ID(J - 1) + K1: KK = ID(K)
              C = A(KJ) / A(KK)
              For I = K + 1 To J
                 NI = ID(I) - ID(I - 1)
                 If (I - K + 1) <= NI Then
                    IJ = ID(J) - J + I
                    KI = ID(I) - I + K
                    A(IJ) = A(IJ) - C * A(KI)
                 End If
              Next I
           Next K
        End If
     Next J
     For K = 1 To NQ - 1: KK = ID(K)
        C = F(K) / A(KK)
        For I = K + 1 To NQ
           NI = ID(I) - ID(I - 1)
           If (I - K + 1) <= NI Then
              KI = ID(I) - I + K
              F(I) = F(I) - C * A(KI)
           End If
        Next I
     Next K
     '--- Back-substitution ---
     NS = ID(NQ): F(NQ) = F(NQ) / A(NS)
     For I1 = 1 To NQ - 1
        I = NQ - I1
        II = ID(I)
        C = 1 / A(II): F(I) = C * F(I)
        For J = I + 1 To NQ
           J1 = J - I + 1: NJ = ID(J) - ID(J - 1)
           If J1 <= NJ Then
              IJ = ID(J) - J + I
              F(I) = F(I) - C * A(IJ) * F(J)
           End If
        Next J
     Next I1
End Sub
Private Sub StressCalc()
     ReDim Stress(NE)
     '-----  Stress Calculations
     For I = 1 To NE
        I1 = NOC(I, 1)
        I2 = NOC(I, 2)
        I3 = MAT(I)
        X21 = X(I2, 1) - X(I1, 1): Y21 = X(I2, 2) - X(I1, 2)
        EL = Sqr(X21 * X21 + Y21 * Y21)
        CS = X21 / EL
        SN = Y21 / EL
        J2 = 2 * I1
        J1 = J2 - 1
        K2 = 2 * I2
        K1 = K2 - 1
        DLT = (F(K1) - F(J1)) * CS + (F(K2) - F(J2)) * SN
        Stress(I) = PM(I3, 1) * (DLT / EL - PM(I3, 2) * DT(I))
     Next I
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
     Print #2, "Program Trussky - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE#   X-Displ     Y-Displ"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(F(2 * I - 1), "0.0000E+00  ");
        Print #2, Format(F(2 * I), "0.0000E+00")
     Next I
     '----- Stress Output -----
     Print #2, "Elem#    Stress"
     For N = 1 To NE
        Print #2, Format(N, "@@@@@  "); Format(Stress(N), "0.0000E+00")
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

