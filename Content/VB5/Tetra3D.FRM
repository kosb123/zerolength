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
'*******       PROGRAM TETRA3D       *******
'*    Three Dimensional Stress Analysis    *
'*          Tetrahedral Elements           *
'*   T.R.Chandrupatla and A.D.Belegundu    *
'*******************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), F(), NOC(), MAT(), PM(), S()
Dim NU(), U(), MPC(), BT(), QT()
Dim D(), B(), DB(), SE(), STR(), DT()
Dim Stress(), React(), PrinStress(), PltStress()
Dim CNST, NQ, DJ
Dim File1 As String, File2 As String, File3 As String
Dim Title As String, Dummy As String
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
     NPR = 3   'Material Properties  E, Nu, Alpha
     ReDim X(NN, NDIM), F(NQ), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     ReDim NU(ND), U(ND), MPC(NMPC, 2), BT(NMPC, 3), QT(12)
     ReDim D(6, 6), B(6, 12), DB(6, 12), SE(12, 12), STR(6), DT(NE)
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
        Input #1, MAT(N), DT(N)
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
        For J = 2 To NEN
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
Private Sub Stiffness()
     ReDim S(NQ, NBW)
     '----- Global Stiffness Matrix -----
          For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
        Call DbMat(N, 1)
     '--- Element Stiffness
        For I = 1 To 12
        For J = 1 To 12
        SE(I, J) = 0
        For K = 1 To 6
           SE(I, J) = SE(I, J) + B(K, I) * DB(K, J) * Abs(DJ) / 6
        Next K: Next J: Next I
     '--- Temperature Load Vector QT()
        C = AL * DT(N)
        For I = 1 To 12
           DSUM = DB(1, I) + DB(2, I) + DB(3, I)
           QT(I) = C * Abs(DJ) * DSUM / 6
        Next I
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
              F(NR) = F(NR) + QT(I)
           Next IT
        Next II
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
Private Sub StressCalc()
     ReDim Stress(NE, 6), PrinStress(NE, 3)
'--- Stress Calculations
     PI = 3.141593
     For N = 1 To NE
        Call DbMat(N, 2)
        '--- Principal Stress Calculations
        AI1 = STR(1) + STR(2) + STR(3)
        AI21 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1)
        AI22 = STR(4) * STR(4) + STR(5) * STR(5) + STR(6) * STR(6)
        AI2 = AI21 - AI22
        AI31 = STR(1) * STR(2) * STR(3) + 2 * STR(4) * STR(5) * STR(6)
        AI32 = STR(1) * STR(4) ^ 2 + STR(2) * STR(5) ^ 2 + STR(3) * STR(6) ^ 2
        AI3 = AI31 - AI32
        C1 = AI2 - AI1 ^ 2 / 3
        C2 = -2 * (AI1 / 3) ^ 3 + AI1 * AI2 / 3 - AI3
        C3 = 2 * Sqr(-C1 / 3)
        TH = -3 * C2 / (C1 * C3): TH2 = Abs(1 - TH * TH)
        If TH = 0 Then TH = PI / 2
        If TH > 0 Then TH = Atn(Sqr(TH2) / TH)
        If TH < 0 Then TH = PI - Atn(Sqr(TH2) / TH)
        TH = TH / 3
        '--- Principal Stresses
        P1 = AI1 / 3 + C3 * Cos(TH)
        P2 = AI1 / 3 + C3 * Cos(TH + 2 * PI / 3)
        P3 = AI1 / 3 + C3 * Cos(TH + 4 * PI / 3)
        For I = 1 To 6
           Stress(N, I) = STR(I)
        Next I
        PrinStress(N, 1) = P1
        PrinStress(N, 2) = P2
        PrinStress(N, 3) = P3
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
     Print #2, "Program Tetra3D - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE#   X-Displ     Y-Displ     Z-Displ"
     For I = 1 To NN
        II = 3 * (I - 1)
        Print #2, Format(I, "@@@@@  "); Format(F(II + 1), "0.0000E+00  ");
        Print #2, Format(F(II + 2), "0.0000E+00  "); Format(F(II + 3), "0.0000E+00")
     Next I
     '----- Stress Output -----
     For N = 1 To NE
        Print #2, "Stresses in Element Number "; N
        Print #2, "   Normal Stresses SX, SY, SZ"
        Print #2, "   "; Format(Stress(N, 1), "0.0000E+00  ");
        Print #2, Format(Stress(N, 2), "0.0000E+00  ");
        Print #2, Format(Stress(N, 3), "0.0000E+00  ")
        Print #2, "  Shear Stresses TYZ, TXZ, TXY"
        Print #2, "   "; Format(Stress(N, 4), "0.0000E+00  ");
        Print #2, Format(Stress(N, 5), "0.0000E+00  ");
        Print #2, Format(Stress(N, 6), "0.0000E+00")
        Print #2, "   Principal Stresses"
        Print #2, "   "; Format(PrinStress(N, 1), "0.0000E+00  ");
        Print #2, Format(PrinStress(N, 2), "0.0000E+00  ");
        Print #2, Format(PrinStress(N, 3), "0.0000E+00")
     Next N
     '----- Reactions -----
     Print #2, "DOF#  REACTION"
     For I = 1 To ND
        N = NU(I)
        Print #2, Format(N, "@@@@@  "); Format(React(I), "0.0000E+00")
     Next I
     Close #2
     picBox.Print "Complete results are in file "; File2
End Sub
Private Sub DbMat(N, ISTR)
     '----- D(), B() and DB() matrices
     '--- First the D-Matrix
     M = MAT(N): E = PM(M, 1): PNU = PM(M, 2): AL = PM(M, 3)
     C4 = E / ((1 + PNU) * (1 - 2 * PNU))
     C1 = C4 * (1 - PNU): C2 = C4 * PNU: C3 = 0.5 * E / (1 + PNU)
     For I = 1 To 6: For J = 1 To 6: D(I, J) = 0: Next J: Next I
     D(1, 1) = C1: D(1, 2) = C2: D(1, 3) = C2
     D(2, 1) = C2: D(2, 2) = C1: D(2, 3) = C2
     D(3, 1) = C2: D(3, 2) = C2: D(3, 3) = C1
     D(4, 4) = C3: D(5, 5) = C3: D(6, 6) = C3
     '--- Strain-Displacement Matrix B()
     I1 = NOC(N, 1): I2 = NOC(N, 2): I3 = NOC(N, 3): I4 = NOC(N, 4)
     X14 = X(I1, 1) - X(I4, 1): X24 = X(I2, 1) - X(I4, 1): X34 = X(I3, 1) - X(I4, 1)
     Y14 = X(I1, 2) - X(I4, 2): Y24 = X(I2, 2) - X(I4, 2): Y34 = X(I3, 2) - X(I4, 2)
     Z14 = X(I1, 3) - X(I4, 3): Z24 = X(I2, 3) - X(I4, 3): Z34 = X(I3, 3) - X(I4, 3)
     DJ1 = X14 * (Y24 * Z34 - Z24 * Y34)
     DJ2 = Y14 * (Z24 * X34 - X24 * Z34)
     DJ3 = Z14 * (X24 * Y34 - Y24 * X34)
     DJ = DJ1 + DJ2 + DJ3
     A11 = (Y24 * Z34 - Z24 * Y34) / DJ
     A21 = (Z24 * X34 - X24 * Z34) / DJ
     A31 = (X24 * Y34 - Y24 * X34) / DJ
     A12 = (Y34 * Z14 - Z34 * Y14) / DJ
     A22 = (Z34 * X14 - X34 * Z14) / DJ
     A32 = (X34 * Y14 - Y34 * X14) / DJ
     A13 = (Y14 * Z24 - Z14 * Y24) / DJ
     A23 = (Z14 * X24 - X14 * Z24) / DJ
     A33 = (X14 * Y24 - Y14 * X24) / DJ
     '---  B Matrix
     For I = 1 To 6: For J = 1 To 12: B(I, J) = 0: Next J: Next I
     B(1, 1) = A11: B(1, 4) = A12: B(1, 7) = A13: B(1, 10) = -A11 - A12 - A13
     B(2, 2) = A21: B(2, 5) = A22: B(2, 8) = A23: B(2, 11) = -A21 - A22 - A23
     B(3, 3) = A31: B(3, 6) = A32: B(3, 9) = A33: B(3, 12) = -A31 - A32 - A33
     B(4, 2) = A31: B(4, 3) = A21: B(4, 5) = A32: B(4, 6) = A22: B(4, 8) = A33
                B(4, 9) = A23: B(4, 11) = B(3, 12): B(4, 12) = B(2, 11)
     B(5, 1) = A31: B(5, 3) = A11: B(5, 4) = A32: B(5, 6) = A12: B(5, 7) = A33
                B(5, 9) = A13: B(5, 10) = B(3, 12): B(5, 12) = B(1, 10)
     B(6, 1) = A21: B(6, 2) = A11: B(6, 4) = A22: B(6, 5) = A12: B(6, 7) = A23
                B(6, 8) = A13: B(6, 10) = B(2, 11): B(6, 11) = B(1, 10)
     '--- DB Matrix DB = D*B
     For I = 1 To 6
     For J = 1 To 12
     DB(I, J) = 0
     For K = 1 To 6
     DB(I, J) = DB(I, J) + D(I, K) * B(K, J)
     Next K: Next J: Next I
     If ISTR = 1 Then Exit Sub
     'Stress Evaluation (Element Nodal Displacements stored in QT() )
     For I = 1 To 4
        NII = 3 * (NOC(N, I) - 1): II = 3 * (I - 1)
        For J = 1 To 3
           QT(II + J) = F(NII + J)
        Next J
     Next I
     C1 = AL * DT(N)
     For I = 1 To 6
        STR(I) = 0
        For K = 1 To 12
           STR(I) = STR(I) + DB(I, K) * QT(K)
        Next K
        STR(I) = STR(I) - C1 * (D(I, 1) + D(I, 2) + D(I, 3))
     Next I
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

