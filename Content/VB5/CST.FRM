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
'*           PROGRAM CST               *
'*     CONSTANT STRAIN TRIANGLE        *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), MAT(), PM()
Dim TH(), DT(), NU(), U(), F(), MPC(), BT()
Dim D(), B(), DB(), SE(), Q(), STR(), TL(), S()
Dim Stress(), React(), PrinStress(), PltStress()
Dim CNST, LC, IPL, NQ, DJ
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
     msg = "1. Plane Stress Analysis" & Chr(13)
     msg = msg + "2. Plane Strain Analysis" & Chr(13)
     msg = msg + "   Your Choice <1 or 2> ?"
     LC = InputBox(msg, "Type of Problem", 1)
     NPR = 3   'Material Properties E, Nu, Alpha
     msg = " 1) No Plot Data" & Chr(13)
     msg = msg + " 2) Create Data File for in-plane Shear Stress" & Chr(13)
     msg = msg + " 3) Create Data File for Von Mises Stress" & Chr(13)
     msg = msg + "    Choose 1 or 2 or 3"
     IPL = InputBox(msg, "Plot Choice", 1)
     Open File1 For Input As #1
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     '----- Total dof is  NQ
     NQ = NDN * NN
     ReDim X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     ReDim TH(NE), DT(NE), NU(ND), U(ND), F(NQ), MPC(NMPC, 2), BT(NMPC, 3)
     ReDim D(3, 3), B(3, 6), DB(3, 6), SE(6, 6), Q(6), STR(3), TL(6)
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
Private Sub Stiffness()
     ReDim S(NQ, NBW)
     '----- Global Stiffness Matrix -----
          For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
     Call DbMat(N, 1)
     '--- Element Stiffness
        For I = 1 To 6
           For J = 1 To 6
              C = 0
              For K = 1 To 3
                 C = C + 0.5 * Abs(DJ) * B(K, I) * DB(K, J) * TH(N)
              Next K
              SE(I, J) = C
           Next J
        Next I
     '--- Temperature Load Vector
        AL = PM(MAT(N), 3)
        C = AL * DT(N): If LC = 2 Then C = C * (1 + PNU)
        For I = 1 To 6
           TL(I) = 0.5 * C * TH(N) * Abs(DJ) * (DB(1, I) + DB(2, I))
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
              F(NR) = F(NR) + TL(I)
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
     ReDim Stress(NE, 3), PrinStress(NE, 3), PltStress(NE)
     '-----  Stress Calculations
     For N = 1 To NE
        Call DbMat(N, 2)
     '--- Principal Stress Calculations
        If STR(3) = 0 Then
           S1 = STR(1): S2 = STR(2): ANG = 0
           If S2 > S1 Then
              S1 = STR(2): S2 = STR(1): ANG = 90
           End If
        Else
           C = 0.5 * (STR(1) + STR(2))
           R = Sqr(0.25 * (STR(1) - STR(2)) ^ 2 + (STR(3)) ^ 2)
           S1 = C + R: S2 = C - R
           If C > STR(1) Then
              ANG = 57.2957795 * Atn(STR(3) / (S1 - STR(1)))
              If STR(3) > 0 Then ANG = 90 - ANG
              If STR(3) < 0 Then ANG = -90 - ANG
           Else
              ANG = 57.29577951 * Atn(STR(3) / (STR(1) - S2))
           End If
        End If
        Stress(N, 1) = STR(1)
        Stress(N, 2) = STR(2)
        Stress(N, 3) = STR(3)
        PrinStress(N, 1) = S1
        PrinStress(N, 2) = S2
        PrinStress(N, 3) = ANG
        If IPL = 2 Then PltStress(N) = 0.5 * (S1 - S2)
        If IPL = 3 Then
           S3 = 0: If LC = 2 Then S3 = PNU * (S1 + S2)
           C = (S1 - S2) ^ 2 + (S2 - S3) ^ 2 + (S3 - S1) ^ 2
           PltStress(N) = Sqr(0.5 * C)
        End If
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
     Print #2, "Program Cst2D - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     If LC = 1 Then Print #2, "Plane Stress Analysis"
     If LC = 2 Then Print #2, "Plane Strain Analysis"
     '----- Displacements -----
     Print #2, "NODE#   X-Displ     Y-Displ"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(F(2 * I - 1), "0.0000E+00  ");
        Print #2, Format(F(2 * I), "0.0000E+00")
     Next I
     '----- Stress Output -----
     Print #2, "ELEM#    SX          SY          TXY";
     Print #2, "         S1          S2        ANGLE SX->S1"
     For N = 1 To NE
        Print #2, Format(N, "@@@@@  ");
        Print #2, Format(Stress(N, 1), "0.0000E+00  ");
        Print #2, Format(Stress(N, 2), "0.0000E+00  ");
        Print #2, Format(Stress(N, 3), "0.0000E+00  ");
        Print #2, Format(PrinStress(N, 1), "0.0000E+00  ");
        Print #2, Format(PrinStress(N, 2), "0.0000E+00  ");
        Print #2, Format(PrinStress(N, 3), "0.0000E+00  ")
     Next N
     '----- Reactions -----
     Print #2, "DOF#  REACTION"
     For I = 1 To ND
        N = NU(I)
        Print #2, Format(N, "@@@@@  "); Format(React(I), "0.0000E+00")
     Next I
     Close #2
     picBox.Print "Complete results are in file "; File2
  If IPL > 1 Then
     cdlDialog.Filter = "All Files (*.*)|*.*|<Element values> FileName.e|*.e"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File3 = cdlDialog.FileName
     Open File3 For Output As #3
     Print #3, "Element Stress Values"
     For N = 1 To NE
        Print #3, PltStress(N)
     Next N
     picBox.Print "Element Stess Data in "; File3
     picBox.Print "Run BESTFIT and then CONTOURA or CONTOURB to plot stresses"
     Close #3
  End If
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

