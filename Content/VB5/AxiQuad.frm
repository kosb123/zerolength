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
'********     AXISYMMETRIC QUAD     **********
'*     2-D STRESS ANALYSIS USING 4-NODE      *
'*  QUADRILATERAL ELEMENTS WITH TEMPERATURE  *
'*     T.R.Chandrupatla and A.D.Belegundu    *
'*********************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), MAT(), PM()
Dim DT(), NU(), U(), F(), MPC(), BT()
Dim D(), B(), DB(), SE(), Q(), STR(), TL(), S()
Dim XNI(), A(), G()
Dim React(), vonMisesStress() As Double
Dim CNST, IPL, NQ, DJ, PI
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
     NPR = 3  'Material Properties E, Nu, Alpha
     msg = " 1) No Plot Data" & Chr(13)
     msg = msg + " 2) Create Data File for Von Mises Stress" & Chr(13)
     msg = msg + "    Choose 1 or 2"
     IPL = InputBox(msg, "Plot Choice", 1)
     Open File1 For Input As #1
     PI = 3.14159
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     '----- Total dof is  NQ
     NQ = NDN * NN
     ReDim X(NN, NDIM), NOC(NE, NEN), MAT(NE), PM(NM, NPR)
     ReDim DT(NE), NU(ND), U(ND), F(NQ)
     ReDim D(4, 4), B(4, 8), DB(4, 8), SE(8, 8), Q(8), STR(4)
     ReDim TL(8), XNI(4, 2), A(3, 4), G(4, 8), MPC(NMPC, 2), BT(NMPC, 3)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        For J = 1 To NDIM
           Input #1, X(N, J)
        Next J
     Next I
     '----- Connectivity, Material, Temp-change
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
     Call IntegPoints
     For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
        Call DMatrix(N)
        Call ElemStiffness(N)
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
     ReDim vonMisesStress(NE, 4)
     '-----  Stress Calculations
     For N = 1 To NE
        Call DMatrix(N)
        For IP = 1 To 4
           Call DbMat(N, 2, IP, RAD) '--- Get DB Matrix with Stress calculation
        '--- Von Mises Stress at Integration Point
           '--- C1 and C2 are invariants of the Stress Tensor
           C1 = STR(1) + STR(2) + STR(4)
           C2 = STR(1) * STR(2) + STR(2) * STR(4) + STR(4) * STR(1)
           C2 = C2 - STR(3) * STR(3)
           SV = Sqr(C1 * C1 - 3 * C2)
           vonMisesStress(N, IP) = SV
        Next IP
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
     cdlDialog.Filter =  "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program AxiQuad - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE#   R-Displ     Z-Displ"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(F(2 * I - 1), "0.0000E+00  ");
        Print #2, Format(F(2 * I), "0.0000E+00")
     Next I
     '----- Stress Output -----
     Print #2, "ELEM#    vonMises Stresses at 4 Integration points"
     For N = 1 To NE
        Print #2, Format(N, "@@@@@  ");
        For IP = 1 To 4
           Print #2, Format(vonMisesStress(N, IP), "0.0000E+00  ");
        Next IP
        Print #2,
     Next N
     picBox.Print "Complete results are in file "; File2
  If IPL > 1 Then
     cdlDialog.Filter = "All Files (*.*)|*.*|<Elem values> FileName.e|*.e"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File3 = cdlDialog.FileName
     Open File3 For Output As #3
     Print #3, "Von Mises Stress "
     For N = 1 To NE
        For IP = 1 To 4
           Print #3, vonMisesStress(N, IP);
        Next IP
        Print #3,
     Next N
     picBox.Print "Element Stess Data in "; File3
     picBox.Print "Run BESTFITQ and then CONTOURA or CONTOURB to plot stresses"
     Close #3
  End If
     '----- Reactions -----
     Print #2, "DOF#  REACTION"
     For I = 1 To ND
        N = NU(I)
        Print #2, Format(N, "@@@@@  "); Format(React(I), "0.0000E+00")
     Next I
     Close #2
     picBox.Print "Complete results are in file "; File2
End Sub
Private Sub IntegPoints()
'------- Integration Points XNI() --------
     C = 0.57735026919
     XNI(1, 1) = -C: XNI(1, 2) = -C
     XNI(2, 1) = C: XNI(2, 2) = -C
     XNI(3, 1) = C: XNI(3, 2) = C
     XNI(4, 1) = -C: XNI(4, 2) = C
End Sub
Private Sub DMatrix(N)
'-----  D() Matrix and Element Nodal Coordinates -----
     '--- Material Properties
     MATN = MAT(N): E = PM(MATN, 1)
     PNU = PM(MATN, 2): AL = PM(MATN, 3)
     C1 = E * (1 - PNU) / ((1 + PNU) * (1 - 2 * PNU))
     C2 = PNU / (1 - PNU)
     '-----  D() Matrix  -----
     For I = 1 To 4: For J = 1 To 4: D(I, J) = 0: Next J: Next I
     D(1, 1) = C1: D(1, 2) = C1 * C2: D(1, 4) = C1 * C2
     D(2, 1) = D(1, 2): D(2, 2) = C1: D(2, 4) = C1 * C2
     D(3, 3) = 0.5 * E / (1 + PNU)
     D(4, 1) = D(1, 4): D(4, 2) = D(2, 4): D(4, 4) = C1
End Sub
Private Sub ElemStiffness(N)
'--------  Element Stiffness and Temperature Load  -----
     For I = 1 To 8: For J = 1 To 8: SE(I, J) = 0: Next J: TL(I) = 0: Next I
     DTE = DT(N)
     '--- Weight Factor is ONE
     '--- Loop on Integration Points
     For IP = 1 To 4
        '---  Get DB Matrix at Integration Point IP
        Call DbMat(N, 1, IP, RAD)
        '--- Element Stiffness Matrix  SE
        C1 = 2 * PI * RAD * DJ
        For I = 1 To 8
           For J = 1 To 8
              C = 0
              For K = 1 To 4
                 C = C + C1 * B(K, I) * DB(K, J)
              Next K
              SE(I, J) = SE(I, J) + C
           Next J
        Next I
        '--- Determine Temperature Load TL
        C1 = AL * DTE * C1
        For I = 1 To 8
           TL(I) = C1 * (DB(1, I) + DB(2, I) + DB(4, I))
        Next I
     Next IP
End Sub
Private Sub DbMat(N, ISTR, IP, RAD)
'-------  DB()  MATRIX  ------
     XI = XNI(IP, 1): ETA = XNI(IP, 2)
     '--- Nodal Coordinates
     N1 = NOC(N, 1): N2 = NOC(N, 2)
     N3 = NOC(N, 3): N4 = NOC(N, 4)
     X1 = X(N1, 1): Y1 = X(N1, 2)
     X2 = X(N2, 1): Y2 = X(N2, 2)
     X3 = X(N3, 1): Y3 = X(N3, 2)
     X4 = X(N4, 1): Y4 = X(N4, 2)
     '--- Formation of Jacobian  TJ
     TJ11 = ((1 - ETA) * (X2 - X1) + (1 + ETA) * (X3 - X4)) / 4
     TJ12 = ((1 - ETA) * (Y2 - Y1) + (1 + ETA) * (Y3 - Y4)) / 4
     TJ21 = ((1 - XI) * (X4 - X1) + (1 + XI) * (X3 - X2)) / 4
     TJ22 = ((1 - XI) * (Y4 - Y1) + (1 + XI) * (Y3 - Y2)) / 4
     '--- Determinant of the JACOBIAN
     DJ = TJ11 * TJ22 - TJ12 * TJ21
     '--- A(3,4) Matrix relates 3 Strains eR, eZ, eRZ to
     '--- Local Derivatives of u
     For I = 1 To 3: For J = 1 To 4: A(I, J) = 0: Next J: Next I
     A(1, 1) = TJ22 / DJ: A(3, 1) = -TJ21 / DJ
     A(1, 2) = -TJ12 / DJ:  A(3, 2) = TJ11 / DJ
     A(2, 3) = -TJ21 / DJ: A(3, 3) = TJ22 / DJ
     A(2, 4) = TJ11 / DJ: A(3, 4) = -TJ12 / DJ
     '--- G(4,8) Matrix relates Local Derivatives of u
     '--- to Local Nodal Displacements q(8)
     For I = 1 To 4: For J = 1 To 8
     G(I, J) = 0: Next J: Next I
     G(1, 1) = -(1 - ETA) / 4: G(2, 1) = -(1 - XI) / 4
     G(3, 2) = -(1 - ETA) / 4: G(4, 2) = -(1 - XI) / 4
     G(1, 3) = (1 - ETA) / 4: G(2, 3) = -(1 + XI) / 4
     G(3, 4) = (1 - ETA) / 4: G(4, 4) = -(1 + XI) / 4
     G(1, 5) = (1 + ETA) / 4: G(2, 5) = (1 + XI) / 4
     G(3, 6) = (1 + ETA) / 4: G(4, 6) = (1 + XI) / 4
     G(1, 7) = -(1 + ETA) / 4: G(2, 7) = (1 - XI) / 4
     G(3, 8) = -(1 + ETA) / 4: G(4, 8) = (1 - XI) / 4
     '--- Shape Function Values
     SH1 = 0.25 * (1 - XI) * (1 - ETA)
     SH2 = 0.25 * (1 + XI) * (1 - ETA)
     SH3 = 0.25 * (1 + XI) * (1 + ETA)
     SH4 = 0.25 * (1 - XI) * (1 + ETA)
     RAD = SH1 * X1 + SH2 * X2 + SH3 * X3 + SH4 * X4
     '--- B(4,8) Matrix Relates Strains to q
     For I = 1 To 3
        For J = 1 To 8
           C = 0
           For K = 1 To 4
              C = C + A(I, K) * G(K, J)
           Next K
           B(I, J) = C
        Next J
     Next I
     B(4, 1) = SH1 / RAD: B(4, 2) = 0
     B(4, 3) = SH2 / RAD: B(4, 4) = 0
     B(4, 5) = SH3 / RAD: B(4, 6) = 0
     B(4, 7) = SH4 / RAD: B(4, 8) = 0
     '--- DB(4,8) Matrix relates Stresses to q(8)
     For I = 1 To 4
        For J = 1 To 8
           C = 0
           For K = 1 To 4
              C = C + D(I, K) * B(K, J)
           Next K:
           DB(I, J) = C
        Next J
     Next I
     If ISTR = 2 Then
           '--- Stress Evaluation
           For I = 1 To NEN
              IIN = NDN * (NOC(N, I) - 1)
              II = NDN * (I - 1)
              For J = 1 To NDN
                 Q(II + J) = F(IIN + J)
              Next J
           Next I
           C1 = AL * DT(N)
           For I = 1 To 4
              C = 0
              For K = 1 To 8
                 C = C + DB(I, K) * Q(K)
              Next K
              STR(I) = C - C1 * (D(I, 1) + D(I, 2) + D(I, 4))
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

