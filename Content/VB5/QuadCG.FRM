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
'********     PROGRAM QuadCG        **********
'*     2-D STRESS ANALYSIS USING 4-NODE      *
'*  QUADRILATERAL ELEMENTS WITH TEMPERATURE  *
'*         CONJUGATE GRADIENT APPROACH       *
'*    T.R.Chandrupatla and A.D.Belegundu     *
'*********************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), MAT(), PM()
Dim TH(), DT(), NU(), U(), F(), DD(), AD(), GG()
Dim D(), B(), DB(), S(), Q(), QE(), STR()
Dim TL(), XNI(), A(), G(), MPC(), BT()
Dim React(), vonMisesStress(), maxShearStress() As Double
Dim CNST, LC, IPL, NQ, DJ
Dim File1 As String, File2 As String, File3 As String
Dim Title As String, Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Stiffness
     Call ModifyForceForBC
     Call CGSolve
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
     ReDim TH(NE), DT(NE), NU(ND), U(ND), F(NQ), DD(NQ), AD(NQ), GG(NQ)
     ReDim D(3, 3), B(3, 8), DB(3, 8), S(NE, 8, 8), Q(NQ), QE(8), STR(3)
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
Private Sub Stiffness()
     '----- Stiffness Matrix of Each Element -----
     Call IntegPoints
     For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
        Call DMatrix(N)
        Call ElemStiffness(N)
        picBox.Print ".... Placing in Global Locations"
        '--- Add temperature loads
        For II = 1 To NEN
           NRT = NDN * (NOC(N, II) - 1)
           For IT = 1 To NDN
              NR = NRT + IT
              I = NDN * (II - 1) + IT
              F(NR) = F(NR) + TL(I)
              GG(NR) = GG(NR) + S(N, I, I)
           Next IT
        Next II
     Next N
End Sub
Private Sub ModifyForceForBC()
'----- GG() diagonal stiffness summation
'----- Decide Penalty Parameter CNST -----
     CNST = 0
     For I = 1 To NQ
        If CNST < GG(I) Then CNST = GG(I)
     Next I
     CNST = CNST * 10000
'----- Modify right hand side F() for Boundary Conditions -----
     '--- Displacement BC ---
     For I = 1 To ND
        N = NU(I)
        F(N) = F(N) + CNST * U(I)
     Next I
     '--- Multi-point Constraints ---
     For I = 1 To NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        F(I1) = F(I1) + CNST * BT(I, 1) * BT(I, 3)
        F(I2) = F(I2) + CNST * BT(I, 2) * BT(I, 3)
     Next I
End Sub
Private Sub CGSolve()
     '----- Conjugate Gradient Method starts here
     GG1 = 0
     For I = 1 To NQ
        GG(I) = -F(I): DD(I) = F(I): Q(I) = 0
        GG1 = GG1 + GG(I) * GG(I)
     Next I
     '-----  ITERATION LOOP
     ITER = 0
    Do
     ITER = ITER + 1
     Print "ITERATION# = "; ITER
     '=====  ELEMENT LOOP  =====
     For N = 1 To NE
        For I = 1 To 4
           IGT = 2 * (NOC(N, I) - 1): ILT = 2 * (I - 1)
           For II = 1 To 2
              IG = IGT + II: IL = ILT + II
              For J = 1 To 4
                 JGT = 2 * (NOC(N, J) - 1): JLT = 2 * (J - 1)
                 For JJ = 1 To 2
                    JG = JGT + JJ: JL = JLT + JJ
                    AD(IG) = AD(IG) + S(N, IL, JL) * DD(JG)
                 Next JJ
              Next J
           Next II
        Next I
     Next N
     '--- Displacement BC ---
     For I = 1 To ND
        N = NU(I)
        AD(N) = AD(N) + CNST * DD(N)
     Next I
     '--- Multi-point Constraints ---
     For I = 1 To NMPC
        I1 = MPC(I, 1): I2 = MPC(I, 2)
        C = BT(I, 1) * DD(I1) + BT(I, 2) * DD(I2)
        AD(I1) = AD(I1) + CNST * BT(I, 1) * C
        AD(I2) = AD(I2) + CNST * BT(I, 2) * C
     Next I
        DAD = 0
        For I = 1 To NQ
           DAD = DAD + DD(I) * AD(I)
        Next I
        AL = GG1 / DAD
        GG2 = 0
        For I = 1 To NQ
           GG(I) = GG(I) + AL * AD(I)
           Q(I) = Q(I) + AL * DD(I)
           GG2 = GG2 + GG(I) * GG(I)
        Next I
        If GG2 < 0.00000001 Then Exit Do
        BTA = GG2 / GG1
        GG1 = GG2
        For I = 1 To NQ
           DD(I) = -GG(I) + BTA * DD(I)
        Next I
        For I = 1 To NQ: AD(I) = 0: Next I
    Loop
End Sub
Private Sub StressCalc()
     ReDim vonMisesStress(NE, 4), maxShearStress(NE, 4)
     '-----  Stress Calculations
     For N = 1 To NE
        Call DMatrix(N)
        For IP = 1 To 4
           Call DbMat(N, 2, IP) '--- Get DB Matrix with Stress calculation
        '--- Von Mises Stress at Integration Point
           C = 0: If LC = 2 Then C = PNU * (STR(1) + STR(2))
           C1 = (STR(1) - STR(2)) ^ 2 + (STR(2) - C) ^ 2 + (C - STR(1)) ^ 2
           SV = Sqr(0.5 * C1 + 3 * STR(3) ^ 2)
        '--- Maximum Shear Stress R
           R = Sqr(0.25 * (STR(1) - STR(2)) ^ 2 + (STR(3)) ^ 2)
           maxShearStress(N, IP) = R
           vonMisesStress(N, IP) = SV
        Next IP
     Next N
End Sub
Private Sub ReactionCalc()
     ReDim React(ND)
     '----- Reaction Calculation -----
     For I = 1 To ND
        N = NU(I)
        React(I) = CNST * (U(I) - Q(N))
     Next I
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Quad - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     If LC = 1 Then Print #2, "Plane Stress Analysis"
     If LC = 2 Then Print #2, "Plane Strain Analysis"
     '----- Displacements -----
     Print #2, "NODE#   X-Displ     Y-Displ"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(Q(2 * I - 1), "0.0000E+00  ");
        Print #2, Format(Q(2 * I), "0.0000E+00")
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
     cdlDialog.Filter = "All Files (*.*)|*.*|<Element values> FileName.e|*.e"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File3 = cdlDialog.FileName
     Open File3 For Output As #3
     If IPL = 2 Then Print #3, "Max. in-plane Shear Stress "
     If IPL = 3 Then Print #3, "Von Mises Stress "
     For N = 1 To NE
        For IP = 1 To 4
           If IPL = 2 Then Print #3, maxShearStress(N, IP);
           If IPL = 3 Then Print #3, vonMisesStress(N, IP);
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
'-----  D() Matrix  -----
     '--- Material Properties
     MATN = MAT(N)
     E = PM(MATN, 1): PNU = PM(MATN, 2)
     AL = PM(MATN, 3)
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
End Sub
Private Sub ElemStiffness(N)
'--------  Element Stiffness and Temperature Load  -----
     For I = 1 To 8: TL(I) = 0: Next I
     DTE = DT(N)
     '--- Weight Factor is ONE
     '--- Loop on Integration Points
     For IP = 1 To 4
        '---  Get DB Matrix at Integration Point IP
        Call DbMat(N, 1, IP)
        '--- Element Stiffness Matrix  SE
        For I = 1 To 8
           For J = 1 To 8
              C = 0
              For K = 1 To 3
                 C = C + B(K, I) * DB(K, J) * DJ * TH(N)
              Next K
              S(N, I, J) = S(N, I, J) + C
           Next J
        Next I
        '--- Determine Temperature Load TL
        AL = PM(MAT(N), 3)
        C = AL * DTE: If LC = 2 Then C = (1 + PNU) * C
        For I = 1 To 8
           TL(I) = TL(I) + TH(N) * DJ * C * (DB(1, I) + DB(2, I))
        Next I
     Next IP
End Sub
Private Sub DbMat(N, ISTR, IP)
'-------  DB()  MATRIX  ------
     XI = XNI(IP, 1): ETA = XNI(IP, 2)
     '--- Nodal Coordinates
     THICK = TH(N)
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
     '--- A(3,4) Matrix relates Strains to
     '--- Local Derivatives of u
     A(1, 1) = TJ22 / DJ: A(2, 1) = 0: A(3, 1) = -TJ21 / DJ
     A(1, 2) = -TJ12 / DJ: A(2, 2) = 0: A(3, 2) = TJ11 / DJ
     A(1, 3) = 0: A(2, 3) = -TJ21 / DJ: A(3, 3) = TJ22 / DJ
     A(1, 4) = 0: A(2, 4) = TJ11 / DJ: A(3, 4) = -TJ12 / DJ
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
     '--- B(3,8) Matrix Relates Strains to q
     For I = 1 To 3
        For J = 1 To 8
           C = 0
           For K = 1 To 4
              C = C + A(I, K) * G(K, J)
           Next K
           B(I, J) = C
        Next J
     Next I
     '--- DB(3,8) Matrix relates Stresses to q(8)
     For I = 1 To 3
        For J = 1 To 8
           C = 0
           For K = 1 To 3
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
                 QE(II + J) = Q(IIN + J)
              Next J
           Next I
           AL = PM(MAT(N), 3)
           C1 = AL * DT(N): If LC = 2 Then C1 = C1 * (1 + PNU)
           For I = 1 To 3
              C = 0
              For K = 1 To 8
                 C = C + DB(I, K) * QE(K)
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

