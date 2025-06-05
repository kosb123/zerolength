VERSION 5.00
Begin VB.Form frmFemCB 
   Caption         =   "FemCB"
   ClientHeight    =   6165
   ClientLeft      =   60
   ClientTop       =   345
   ClientWidth     =   8715
   LinkTopic       =   "Form1"
   ScaleHeight     =   6165
   ScaleWidth      =   8715
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
      ScaleHeight     =   5235
      ScaleWidth      =   7875
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
'*****        PROGRAM HEXAFNT          *****
'*   3-D STRESS ANALYSIS USING  8-NODE     *
'*    ISOPARAMETRIC HEXAHEDRAL ELEMENT     *
'*          USING FRONTAL SOLVER           *
'*   T.R.Chandrupatla and A.D.Belegundu    *
'*******************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW
Dim X(), NOC(), NU(), U(), MAT(), F()
Dim SE(), PM(), MPC(), BT()
Dim XI(), XNI(), D(), GN(), H(), TJ()
Dim AJ(), G(), B(), DB(), QT(), STR(), DT()
Dim IDE(), ISBL(), S(), IEBL(), INDX()
Dim React(), vonMisesStress(), maxShearStress() As Double
Dim CNST, NQ, DJ, MTN1, MTN, NEDF, NFRON, NTOGO, NDCNT, ICOUNT
Dim File1 As String, File2 As String
Dim Title As String, Dummy As String
Private Type Record
   VarNum As Integer
   Coeff As Double
End Type
Dim Adat As Record
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call PreFront
     RecordLen = Len(Adat)
     Open "SCRATCH.DAT" For Random As #3 Len = RecordLen     'Scratch file for writing
     Call Stiffness
     Call BackSub
     Close #3
     Kill "SCRATCH.DAT"
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
     NPR = 3   ' Material Properties E, Nu, Alpha
     ReDim X(NN, NDIM), NOC(NE, NEN), NU(ND), U(ND), MAT(NE), F(NQ)
     ReDim SE(24, 24), PM(NM, NPR), MPC(NMPC, 2), BT(NMPC, 3)
     ReDim XI(3, 8), XNI(3, 8), D(6, 6), GN(3, 8), H(9, 24), TJ(3, 3)
     ReDim AJ(3, 3), G(6, 9), B(6, 24), DB(6, 24), QT(24), STR(6), DT(NE)
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
Private Sub PreFront()
        '----- Mark Last Appearance of Node / Make it negative in NOC()
        ' Last appearance is first appearance for reverse element order
        NEDF = NEN * NDN
        For I = 1 To NN
           II = 0
           For J = NE To 1 Step -1
              For K = 1 To NEN
                 If I = NOC(J, K) Then
                     II = 1
                     Exit For
                 End If
              Next K
              If II = 1 Then Exit For
           Next J
           NOC(J, K) = -I
        Next I
        '===== Block Size Determination
        NQ = NN * NDN
        ReDim IDE(NQ)
        For I = 1 To NQ: IDE(I) = 0: Next I
        For I = 1 To NMPC: For J = 1 To 2: IDE(MPC(I, J)) = 1: Next J: Next I
        IFRON = 0: For I = 1 To NQ: IFRON = IFRON + IDE(I): Next I
        IBL = IFRON
        For N = 1 To NE
           INEG = 0
           For I = 1 To NEN
              I1 = NOC(N, I): IA = NDN * (Abs(I1) - 1)
              For J = 1 To NDN
                 IA = IA + 1
                 If IDE(IA) = 0 Then
                    IFRON = IFRON + 1: IDE(IA) = 1
                 End If
              Next J
              If I1 < 0 Then INEG = INEG + 1
           Next I
           If IBL < IFRON Then IBL = IFRON
           IFRON = IFRON - NDN * INEG
        Next N
        Erase IDE
        ReDim ISBL(IBL), S(IBL, IBL), IEBL(NEDF), INDX(IBL)
        NFRON = 0: NTOGO = 0: NDCNT = 0
        For I = 1 To IBL: INDX(I) = I: Next I
End Sub
Private Sub Stiffness()
     '----- Global Stiffness Matrix -----
     Call IntegPoints
     MTN1 = 0
     For N = 1 To NE
        picBox.Print "Forming Stiffness Matrix of Element "; N
        MTN = MAT(N)
        If MTN <> MTN1 Then
           Call DMatrix(N)
        End If
        Call ElemStiffness(N)
        If N = 1 Then
           CNST = 0
           For I = 1 To NEDF: CNST = CNST + SE(I, I): Next I
           CNST = 100000000000# * CNST
           Call MpcFron
        End If
        '----- Account for temperature loads QT()
        For I = 1 To NEN
          IL = 3 * (I - 1): IG = 3 * (Abs(NOC(N, I)) - 1)
          For J = 1 To 3
             IL = IL + 1: IG = IG + 1
             F(IG) = F(IG) + QT(IL)
          Next J
        Next I
        Call Front(N)         'Frontal assembly  and Forward Elimination
     Next N
End Sub
Private Sub StressCalc()
     ReDim vonMisesStress(NE, 8)
     '-----  Stress Calculations
     MTN1 = 0
     For N = 1 To NE
        MTN = MAT(N)
        If MTN <> MTN1 Then
           Call DMatrix(N)
        End If
        For IP = 1 To 8
        '--- Von Mises Stress at Integration Points
           Call DbMat(N, 2, IP) '--- Get DB Matrix with Stress calculation
           '--- Calculation of Von Mises Stress at IP
           SIV1 = STR(1) + STR(2) + STR(3)
           SIV2 = STR(1) * STR(2) + STR(2) * STR(3) + STR(3) * STR(1)
           SIV2 = SIV2 - STR(4) ^ 2 - STR(5) ^ 2 - STR(6) ^ 2
           vonMisesStress(N, IP) = Sqr(SIV1 * SIV1 - 3 * SIV2)
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
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program HexaFront - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Displacements -----
     Print #2, "NODE#   X-Displ     Y-Displ     Z-Displ"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(F(3 * I - 2), "0.0000E+00  ");
        Print #2, Format(F(3 * I - 1), "0.0000E+00  "); Format(F(3 * I), "0.0000E+00")
     Next I
     '----- Stress Output -----
     For N = 1 To NE
     Print #2, "vonMises Stresses at 8 Integration points in ELEM#  ";
        Print #2, N
        For IP = 1 To 8
           Print #2, Format(vonMisesStress(N, IP), "0.0000E+00  ");
           If IP Mod 4 = 0 Then Print #2,
        Next IP
     Next N
     picBox.Print "Complete results are in file "; File2
     '----- Reactions -----
     Print #2, "NODE#  REACTION"
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
     XI(1, 1) = -1: XI(2, 1) = -1: XI(3, 1) = -1
     XI(1, 2) = 1: XI(2, 2) = -1: XI(3, 2) = -1
     XI(1, 3) = 1: XI(2, 3) = 1: XI(3, 3) = -1
     XI(1, 4) = -1: XI(2, 4) = 1: XI(3, 4) = -1
     XI(1, 5) = -1: XI(2, 5) = -1: XI(3, 5) = 1
     XI(1, 6) = 1: XI(2, 6) = -1: XI(3, 6) = 1
     XI(1, 7) = 1: XI(2, 7) = 1: XI(3, 7) = 1
     XI(1, 8) = -1: XI(2, 8) = 1: XI(3, 8) = 1
     For I = 1 To 8
        XNI(1, I) = C * XI(1, I)
        XNI(2, I) = C * XI(2, I)
        XNI(3, I) = C * XI(3, I)
     Next I
End Sub
Private Sub DMatrix(N)
     '--- D() Matrix relating Stresses to Strains
     E = PM(MTN, 1): PNU = PM(MTN, 2): AL = PM(MTN, 3)
     C1 = E / ((1 + PNU) * (1 - 2 * PNU))
     C2 = 0.5 * E / (1 + PNU)
     For I = 1 To 6: For J = 1 To 6: D(I, J) = 0: Next J: Next I
     D(1, 1) = C1 * (1 - PNU): D(1, 2) = C1 * PNU: D(1, 3) = D(1, 2)
     D(2, 1) = D(1, 2): D(2, 2) = D(1, 1): D(2, 3) = D(1, 2)
     D(3, 1) = D(1, 3): D(3, 2) = D(2, 3): D(3, 3) = D(1, 1)
     D(4, 4) = C2: D(5, 5) = C2: D(6, 6) = C2
     MTN1 = MTN
End Sub
Private Sub ElemStiffness(N)
'--------  Element Stiffness  -----
     For I = 1 To 24: For J = 1 To 24
     SE(I, J) = 0: Next J: QT(I) = 0: Next I
     DTE = DT(N)
     '--- Weight Factor is ONE
     '--- Loop on Integration Points
     For IP = 1 To 8
        '---  Get DB Matrix at Integration Point IP
        Call DbMat(N, 1, IP)
        '--- Element Stiffness Matrix  SE
        For I = 1 To 24
           For J = 1 To 24
              For K = 1 To 6
                 SE(I, J) = SE(I, J) + B(K, I) * DB(K, J) * DJ
              Next K
           Next J
        Next I
        '--- Determine Temperature Load QT()
        C = AL * DTE
        For I = 1 To 24
           DSUM = DB(1, I) + DB(2, I) + DB(3, I)
           QT(I) = QT(I) + C * Abs(DJ) * DSUM / 6
        Next I
     Next IP
End Sub
Private Sub DbMat(N, ISTR, IP)
'-------  DB()  MATRIX  ------
     '--- Gradient of Shape Functions - The GN() Matrix
     For I = 1 To 3
        For J = 1 To 8
           C = 1
           For K = 1 To 3
              If K <> I Then
                 C = C * (1 + XI(K, J) * XNI(K, IP))
              End If
           Next K
           GN(I, J) = 0.125 * XI(I, J) * C
        Next J
     Next I
     '--- Formation of Jacobian  TJ
     For I = 1 To 3
        For J = 1 To 3
           TJ(I, J) = 0
           For K = 1 To 8
              KN = Abs(NOC(N, K))
              TJ(I, J) = TJ(I, J) + GN(I, K) * X(KN, J)
           Next K
        Next J
     Next I
     '--- Determinant of the JACOBIAN
     DJ1 = TJ(1, 1) * (TJ(2, 2) * TJ(3, 3) - TJ(3, 2) * TJ(2, 3))
     DJ2 = TJ(1, 2) * (TJ(2, 3) * TJ(3, 1) - TJ(3, 3) * TJ(2, 1))
     DJ3 = TJ(1, 3) * (TJ(2, 1) * TJ(3, 2) - TJ(3, 1) * TJ(2, 2))
     DJ = DJ1 + DJ2 + DJ3
     '--- Inverse of the Jacobian AJ()
     AJ(1, 1) = (TJ(2, 2) * TJ(3, 3) - TJ(2, 3) * TJ(3, 2)) / DJ
     AJ(1, 2) = (TJ(3, 2) * TJ(1, 3) - TJ(3, 3) * TJ(1, 2)) / DJ
     AJ(1, 3) = (TJ(1, 2) * TJ(2, 3) - TJ(1, 3) * TJ(2, 2)) / DJ
     AJ(2, 1) = (TJ(2, 3) * TJ(3, 1) - TJ(2, 1) * TJ(3, 3)) / DJ
     AJ(2, 2) = (TJ(1, 1) * TJ(3, 3) - TJ(1, 3) * TJ(3, 1)) / DJ
     AJ(2, 3) = (TJ(1, 3) * TJ(2, 1) - TJ(1, 1) * TJ(2, 3)) / DJ
     AJ(3, 1) = (TJ(2, 1) * TJ(3, 2) - TJ(2, 2) * TJ(3, 1)) / DJ
     AJ(3, 2) = (TJ(1, 2) * TJ(3, 1) - TJ(1, 1) * TJ(3, 2)) / DJ
     AJ(3, 3) = (TJ(1, 1) * TJ(2, 2) - TJ(1, 2) * TJ(2, 1)) / DJ
     '--- H() Matrix relates local derivatives of  u  to local
     '    displacements  q
     For I = 1 To 9
        For J = 1 To 24
           H(I, J) = 0
        Next J
     Next I
     For I = 1 To 3
        For J = 1 To 3
           IR = 3 * (I - 1) + J
           For K = 1 To 8
              IC = 3 * (K - 1) + I
              H(IR, IC) = GN(J, K)
           Next K
        Next J
     Next I
     '--- G() Matrix relates strains to local derivatives of  u
     For I = 1 To 6
        For J = 1 To 9
           G(I, J) = 0
        Next J
     Next I
     G(1, 1) = AJ(1, 1): G(1, 2) = AJ(1, 2): G(1, 3) = AJ(1, 3)
     G(2, 4) = AJ(2, 1): G(2, 5) = AJ(2, 2): G(2, 6) = AJ(2, 3)
     G(3, 7) = AJ(3, 1): G(3, 8) = AJ(3, 2): G(3, 9) = AJ(3, 3)
     G(4, 4) = AJ(3, 1): G(4, 5) = AJ(3, 2): G(4, 6) = AJ(3, 3)
          G(4, 7) = AJ(2, 1): G(4, 8) = AJ(2, 2): G(4, 9) = AJ(2, 3)
     G(5, 1) = AJ(3, 1): G(5, 2) = AJ(3, 2): G(5, 3) = AJ(3, 3)
          G(5, 7) = AJ(1, 1): G(5, 8) = AJ(1, 2): G(5, 9) = AJ(1, 3)
     G(6, 1) = AJ(2, 1): G(6, 2) = AJ(2, 2): G(6, 3) = AJ(2, 3)
          G(6, 4) = AJ(1, 1): G(6, 5) = AJ(1, 2): G(6, 6) = AJ(1, 3)
     '--- B() Matrix relates strains to  q
     For I = 1 To 6
        For J = 1 To 24
           B(I, J) = 0
           For K = 1 To 9
              B(I, J) = B(I, J) + G(I, K) * H(K, J)
           Next K
        Next J
     Next I
     '--- DB() Matrix relates stresses to  q
     For I = 1 To 6
        For J = 1 To 24
           DB(I, J) = 0
           For K = 1 To 6
              DB(I, J) = DB(I, J) + D(I, K) * B(K, J)
           Next K
        Next J
     Next
     If ISTR = 1 Then Exit Sub
           '--- Element Nodal Displacements stored in QT()
           For I = 1 To 8
              IIN = 3 * (Abs(NOC(N, I)) - 1)
              II = 3 * (I - 1)
              For J = 1 To 3
                 QT(II + J) = F(IIN + J)
              Next J
           Next I
           '--- Stress Calculation STR = DB * Q
           For I = 1 To 6
              STR(I) = 0
              For J = 1 To 24
                 STR(I) = STR(I) + DB(I, J) * QT(J)
              Next J
              STR(I) = STR(I) - CAL * (D(I, 1) + D(I, 2) + D(I, 3))
           Next I
End Sub
Private Sub MpcFron()
        '----- Modifications for Multipoint Constraints by Penalty Method
        For I = 1 To NMPC
           I1 = MPC(I, 1)
           IFL = 0
           For J = 1 To NFRON
              J1 = INDX(J)
              If I1 = ISBL(J1) Then
                 IFL = 1: Exit For
              End If
           Next J
           If IFL = 0 Then
              NFRON = NFRON + 1: J1 = INDX(NFRON): ISBL(J1) = I1
           End If
           I2 = MPC(I, 2)
           IFL = 0
           For K = 1 To NFRON
              K1 = INDX(K)
              If K1 = ISBL(K1) Then
                 IFL = 1: Exit For
              End If
           Next K
           If IFL = 0 Then
              NFRON = NFRON + 1: K1 = INDX(NFRON): ISBL(K1) = I2
           End If
           '----- Stiffness Modification
           S(J1, J1) = S(J1, J1) + CNST * BT(I, 1) ^ 2
           S(K1, K1) = S(K1, K1) + CNST * BT(I, 2) ^ 2
           S(J1, K1) = S(J1, K1) + CNST * BT(I, 1) * BT(I, 2)
           S(K1, J1) = S(J1, K1)
           '----- Force Modification
           F(I1) = F(I1) + CNST * BT(I, 3) * BT(I, 1)
           F(I2) = F(I2) + CNST * BT(I, 3) * BT(I, 2)
        Next I
End Sub
Private Sub Front(N)
'----- Frontal Method Assembly and Elimination -----
'----------------  Assembly of Element N  --------------------
        For I = 1 To NEN
           I1 = NOC(N, I): IA = Abs(I1): IS1 = Sgn(I1)
           IDF = NDN * (IA - 1): IE1 = NDN * (I - 1)
           For J = 1 To NDN
              IDF = IDF + 1: IE1 = IE1 + 1: IFL = 0
              If NFRON > NTOGO Then
                 For II = NTOGO + 1 To NFRON
                    IX = INDX(II)
                    If IDF = ISBL(IX) Then
                       IFL = 1: Exit For
                    End If
                 Next II
              End If
              If IFL = 0 Then
                 NFRON = NFRON + 1: II = NFRON: IX = INDX(II)
              End If
              ISBL(IX) = IDF: IEBL(IE1) = IX
              If IS1 = -1 Then
                 NTOGO = NTOGO + 1
                 ITEMP = INDX(NTOGO)
                 INDX(NTOGO) = INDX(II)
                 INDX(II) = ITEMP
              End If
           Next J
        Next I
        For I = 1 To NEDF
           I1 = IEBL(I)
           For J = 1 To NEDF
              J1 = IEBL(J)
              S(I1, J1) = S(I1, J1) + SE(I, J)
           Next J
        Next I
'------------------------------------------------------------------
     If NDCNT < ND Then
'-----  Modification for displacement BCs / Penalty Approach  -----
        For I = 1 To NTOGO
           I1 = INDX(I)
           IG = ISBL(I1)
              For J = 1 To ND
                 If IG = NU(J) Then
                    S(I1, I1) = S(I1, I1) + CNST
                    F(IG) = F(IG) + CNST * U(J)
                    NDCNT = NDCNT + 1       'Counter for check
                    Exit For
                 End If
              Next J
        Next I
     End If
'------------   Elimination of completed variables   ---------------
        NTG1 = NTOGO
        For II = 1 To NTG1
           IPV = INDX(1): IPG = ISBL(IPV)
           Pivot = S(IPV, IPV)
        '-----  Write separator "0" and PIVOT value to disk  -----
           Adat.VarNum = 0
           Adat.Coeff = Pivot
           ICOUNT = ICOUNT + 1
           Put #3, ICOUNT, Adat
           S(IPV, IPV) = 0
           For I = 2 To NFRON
              I1 = INDX(I): IG = ISBL(I1)
              If S(I1, IPV) <> 0 Then
                  C = S(I1, IPV) / Pivot: S(I1, IPV) = 0
                  For J = 2 To NFRON
                     J1 = INDX(J)
                     If S(IPV, J1) <> 0 Then
                        S(I1, J1) = S(I1, J1) - C * S(IPV, J1)
                     End If
                  Next J
                  F(IG) = F(IG) - C * F(IPG)
              End If
           Next I
           For J = 2 To NFRON
        '-----  Write Variable# and Reduced Coeff/PIVOT to disk  -----
              J1 = INDX(J)
              If S(IPV, J1) <> 0 Then
                 ICOUNT = ICOUNT + 1: IBA = ISBL(J1)
                 Adat.VarNum = IBA
                 Adat.Coeff = S(IPV, J1) / Pivot
                 Put #3, ICOUNT, Adat
                 S(IPV, J1) = 0
              End If
           Next J
           ICOUNT = ICOUNT + 1
        '-----  Write Eliminated Variable# and RHS/PIVOT to disk  -----
           Adat.VarNum = IPG
           Adat.Coeff = F(IPG) / Pivot
           F(IPG) = 0
           Put #3, ICOUNT, Adat
        '----- (NTOGO) into (1); (NFRON) into (NTOGO)
        '----- IPV into (NFRON) and reduce front & NTOGO sizes by 1
           If NTOGO > 1 Then
              INDX(1) = INDX(NTOGO)
           End If
           INDX(NTOGO) = INDX(NFRON): INDX(NFRON) = IPV
           NFRON = NFRON - 1: NTOGO = NTOGO - 1
        Next II
End Sub
Private Sub BackSub()
        '===== Backsubstitution
        Do While ICOUNT > 0
           Get #3, ICOUNT, Adat
           ICOUNT = ICOUNT - 1
           N1 = Adat.VarNum
           F(N1) = Adat.Coeff
           Do
              Get #3, ICOUNT, Adat
              ICOUNT = ICOUNT - 1
              N2 = Adat.VarNum
              If N2 = 0 Then Exit Do
              F(N1) = F(N1) - Adat.Coeff * F(N2)
           Loop
        Loop
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

