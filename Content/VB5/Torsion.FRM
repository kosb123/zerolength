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
'*          PROGRAM  TORSION           *
'*   TORSION WITH 3-NODED TRIANGLES    *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW, NHF, NCV
Dim X(), NOC(), MAT(), PM(), F(), NU(), U()
Dim NFLUX(), NCONV(), FLUX(), H(), TINF(), TAU()
Dim BT(), S()
Dim CNST, IPL, ALPHA, TORQUE, SFAC
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
     Call AngleOfTwist
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
     msg = " 1) No Plot Data" & Chr(13)
     msg = msg + " 2) Create Data File Containing Stress Function Values" & Chr(13)
     msg = msg + "    Choose 1 or 2"
     IPL = InputBox(msg, "Plot Choice", 1)       '--- default is no data
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     NM = 1
     '--- ND = NO. OF SPECIFIED STRESS FUNCTION VALUES
     '--- NL = NO. OF GENERALIZED NODAL FORCES "0" HERE
     NPR = 1 ' Material Property(SHEAR MODULUS)
     NMPC = 0
     '--- ELEMENT CHARACTERISTIC NOT USED
     '--- NO. OF MATERIALS = 1 FOR THIS PROGRAM
     ReDim X(NN, 2), NOC(NE, 3), MAT(NE), PM(NM, NPR), F(NN)
     ReDim NU(ND), U(ND)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        For J = 1 To NDIM
           Input #1, X(N, J)
        Next J
     Next I
     '----- Connectivity, Material#
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
        For J = 1 To NEN
           Input #1, NOC(N, J)
        Next J
        Input #1, MAT(N)
     Next I
     '----- Boundary Conditions
     Line Input #1, Dummy
     For I = 1 To ND
        Input #1, NU(I), U(I)
     Next I
     '----- Dummy Read
     Line Input #1, Dummy
     '----- Shear Modulus of Material
     Line Input #1, Dummy
     For I = 1 To NM
        Input #1, N
        For J = 1 To NPR
           Input #1, PM(N, J)
        Next J
     Next I
     '----- Applied Torque
     Line Input #1, Dummy
     Input #1, TORQUE
     '----- SYMMETRY FACTOR (eg. if 1/4 symmetry, then =4.0)
     Line Input #1, Dummy
     Input #1, SFAC
     Close #1
End Sub
Private Sub Bandwidth()
     '----- Bandwidth Evaluation -----
     IDIFF = 0
     For K = 1 To NE
        For I = 1 To 2
           For J = I + 1 To 3
              II = Abs(NOC(K, J) - NOC(K, I))
              If II > IDIFF Then IDIFF = II
           Next J
        Next I
     Next K
     NBW = IDIFF + 1
     picBox.Print "The Bandwidth is"; NBW
End Sub
Private Sub Stiffness()
     '----- Initialization of Stiffness Matrix
     ReDim S(NN, NBW)
     For I = 1 To NN
        For J = 1 To NBW
           S(I, J) = 0
        Next J
     Next I
     '----- Stiffness Matrix
     ReDim BT(2, 3)
     For I = 1 To NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        AREA = 0.5 * Abs(DETJ)
        '--- Load Calculation
           C = 2 * AREA / 3
           F(I1) = F(I1) + C: F(I2) = F(I2) + C: F(I3) = F(I3) + C
        '--- Stiffness Formation
        BT(1, 1) = Y23: BT(1, 2) = Y31: BT(1, 3) = Y12
        BT(2, 1) = X32: BT(2, 2) = X13: BT(2, 3) = X21
        For II = 1 To 3
           For JJ = 1 To 3
              II1 = NOC(I, II): II2 = NOC(I, JJ)
              If II1 <= II2 Then
                 Sum = 0
                 For J = 1 To 2
                    Sum = Sum + BT(J, II) * BT(J, JJ)
                 Next J
                 S(II1, II2 - II1 + 1) = S(II1, II2 - II1 + 1) + Sum * AREA
              End If
           Next JJ
        Next II
     Next I
End Sub
Private Sub ModifyForBC()
     If ND = 0 Then Exit Sub
       '--- Modify for Temp. Boundary Conditions
     '--- MODIFY FOR BOUNDARY CONDITIONS
     CNST = S(1, 1)
     For I = 2 To NN
       If CNST < S(I, 1) Then CNST = S(I, 1)
     Next I
     CNST = CNST * 1000000
     For I = 1 To ND
       N = NU(I)
       S(N, 1) = S(N, 1) + CNST
       F(N) = F(N) + CNST * U(I)
     Next I
End Sub
Private Sub BandSolver()
     '----- Band Solver -----
     N1 = NN - 1
     '--- Forward Elimination
     For K = 1 To N1
        NK = NN - K + 1
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
     F(NN) = F(NN) / S(NN, 1)
     For KK = 1 To N1
        K = NN - KK
        C1 = 1 / S(K, 1)
        F(K) = C1 * F(K)
        NK = NN - K + 1
        If NK > NBW Then NK = NBW
        For J = 2 To NK
           F(K) = F(K) - C1 * S(K, J) * F(K + J - 1)
        Next J
     Next KK
End Sub
Private Sub AngleOfTwist()
     ReDim TAU(NE, 2)
     Dim msg1 As String, msg2 As String
     '---- Angle of Twist per Unit Length
     Sum = 0
     For I = 1 To NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        Sum = Sum + Abs(DETJ) / 3 * (F(I1) + F(I2) + F(I3))
     Next I
     SMOD = PM(1, 1)
     ALPHA = TORQUE / SMOD / Sum / SFAC
     '----- Shear Stress Calculations
     For I = 1 To NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        BT(1, 1) = Y23: BT(1, 2) = Y31: BT(1, 3) = Y12
        BT(2, 1) = X32: BT(2, 2) = X13: BT(2, 3) = X21
        For II = 1 To 3
           For JJ = 1 To 2
              BT(JJ, II) = BT(JJ, II) / DETJ
           Next JJ
        Next II
        TAUYZ = -(BT(1, 1) * F(I1) + BT(1, 2) * F(I2) + BT(1, 3) * F(I3))
        TAUXZ = BT(2, 1) * F(I1) + BT(2, 2) * F(I2) + BT(2, 3) * F(I3)
        TAUYZ = TAUYZ * SMOD * ALPHA
        TAUXZ = TAUXZ * SMOD * ALPHA
        TAU(I, 1) = TAUXZ
        TAU(I, 2) = TAUXZ
     Next I
End Sub
Private Sub Output()
     '===== Print Temperature, Heat Flow Data
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Torsion - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Stress Function Values -----
     Print #2, "NODE#   Stress Function Value"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(F(I), "0.0000E+00  ")
     Next I
     '----- Torque and Shear Stresses -----
     Print #2, "TWIST PER UNIT LENGTH = "; ALPHA
     Print #2, "===== Shear Stresses TauYZ, TauXZ in Each Element ====="
     Print #2, "ELEM#    TauYZ           TauXZ"
     For N = 1 To NE
        Print #2, Format(N, "@@@@@  ");
        Print #2, Format(TAU(N, 1), "0.0000E+00  ");
        Print #2, Format(TAU(N, 2), "0.0000E+00  ")
     Next N
     Close #2
     picBox.Print "Complete results are in file "; File2
  If IPL > 1 Then
     cdlDialog.Filter = "All Files (*.*)|*.*|<Element values> FileName.e|*.e"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File3 = cdlDialog.FileName
     Open File3 For Output As #3
       Print #3, "Stress Function Values"
     For N = 1 To NN
        Print #3, F(N)
     Next N
     Close #3
     picBox.Print "Stress Function Values in Data File "; File3
     picBox.Print "Run CONTOURA or CONTOURB to plot Stress Function Contours"
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

