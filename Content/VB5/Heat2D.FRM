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
'*        PROGRAM  HEAT2D            *
'*   HEAT 2-D  WITH 3-NODED TRIANGLES  *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NN, NE, NM, NDIM, NEN, NDN
Dim ND, NL, NPR, NMPC, NBW, NHF, NCV
Dim X(), NOC(), MAT(), PM(), F(), NU(), U(), EHS()
Dim NFLUX(), NCONV(), FLUX(), H(), TINF(), Q()
Dim BT(), S()
Dim CNST, IPL
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
     Call HeatFlowCalc
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
     msg = msg + " 2) Create Data File Containing Nodal Temperatures" & Chr(13)
     msg = msg + "    Choose 1 or 2"
     IPL = InputBox(msg, "Plot Choice", 1)       '--- default is no data
     Line Input #1, Dummy: Input #1, Title
     Line Input #1, Dummy: Input #1, NN, NE, NM, NDIM, NEN, NDN
     Line Input #1, Dummy: Input #1, ND, NL, NMPC
     NPR = 1   'Material Property Thermal Conductivity
     NMPC = 0
     '--- ND = NO. OF SPECIFIED TEMPERATURES
     '--- NL = NO. OF NODAL HEAT SOURCES
     'NOTE!!  NPR =1 (THERMAL CONDUCTIVITY) AND NMPC = 0 FOR THIS PROGRAM
     '--- EHS(I) = ELEMENT HEAT SOURCE, I = 1,...,NE
     ReDim X(NN, 2), NOC(NE, 3), MAT(NE), PM(NM, NPR), F(NN)
     ReDim NU(ND), U(ND), EHS(NE)
     '=============  READ DATA  ===============
     '----- Coordinates
     Line Input #1, Dummy
     For I = 1 To NN
        Input #1, N
        For J = 1 To NDIM
           Input #1, X(N, J)
        Next J
     Next I
     '----- Connectivity, Material#, Element Heat Source
     Line Input #1, Dummy
     For I = 1 To NE
        Input #1, N
        For J = 1 To NEN
           Input #1, NOC(N, J)
        Next J
        Input #1, MAT(N), EHS(N)
     Next I
     '----- Temperature BC
     Line Input #1, Dummy
     For I = 1 To ND
        Input #1, NU(I), U(I)
     Next I
     '----- Nodal Heat Sources
     Line Input #1, Dummy
     For I = 1 To NN: F(I) = 0: Next I
     For I = 1 To NL
        Input #1, N
        Input #1, F(N)
     Next I
     '----- Thermal Conductivity of Material
     Line Input #1, Dummy
     For I = 1 To NM
        Input #1, N
        For J = 1 To NPR
           Input #1, PM(N, J)
        Next J
     Next I
     'No. of edges with Specified Heat flux FOLLOWED BY two edges & q0 (positive if out)
     Line Input #1, Dummy
     Input #1, NHF
     If NHF > 0 Then
        ReDim NFLUX(NHF, 2), FLUX(NHF)
        For I = 1 To NHF
            Input #1, NFLUX(I, 1), NFLUX(I, 2), FLUX(I)
        Next I
     End If
     'No. of Edges with Convection FOLLOWED BY edge(2 nodes) &  h & Tinf
     Line Input #1, Dummy
     Input #1, NCV
     If NCV > 0 Then
        ReDim NCONV(NCV, 2), H(NCV), TINF(NCV)
        For I = 1 To NCV
            Input #1, NCONV(I, 1), NCONV(I, 2), H(I), TINF(I)
        Next I
     End If
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
     '----- Initialization of Conductivity Matrix and Heat Rate Vector
     ReDim S(NN, NBW)
     For I = 1 To NN
        For J = 1 To NBW
           S(I, J) = 0
        Next J
     Next I
     If NHF > 0 Then
       For I = 1 To NHF
          N1 = NFLUX(I, 1): N2 = NFLUX(I, 2)
          V = FLUX(I)
          ELEN = Sqr((X(N1, 1) - X(N2, 1)) ^ 2 + (X(N1, 2) - X(N2, 2)) ^ 2)
          F(N1) = F(N1) - ELEN * V / 2
          F(N2) = F(N2) - ELEN * V / 2
       Next I
     End If
     If NCV > 0 Then
       For I = 1 To NCV
          N1 = NCONV(I, 1): N2 = NCONV(I, 2)
          ELEN = Sqr((X(N1, 1) - X(N2, 1)) ^ 2 + (X(N1, 2) - X(N2, 2)) ^ 2)
          F(N1) = F(N1) + ELEN * H(I) * TINF(I) / 2
          F(N2) = F(N2) + ELEN * H(I) * TINF(I) / 2
          S(N1, 1) = S(N1, 1) + H(I) * ELEN / 3
          S(N2, 1) = S(N2, 1) + H(I) * ELEN / 3
          If N1 >= N2 Then
             N3 = N1: N1 = N2: N2 = N3
          End If
          S(N1, N2 - N1 + 1) = S(N1, N2 - N1 + 1) + H(I) * ELEN / 6
       Next I
     End If
     '----- Conductivity Matrix
     ReDim BT(2, 3)
     For I = 1 To NE
        I1 = NOC(I, 1): I2 = NOC(I, 2): I3 = NOC(I, 3)
        X32 = X(I3, 1) - X(I2, 1): X13 = X(I1, 1) - X(I3, 1)
        X21 = X(I2, 1) - X(I1, 1)
        Y23 = X(I2, 2) - X(I3, 2): Y31 = X(I3, 2) - X(I1, 2)
        Y12 = X(I1, 2) - X(I2, 2)
        DETJ = X13 * Y23 - X32 * Y31
        AREA = 0.5 * Abs(DETJ)
        '--- Element Heat Sources
        If EHS(I) <> 0 Then
           C = EHS(I) * AREA / 3
           F(I1) = F(I1) + C: F(I2) = F(I2) + C: F(I3) = F(I3) + C
        End If
        BT(1, 1) = Y23: BT(1, 2) = Y31: BT(1, 3) = Y12
        BT(2, 1) = X32: BT(2, 2) = X13: BT(2, 3) = X21
        For II = 1 To 3
           For JJ = 1 To 2
              BT(JJ, II) = BT(JJ, II) / DETJ
           Next JJ
        Next II
        For II = 1 To 3
           For JJ = 1 To 3
              II1 = NOC(I, II): II2 = NOC(I, JJ)
              If II1 <= II2 Then
                 Sum = 0
                 For J = 1 To 2
                    Sum = Sum + BT(J, II) * BT(J, JJ)
                 Next J
                 IC = II2 - II1 + 1
                 S(II1, IC) = S(II1, IC) + Sum * AREA * PM(MAT(I), 1)
              End If
           Next JJ
        Next II
     Next I
End Sub
Private Sub ModifyForBC()
     If ND = 0 Then Exit Sub
       '--- Modify for Temp. Boundary Conditions
       Sum = 0
       For I = 1 To NN
          Sum = Sum + S(I, 1)
       Next I
       Sum = Sum / NN
       CNST = Sum * 1000000
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
Private Sub HeatFlowCalc()
     ReDim Q(NE, 2)
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
        QX = BT(1, 1) * F(I1) + BT(1, 2) * F(I2) + BT(1, 3) * F(I3)
        QX = -QX * PM(MAT(I), 1)
        QY = BT(2, 1) * F(I1) + BT(2, 2) * F(I2) + BT(2, 3) * F(I3)
        QY = -QY * PM(MAT(I), 1)
        Q(I, 1) = QX
        Q(I, 2) = QY
     Next I
End Sub
Private Sub Output()
     '===== Print Temperature, Heat Flow Data
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Heat2D - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Output for Data from " + File1
     Print #2, Title
     '----- Temperatures -----
     Print #2, "NODE#   Temperature"
     For I = 1 To NN
        Print #2, Format(I, "@@@@@  "); Format(F(I), "0.0000E+00  ")
     Next I
     '----- Stress Output -----
     Print #2, "===== Conduction Heat Flow per Unit Area in Each Element ====="
     Print #2, "ELEM#   Qx= -K*DT/Dx    Qy= -K*DT/Dy"
     For N = 1 To NE
        Print #2, Format(N, "@@@@@  ");
        Print #2, Format(Q(N, 1), "0.0000E+00  ");
        Print #2, Format(Q(N, 2), "0.0000E+00  ")
     Next N
     Close #2
     picBox.Print "Complete results are in file "; File2
  If IPL > 1 Then
     cdlDialog.Filter = "All Files (*.*)|*.*|<Element values> FileName.e|*.e"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File3 = cdlDialog.FileName
     Open File3 For Output As #3
       Print #3, "Nodal Temperatures"
     For N = 1 To NN
        Print #3, F(N)
     Next N
     picBox.Print "Nodal Temperature Data in "; File3
     picBox.Print "Run CONTOURA or CONTOURB to plot Isotherms"
     Close #3
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

