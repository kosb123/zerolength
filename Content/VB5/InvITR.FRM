VERSION 5.00
Begin VB.Form frmFemCB 
   Caption         =   "FemCB"
   ClientHeight    =   6228
   ClientLeft      =   60
   ClientTop       =   348
   ClientWidth     =   8616
   LinkTopic       =   "Form1"
   ScaleHeight     =   6228
   ScaleWidth      =   8616
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
'*****       PROGRAM INVITR        *****
'*      Inverse Iteration Method       *
'*  for Eigenvalues and Eigenvectors   *
'*        Searching in Subspace        *
'*         for Banded Matrices         *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim NQ, NBW
Dim S(), GM(), EV1(), EV2(), EVC(), EVL()
Dim EVT(), EVS(), ST(), NITER()
Dim TOL, SH, NEV, NEV1, ITMAX, PI
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call BanSolve1        '<----Stiffness to Upper Triangle
     Call InverseIter
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     TOL = InputBox("Enter Value", "Tolerance", 0.000001)
     NEV = InputBox("Enter Number", "Number of Eigenvalues Desired", 1)
     SH = 0
     Open File1 For Input As #1
     Line Input #1, Title
     Line Input #1, Dummy
     Input #1, NQ, NBW
     ReDim S(NQ, NBW), GM(NQ, NBW), EV1(NQ), EV2(NQ), NITER(NEV)
     ReDim EVT(NQ), EVS(NQ), ST(NQ), EVC(NQ, NEV), EVL(NEV)
     '=============  READ DATA  ===============
     '----- Read in Stiffness Matrix -----
     Line Input #1, Dummy
     For I = 1 To NQ
        For J = 1 To NBW
           Input #1, S(I, J)
        Next J
     Next I
     '----- Read in Mass Matrix
     Line Input #1, Dummy
     For I = 1 To NQ
        For J = 1 To NBW
           Input #1, GM(I, J)
        Next J
     Next I
     '----- Starting Vector for Inverse Iteration
     Line Input #1, Dummy
     For I = 1 To NQ
           Input #1, ST(I)
     Next I
     Close #1
     SH = InputBox("SHIFT", "Shift Value for Eigenvalue", 0)
     If SH <> 0 Then
         For I = 1 To NQ: For J = 1 To NBW
             S(I, J) = S(I, J) - SH * GM(I, J)
         Next J: Next I
     End If
End Sub
Private Sub InverseIter()
     ITMAX = 50: NEV1 = NEV
     PI = 3.14159
     For NV = 1 To NEV
        '--- Starting Value for Eigenvector
        For I = 1 To NQ
           EV1(I) = ST(I)
        Next I
        EL2 = 0: ITER = 0
        Do
           EL1 = EL2
           ITER = ITER + 1
           If ITER > ITMAX Then
              picBox.Print "No Convergence for Eigenvalue# "; NV
              NEV1 = NV - 1
              Exit Sub
           End If
           If NV > 1 Then
              '----  Starting Vector Orthogonal to
              '----       Evaluated Vectors
              For I = 1 To NV - 1
                 CV = 0
                 For K = 1 To NQ
                    KA = K - NBW + 1: KZ = K + NBW - 1
                    If KA < 1 Then KA = 1
                    If KZ > NQ Then KZ = NQ
                    For L = KA To KZ
                       If L < K Then
                          K1 = L: L1 = K - L + 1
                       Else
                          K1 = K: L1 = L - K + 1
                       End If
                       CV = CV + EVS(K) * GM(K1, L1) * EVC(L, I)
                    Next L
                 Next K
                 For K = 1 To NQ
                    EV1(K) = EV1(K) - CV * EVC(K, I)
                 Next K
              Next I
           End If
           For I = 1 To NQ
              IA = I - NBW + 1: IZ = I + NBW - 1: EVT(I) = 0
              If IA < 1 Then IA = 1
              If IZ > NQ Then IZ = NQ
              For K = IA To IZ
                 If K < I Then
                    I1 = K: K1 = I - K + 1
                 Else
                    I1 = I: K1 = K - I + 1
                 End If
                 EVT(I) = EVT(I) + GM(I1, K1) * EV1(K)
              Next K
              EV2(I) = EVT(I)
           Next I
           Call BanSolve2          '<--- Reduce Right Side and Solve
           C1 = 0: C2 = 0
           For I = 1 To NQ
              C1 = C1 + EV2(I) * EVT(I)
           Next I
           For I = 1 To NQ
              IA = I - NBW + 1: IZ = I + NBW - 1: EVT(I) = 0
              If IA < 1 Then IA = 1
              If IZ > NQ Then IZ = NQ
              For K = IA To IZ
                 If K < I Then
                    I1 = K: K1 = I - K + 1
                 Else
                    I1 = I: K1 = K - I + 1
                 End If
                 EVT(I) = EVT(I) + GM(I1, K1) * EV2(K)
              Next K
           Next I
           For I = 1 To NQ
              C2 = C2 + EV2(I) * EVT(I)
           Next I
           EL2 = C1 / C2
           C2 = Sqr(C2)
           For I = 1 To NQ
              EV1(I) = EV2(I) / C2
              EVS(I) = EV1(I)
           Next I
        Loop While Abs(EL2 - EL1) / Abs(EL2) > TOL
        For I = 1 To NQ
           EVC(I, NV) = EV1(I)
        Next I
        NITER(NV) = ITER
        EL2 = EL2 + SH
        EVL(NV) = EL2
     Next NV
End Sub
Private Sub BanSolve1()
'-----  Gauss Elimination LDU Approach (for Symmetric Banded Matrices)
     '----- Reduction to Upper Triangular Form
     For K = 1 To NQ - 1
        NK = NQ - K + 1
        If NK > NBW Then NK = NBW
        For I = 2 To NK
           C1 = S(K, I) / S(K, 1)
           I1 = K + I - 1
           For J = I To NK
              J1 = J - I + 1
              S(I1, J1) = S(I1, J1) - C1 * S(K, J)
           Next J
        Next I
     Next K
End Sub
Private Sub BanSolve2()
     '----- Reduction of the right hand side
     For K = 1 To NQ - 1
        NK = NQ - K + 1
        If NK > NBW Then NK = NBW
        For I = 2 To NK: I1 = K + I - 1
           C1 = 1 / S(K, 1)
           EV2(I1) = EV2(I1) - C1 * S(K, I) * EV2(K)
        Next I
     Next K
     '----- Back Substitution
     EV2(NQ) = EV2(NQ) / S(NQ, 1)
     For II = 1 To NQ - 1
        I = NQ - II: C1 = 1 / S(I, 1)
        NI = NQ - I + 1
        If NI > NBW Then NI = NBW
        EV2(I) = C1 * EV2(I)
        For K = 2 To NI
           EV2(I) = EV2(I) - C1 * S(I, K) * EV2(I + K - 1)
        Next K
     Next II
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program InvItr - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Eigenvalues & Eigenvectors for data in file: "; File1
     '----- Eigenvalues and Eigenvectors -----
     If NEV1 < NEV Then
        Print #2, "Convergence for "; NEV1; " Eigenvalues Only."
        NEV = NEV1
     End If
     For NV = 1 To NEV
        Print #2,
        Print #2, "Eigenvalue Number "; NV;
        Print #2, "     Iteration Number "; NITER(NV)
        Print #2, "Eigenvalue = ";
        Print #2, Format(EVL(NV), "0.0000E+00  ");
        OMEGA = Sqr(EVL(NV)): FREQ = 0.5 * OMEGA / PI
        Print #2, "Omega = ";
        Print #2, Format(OMEGA, "0.0000E+00  ");
        Print #2, "Freq Hz = ";
        Print #2, Format(FREQ, "0.0000E+00")
        Print #2, "Eigenvector "
        For I = 1 To NQ
           Print #2, Format(EVC(I, NV), "0.0000E+00  ");
        Next I
        Print #2,
     Next NV
     Close #2
     picBox.Print "RESULTS ARE IN FILE "; File2
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

