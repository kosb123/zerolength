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
'***************************************
'*           PROGRAM JACOBI            *
'*    Generalized Jacobi's Method      *
'*       for symmetric matrices        *
'*  T.R.Chandrupatla and A.D.Belegundu *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim S(), GM(), EVL(), EVC(), NORD()
Dim TOL, PI, NQ, CAB, NSW, NSWMAX
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call Jacobi
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, Title
     Line Input #1, Dummy
     Input #1, NQ, NBW
     ReDim S(NQ, NQ), GM(NQ, NQ), EVL(NQ), EVC(NQ, NQ), NORD(NQ)
     '=============  READ DATA  ===============
     '----- Banded Stiffness Matrix into S(NQ, NQ) -----
     Line Input #1, Dummy
     For I = 1 To NQ
        For JN = 1 To NBW
           Input #1, STIFF
           J = I + JN - 1
           If J <= NQ Then
              S(I, J) = STIFF: S(J, I) = STIFF
           End If
        Next JN
     Next I
     '----- Banded Mass Matrix into GM(NQ,NQ) -----
     Line Input #1, Dummy
     For I = 1 To NQ
        For JN = 1 To NBW
           Input #1, EMASS
           J = I + JN - 1
           If J <= NQ Then
              GM(I, J) = EMASS: GM(J, I) = EMASS
           End If
        Next JN
     Next I
     Close #1
End Sub
Private Sub Jacobi()
     TOL = InputBox("Enter Value", "Tolerance", 0.000001)
     PI = 3.14159
     '----- NORD( ) is for ascending order of Eigenvalues
     For I = 1 To NQ: NORD(I) = I: Next I
     '----- Initialize Eigenvector Matrix -----
     For I = 1 To NQ
        For J = 1 To NQ
           EVC(I, J) = 0
        Next J
        EVC(I, I) = 1
     Next I
     C1 = S(1, 1): C2 = GM(1, 1)
     For I = 1 To NQ
        If C1 > S(I, I) Then C1 = S(I, I)
        If C2 < GM(I, I) Then C2 = GM(I, I)
     Next I
     TOLS = TOL * C1: TOLM = TOL * C2: NSWMAX = 50
     NSWMAX = InputBox("Enter Number", "Maximum Number of Sweeps", 50)
     K1 = 1: I1 = 1: NSW = 0
   Do
     NSW = NSW + 1
     If NSW > NSWMAX Then
        '--- No convergence
        Exit Sub
     End If
     picBox.Print "Sweep Number  "; NSW
     For K = K1 To NQ - 1
        For I = I1 To K
           J = NQ - K + I
           If Abs(S(I, J)) > TOLS Or Abs(GM(I, J)) > TOLM Then
              AA = S(I, I) * GM(I, J) - GM(I, I) * S(I, J)
              BB = S(J, J) * GM(I, J) - GM(J, J) * S(I, J)
              CC = S(I, I) * GM(J, J) - GM(I, I) * S(J, J)
              CAB = 0.25 * CC * CC + AA * BB
              If CAB < 0 Then
                 picBox.Print "Square Root of Negative Term -- Check Matrices"
                 Exit Sub
              End If
              If AA = 0 Then
                 BET = 0: ALP = -S(I, J) / S(I, I)
              ElseIf BB = 0 Then
                 ALP = 0: BET = -S(I, J) / S(J, J)
              Else
                 SQC = Sqr(CAB): If CC < 0 Then SQC = -SQC
                 ALP = (-0.5 * CC + SQC) / AA
                 BET = -AA * ALP / BB
              End If
              '----- Only Upper Triangular Part is used in Diagonalization
              If I > 1 Then
                 For N = 1 To I - 1
                    SI = S(N, I): SJ = S(N, J)
                    EMI = GM(N, I): EMJ = GM(N, J)
                    S(N, I) = SI + BET * SJ
                    S(N, J) = SJ + ALP * SI
                    GM(N, I) = EMI + BET * EMJ
                    GM(N, J) = EMJ + ALP * EMI
                 Next N
              End If
              If J < NQ Then
                 For N = J + 1 To NQ
                    SI = S(I, N): SJ = S(J, N)
                    EMI = GM(I, N): EMJ = GM(J, N)
                    S(I, N) = SI + BET * SJ
                    S(J, N) = SJ + ALP * SI
                    GM(I, N) = EMI + BET * EMJ
                    GM(J, N) = EMJ + ALP * EMI
                 Next N
              End If
              If I < J Then
                 For N = I + 1 To J - 1
                    SI = S(I, N): SJ = S(N, J)
                    EMI = GM(I, N): EMJ = GM(N, J)
                    S(I, N) = SI + BET * SJ
                    S(N, J) = SJ + ALP * SI
                    GM(I, N) = EMI + BET * EMJ
                    GM(N, J) = EMJ + ALP * EMI
                 Next N
              End If
              SII = S(I, I)
              SIJ = S(I, J)
              SJJ = S(J, J)
              S(I, J) = 0
              S(I, I) = SII + 2 * BET * SIJ + BET * BET * SJJ
              S(J, J) = SJJ + 2 * ALP * SIJ + ALP * ALP * SII
              EII = GM(I, I)
              EIJ = GM(I, J)
              EJJ = GM(J, J)
              GM(I, J) = 0
              GM(I, I) = EII + 2 * BET * EIJ + BET * BET * EJJ
              GM(J, J) = EJJ + 2 * ALP * EIJ + ALP * ALP * EII
              '----- EIGENVECTORS -----
              For N = 1 To NQ
                 EVI = EVC(N, I): EVJ = EVC(N, J)
                 EVC(N, I) = EVI + BET * EVJ
                 EVC(N, J) = EVJ + ALP * EVI
              Next N
           End If
        Next I
     Next K
     For K = 1 To NQ - 1
        For I = 1 To K
           J = NQ - K + I
           IFL = 0
           If Abs(S(I, J)) > TOLS Or Abs(GM(I, J)) > TOLM Then
              K1 = K: I1 = I: IFL = 1
           End If
          If IFL = 1 Then Exit For
        Next I
        If IFL = 1 Then Exit For
     Next K
   Loop While IFL = 1
     '-----  Calculation of Eigenvalues -----
     For I = 1 To NQ
        If Abs(GM(I, I)) < TOLM Then GM(I, I) = TOLM
           EVL(I) = S(I, I) / GM(I, I)
     Next I
     '----- Scaling of Eigenvectors
     For I = 1 To NQ
        GM2 = Sqr(Abs(GM(I, I)))
        For J = 1 To NQ
           EVC(J, I) = EVC(J, I) / GM2
        Next J
     Next I
     '-----   RESULTS   -----
     '--- Ascending Order of Eigenvalues
     For I = 1 To NQ
        II = NORD(I): I1 = II
        C1 = EVL(II): J1 = I
        For J = I To NQ
           IJ = NORD(J)
           If C1 > EVL(IJ) Then
              C1 = EVL(IJ): I1 = IJ: J1 = J
           End If
        Next J
        If I1 <> II Then
           NORD(I) = I1: NORD(J1) = II
        End If
     Next I
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Jacobi - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Eigenvalues & Eigenvectors for data in file: "; File1
     '----- Eigenvalues and Eigenvectors -----
     If NSW > NSWMAX Then
        Print #2, "No convergence for "; NSWMAX; " sweeps"
        Close #2
        Exit Sub
     End If
     If CAB < 0 Then
         Print #2, "Square Root of Negative Term -- Check Matrices"
         Close #2
         Exit Sub
     End If
     For I = 1 To NQ
        II = NORD(I)
        Print #2,
        Print #2, "Eigenvalue Number "; I
        C = EVL(II): OMEGA = Sqr(C): FREQ = 0.5 * OMEGA / PI
        Print #2, "Eigenvalue = ";
        Print #2, Format(C, "0.0000E+00  ");
        Print #2, "Omega = ";
        Print #2, Format(OMEGA, "0.0000E+00  ");
        Print #2, "Freq Hz = ";
        Print #2, Format(FREQ, "0.0000E+00")
        Print #2, "Eigenvector "
        For J = 1 To NQ
           Print #2, Format(EVC(J, II), "0.0000E+00  ");
        Next J
        Print #2,
     Next I
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

