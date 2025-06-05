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
'********************************************
'*      GENERALIZED EIGENVALUE PROBLEM      *
'*              Kx = lambda Mx              *
'*    T.R.Chandrupatla and A.D.Belegundu    *
'********************************************
DefInt I-N
DefDbl A-H, O-Z
Dim S(), GM(), D(), B(), NORD()
Dim NQ, ITER, PI
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData(S(), GM(), D(), B(), NORD(), NQ)
'----- Cholesky Factorization of Mass Matrix
     Call Cholesky(GM(), NQ)
'----- Update of Stiffness Matrix - Standard Form Ax=(lambda)x
     Call UpdateStiff(S(), GM(), NQ)
'----- Tri-diagonalize D() Diagonal, B() Sub-diagonal
'----- S() has the rotation matrix
     Call TriDiag(S(), D(), B(), NQ)
'----- Find Eigenvalues and Eigenvectors
     Call EigenTD(D(), B(), S(), NQ, ITER)
'----- Determine Eigenvectors
     Call EigenVEC(S(), GM(), NQ)
'----- Ascending Order of Eigenvalues
     Call AscendEigen(D(), NORD(), NQ)
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData(S(), GM(), D(), B(), NORD(), NQ)
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, Title
     Line Input #1, Dummy
     Input #1, NQ, NBW
     PI = 3.14159
     ReDim S(NQ, NQ), GM(NQ, NQ)
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
Private Sub Cholesky(A(), N)
     '----- L into lower left triangle of A
     For K = 1 To N
        A(K, K) = Sqr(A(K, K))
        For I = K + 1 To N
           A(I, K) = A(I, K) / A(K, K)
        Next I
        For J = K + 1 To N
           For I = J To N
              A(I, J) = A(I, J) - A(I, K) * A(J, K)
           Next I
        Next J
     Next K
End Sub
Sub UpdateStiff(A(), B(), N)
     '--- Forward Substitution I  (invL)*A
     For J = 1 To N
        A(1, J) = A(1, J) / B(1, 1)
        For I = 2 To N
           For K = 1 To I - 1
              A(I, J) = A(I, J) - B(I, K) * A(K, J)
           Next K
           A(I, J) = A(I, J) / B(I, I)
        Next I
     Next J
     '--- Forward Substitution II   (invL)*A*(invL)'
     For J = 1 To N
        A(J, 1) = A(J, 1) / B(1, 1)
        For I = 2 To N
           For K = 1 To I - 1
              A(J, I) = A(J, I) - B(I, K) * A(J, K)
           Next K
           A(J, I) = A(J, I) / B(I, I)
        Next I
     Next J
End Sub
Sub TriDiag(A(), D(), B(), N)
    ReDim D(N), B(N - 1)
    For I = 1 To N - 2
       AA = 0
       For J = I + 1 To N
          AA = AA + A(J, I) * A(J, I)
       Next J
       AA = Sqr(AA)
       WW = 2 * AA * (AA + Abs(A(I + 1, I)))
       WW = Sqr(WW)
       IA = Sgn(A(I + 1, I))
       '----- Diagonal and Next to Diagonal Term
       D(I) = A(I, I)
       B(I) = -IA * AA
       '----- Unit Vector W() in Column I from Row I+1 to N
       For J = I + 1 To N
          A(J, I) = A(J, I) / WW
       Next J
       A(I + 1, I) = A(I + 1, I) + IA * AA / WW
       '----- W'A in Row I from Col I+1 to N
       BET = 0
       For J = I + 1 To N
          A(I, J) = 0
          For K = I + 1 To N
              A(I, J) = A(I, J) + A(K, I) * A(K, J)
          Next K
          BET = BET + A(I, J) * A(J, I)
       Next J
       '----- Modified A()
       For J = I + 1 To N
          For K = I + 1 To N
             A(J, K) = A(J, K) - 2 * A(J, I) * A(I, K) - 2 * A(I, J) * A(K, I)
             A(J, K) = A(J, K) + 4 * BET * A(K, I) * A(J, I)
          Next K
       Next J
    Next I
    D(N - 1) = A(N - 1, N - 1)
    B(N - 1) = A(N - 1, N)
    D(N) = A(N, N)
    A(N - 1, N - 1) = 1
    A(N, N) = 1
    A(N - 1, N) = 0
    A(N, N - 1) = 0
    '----- Now Create the Q matrix in A()
    For I = 1 To N - 2
       II = N - I - 1
       A(II, II) = 1
       For J = 1 To I + 1
          IJ = II + J
          A(II, IJ) = 0
          For K = 1 To I + 1
             IK = II + K
             A(II, IJ) = A(II, IJ) + A(IK, IJ) * A(IK, II)
          Next K
          For K = 1 To I + 1
             IK = II + K
             A(IK, IJ) = A(IK, IJ) - 2 * A(IK, II) * A(II, IJ)
          Next K
       Next J
       For J = 1 To I + 1
          IJ = II + J
          A(II, IJ) = 0
          A(IJ, II) = 0
       Next J
    Next I
End Sub
Sub EigenTD(A(), B(), Q(), N, ITER)
     ITER = 0
     M = N
  Do
     ITER = ITER + 1
     dd = 0.5 * (A(M - 1) - A(M))
     BB = B(M - 1) * B(M - 1)
     BOT = dd + Sgn(dd) * Sqr(dd * dd + BB)
     P = A(1) - A(M) + BB / BOT
     X = B(1)
     For I = 1 To M - 1
        PP = Sqr(P * P + X * X)
        cs = -P / PP
        sn = X / PP
        If I > 1 Then B(I - 1) = cs * B(I - 1) - sn * X
        A1 = A(I): A2 = A(I + 1): B1 = B(I)
        A(I) = A1 * cs * cs - 2 * B1 * cs * sn + A2 * sn * sn
        B(I) = (A1 - A2) * cs * sn + B1 * (cs * cs - sn * sn)
        A(I + 1) = A1 * sn * sn + 2 * B1 * cs * sn + A2 * cs * cs
        '----- Update Q()
        For K = 1 To N
           A1 = Q(K, I): A2 = Q(K, I + 1)
           Q(K, I) = cs * A1 - sn * A2
           Q(K, I + 1) = sn * A1 + cs * A2
        Next K
        If I = M - 1 Then Exit For
        X = -B(I + 1) * sn
        B(I + 1) = B(I + 1) * cs
        P = B(I)
     Next I
     Do While M > 1 And Abs(B(M - 1)) < 0.000001
           M = M - 1
     Loop
  Loop While M > 1
  M = M - 1
End Sub
Sub EigenVEC(A(), B(), N)
     '--- Backsubstitution ---
     For J = 1 To N
        A(N, J) = A(N, J) / B(N, N)
        For I = N - 1 To 1 Step -1
           For K = N To I + 1 Step -1
              A(I, J) = A(I, J) - B(K, I) * A(K, J)
           Next K
           A(I, J) = A(I, J) / B(I, I)
        Next I
     Next J
End Sub
Sub AscendEigen(D(), NORD(), NQ)
     ReDim NORD(NQ)
     For I = 1 To NQ: NORD(I) = I: Next I
     '----- Ascending Order of Eigenvalues
     For I = 1 To NQ - 1
        II = NORD(I): I1 = II
        C1 = D(II): J1 = I
        For J = I + 1 To NQ
           IJ = NORD(J)
           If C1 > D(IJ) Then
              C1 = D(IJ): I1 = IJ: J1 = J
           End If
        Next J
        NORD(I) = I1: NORD(J1) = II
     Next I
End Sub
Private Sub Output()
     '===== Print Eigenvalues and Eigenvectors etc
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Geneigen - CHANDRUPATLA & BELEGUNDU"
     Print #2, "Eigenvalues & Eigenvectors for data in file: "; File1
     '----- Eigenvalues and Eigenvectors -----
     For I = 1 To NQ
        II = NORD(I)
        Print #2,
        Print #2, "Eigenvalue Number "; I
        C = D(II): OMEGA = Sqr(C): FREQ = 0.5 * OMEGA / PI
        Print #2, "Eigenvalue = ";
        Print #2, Format(C, "0.0000E+00  ");
        Print #2, "Omega = ";
        Print #2, Format(OMEGA, "0.0000E+00  ");
        Print #2, "Freq Hz = ";
        Print #2, Format(FREQ, "0.0000E+00")
        Print #2, "Eigenvector "
        For J = 1 To NQ
           Print #2, Format(S(J, II), "0.0000E+00  ");
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

