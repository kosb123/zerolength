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
'********  PROGRAM SKYLINE   ***********
'*    SKYLINE METHOD FOR EQUATIONS     *
'* T.R.Chandrupatla and A.D.Belegundu  *
'***************************************
DefInt I-N
DefDbl A-H, O-Z
Dim N, NSUM
Dim A(), B(), ID()
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call SkyLine
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     Line Input #1, Title: Line Input #1, Dummy
     Input #1, N
     Line Input #1, Dummy
     NSUM = 0
     ReDim ID(N)
     '--- Column Heights then convert to Diagonal Pointers
     Input #1, ID(1)
     For I = 2 To N
        Input #1, ID(I)
        ID(I) = ID(I) + ID(I - 1)
     Next I
     NSUM = ID(N)
     ReDim A(NSUM), B(N)
     NT = 0: NI = 1
     Line Input #1, Dummy
     For I = 1 To N
        If I > 1 Then
           NT = ID(I - 1)
           NI = ID(I) - ID(I - 1)
        End If
        For J = 1 To NI
           Input #1, A(NT + J)
        Next J
     Next I
     Line Input #1, Dummy
     For I = 1 To N
        Input #1, B(I)
     Next I
     Close #1
End Sub
Private Sub SkyLine()
  '----- Skyline Solver
     For J = 2 To N
        NJ = ID(J) - ID(J - 1)
        If NJ > 1 Then
           K1 = 0: NJ = J - NJ + 1
           For K = NJ To J - 1
              K1 = K1 + 1: KJ = ID(J - 1) + K1: KK = ID(K)
              C = A(KJ) / A(KK)
              For I = K + 1 To J
                 NI = ID(I) - ID(I - 1)
                 If (I - K + 1) <= NI Then
                    IJ = ID(J) - J + I
                    KI = ID(I) - I + K
                    A(IJ) = A(IJ) - C * A(KI)
                 End If
              Next I
           Next K
        End If
     Next J
     For K = 1 To N - 1
        KK = ID(K)
        C = B(K) / A(KK)
        For I = K + 1 To N
           NI = ID(I) - ID(I - 1)
           If (I - K + 1) <= NI Then
               KI = ID(I) - I + K
               B(I) = B(I) - C * A(KI)
           End If
        Next I
     Next K
     '----- Back-substitution -----
     NS = ID(N): B(N) = B(N) / A(NS)
     For I1 = 1 To N - 1
        I = N - I1
        II = ID(I)
        C = 1 / A(II): B(I) = C * B(I)
        For J = I + 1 To N
           J1 = J - I + 1: NJ = ID(J) - ID(J - 1)
           If J1 <= NJ Then
              IJ = ID(J) - J + I
              B(I) = B(I) - C * A(IJ) * B(J)
           End If
        Next J
     Next I1
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program Gauss - CHANDRUPATLA & BELEGUNDU"
     Print #2, Title
     '----- Solution Vector -----
     For I = 1 To N
        Print #2, "X("; I; ")= "; Format(B(I), "0.0000E+00")
     Next I
     Close #2
     picBox.Print "RESULTS ARE IN FILE "; File2
End Sub
Private Sub cmdView_Click()
   Dim ALine As String, CRLF As String, File1 As String
   CRLF = Chr$(13) + Chr$(10)
   picBox.Visible = False
   txtView.Visible = True
   Open File2 For Input As #1
   txtView.Text = ""
   Do While Not EOF(1)
     Line Input #1, ALine
     txtView.Text = txtView.Text + ALine + CRLF
   Loop
   Close #1
End Sub

