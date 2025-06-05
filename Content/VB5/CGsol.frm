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
'"*********   PROGRAM CGSOL   *********
'*     CONJUGATE GRADIENT METHOD      *
'*   FOR SOLVING AX=B, A Symmetric    *
'* T.R.Chandrupatla and A.D.Belegundu *
'**************************************
DefInt I-N
DefDbl A-H, O-Z
Dim N, ITER
Dim A(), B(), X(), G(), D(), AD()
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call CgSol
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
     ReDim A(N, N), B(N), X(N)
     '=============  READ DATA  ===============
     For I = 1 To N
        For J = 1 To N
           Input #1, A(I, J)
        Next J
     Next I
     Line Input #1, Dummy
     For I = 1 To N
        Input #1, B(I)
     Next I
     Close #1
End Sub
Private Sub CgSol()
        ReDim G(N), D(N), AD(N)
        For I = 1 To N
           X(I) = 0
           G(I) = -B(I)
           D(I) = B(I)
        Next I
        GG1 = 0
        For I = 1 To N
           GG1 = GG1 + G(I) * G(I)
        Next I
Do While GG1 > 0.000001
        ITER = ITER + 1
        DAD = 0
        For I = 1 To N
           C = 0
           For J = 1 To N
              C = C + A(I, J) * D(J)
           Next J
           AD(I) = C
           DAD = DAD + C * D(I)
        Next I
        AL = GG1 / DAD
        GG2 = 0
        For I = 1 To N
           X(I) = X(I) + AL * D(I)
           G(I) = G(I) + AL * AD(I)
           GG2 = GG2 + G(I) * G(I)
        Next I
        BT = GG2 / GG1
        For I = 1 To N
           D(I) = -G(I) + BT * D(I)
        Next I
        GG1 = GG2
    Loop
        Erase G, D, AD
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     File2 = InputBox("Output File d:\dir\fileName.ext", "Name of File")
     Open File2 For Output As #2
     Print #2, "Program Gauss - CHANDRUPATLA & BELEGUNDU"
     Print #2, Title
     '----- Solution Vector -----
     For I = 1 To N
        Print #2, "X("; I; ")= "; Format(X(I), "0.0000E+00")
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

