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
'************************************************
'*               PROGRAM MESHGEN                *
'*  MESH GENERATOR FOR TWO DIMENSIONAL REGIONS  *
'*     (c) T.R.CHANDRUPATLA & A.D.BELEGUNDU     *
'************************************************
DefInt I-N
DefSng A-H, O-Z
Dim NS, NW, NSJ, NSR, NWR, NNS, NNW, NNT, NGN, NODE
Dim NN, NE, NM, NEN
Dim IDBLK(), NSD(), NWD(), NGCN(), SR(), WR(), SH()
Dim X(), XB(), XP(), NOC(), MAT(), MERG(), NNAR()
Dim Title As String, File1 As String, File2 As String
Dim Dummy As String
Private Sub cmdEnd_Click()
   End
End Sub
Private Sub cmdStart_Click()
     Call InputData
     Call GlobalNode
     Call CoordConnect
     Call Output
     cmdView.Enabled = True
     cmdStart.Enabled = False
End Sub
Private Sub InputData()
     cdlDialog.Filter = "All Files (*.*)|*.*|<Input file> FileName.inp|*.inp"
     cdlDialog.ShowOpen
     File1 = cdlDialog.FileName
     Open File1 For Input As #1
     '=============  READ DATA  ===============
     Line Input #1, Dummy: Line Input #1, Title
     Line Input #1, Dummy
     Input #1, NEN          ' NEN = 3 for Triangle  4 for Quad
     If NEN < 3 Then NEN = 3
     If NEN > 4 Then NEN = 4
     'Hints:  A region is divided into 4-cornered blocks viewed as a
     '        mapping from a Checkerboard pattern of S- and W- Sides
     '        S- Side is one with lower number of final divisions
     '        Blocks, Corners, S- and W- Sides are labeled as shown in Fig. 12.2
     '        Make a sketch and identify void blocks and merging sides
     '----- Block Data -----
     '#S-Spans(NS)  #W-Spans(NW)  #PairsOfEdgesMerged(NSJ)
     Line Input #1, Dummy: Line Input #1, Dummy
     Input #1, NS, NW, NSJ
     NSW = NS * NW: NGN = (NS + 1) * (NW + 1): NM = 1
     ReDim IDBLK(NSW), NSD(NS), NWD(NW), NGCN(NGN), SH(8)
     '------------- Span Divisions ---------------
     Line Input #1, Dummy
     NNS = 1: NNW = 1
     '--- Number of divisions for each S-Span
     Line Input #1, Dummy
     For KS = 1 To NS
        Input #1, N
        Input #1, NSD(N)
        NNS = NNS + NSD(N)
     Next KS
     '--- Number of divisions for each W-Span
     Line Input #1, Dummy
     For KW = 1 To NW
        Input #1, N
        Input #1, NWD(N)
        NNW = NNW + NWD(N)
     Next KW
     '--- Block Material Data
     Input #1, Dummy: Input #1, Dummy
     '-------- Block Identifier / Material# (Default# is 1) --------
     For I = 1 To NSW: IDBLK(I) = 1: Next I
     Do
        Input #1, NTMP
        If NTMP = 0 Then Exit Do
        Input #1, IDBLK(NTMP)
        If NM < IDBLK(NTMP) Then NM = IDBLK(NTMP)
     Loop
     '-----------------   Block Corner Data  ---------------
     NSR = NS * (NW + 1): NWR = NW * (NS + 1)
     ReDim XB(NGN, 2), SR(NSR, 2), WR(NWR, 2)
     Input #1, Dummy: Input #1, Dummy
     Do
        Input #1, NTMP
        If NTMP = 0 Then Exit Do
        Input #1, XB(NTMP, 1)
        Input #1, XB(NTMP, 2)
     Loop
     '---------- Evaluate Mid-points of S-Sides -------------
     For I = 1 To NW + 1
        For J = 1 To NS
           IJ = (I - 1) * NS + J
           SR(IJ, 1) = 0.5 * (XB(IJ + I - 1, 1) + XB(IJ + I, 1))
           SR(IJ, 2) = 0.5 * (XB(IJ + I - 1, 2) + XB(IJ + I, 2))
        Next J
     Next I
     '---------- Evaluate Mid-points of W-Sides -------------
     For I = 1 To NW
        For J = 1 To NS + 1
           IJ = (I - 1) * (NS + 1) + J
           WR(IJ, 1) = 0.5 * (XB(IJ, 1) + XB(IJ + NS + 1, 1))
           WR(IJ, 2) = 0.5 * (XB(IJ, 2) + XB(IJ + NS + 1, 2))
        Next J
     Next I
     '------ Mid Points for Sides that are curved or graded ------
     Line Input #1, Dummy: Line Input #1, Dummy
     '--- S-Sides
     Do
        Input #1, NTMP
        If NTMP = 0 Then Exit Do
        Input #1, SR(NTMP, 1)
        Input #1, SR(NTMP, 2)
     Loop
     Line Input #1, Dummy
     '--- W-Sides
     Do
        Input #1, NTMP
        If NTMP = 0 Then Exit Do
        Input #1, WR(NTMP, 1)
        Input #1, WR(NTMP, 2)
     Loop
     '--------- Merging Sides ----------
     If NSJ > 0 Then
        Input #1, Dummy: Input #1, Dummy
        ReDim MERG(NSJ, 4)
        For I = 1 To NSJ
           Input #1, N
           Input #1, L1
           Input #1, L2
           Call SideDiv(L1, L2, IDIV1)
           Input #1, L3
           Input #1, L4
           Call SideDiv(L3, L4, IDIV2)
           If IDIV1 <> IDIV2 Then
              picBox.Print "#Div don't match. Check merge data."
              End
           End If
           MERG(I, 1) = L1: MERG(I, 2) = L2
           MERG(I, 3) = L3: MERG(I, 4) = L4
        Next I
     End If
     Close #1
End Sub
Private Sub GlobalNode()
     '------- Global Node Locations of Corner Nodes ---------
     NTMPI = 1
     For I = 1 To NW + 1
        If I = 1 Then IINC = 0 Else IINC = NNS * NWD(I - 1)
        NTMPI = NTMPI + IINC: NTMPJ = 0
        For J = 1 To NS + 1
           IJ = (NS + 1) * (I - 1) + J
           If J = 1 Then JINC = 0 Else JINC = NSD(J - 1)
           NTMPJ = NTMPJ + JINC
           NGCN(IJ) = NTMPI + NTMPJ
        Next J
     Next I
     '---------------- Node Point Array --------------------
     NNT = NNS * NNW
     ReDim NNAR(NNT)
     For I = 1 To NNT
        NNAR(I) = -1
     Next I
     '--------- Zero Non-Existing Node Locations ---------
     For KW = 1 To NW
        For KS = 1 To NS
           KSW = NS * (KW - 1) + KS
           If IDBLK(KSW) <= 0 Then
              '-------- Operation within an Empty Block --------
              K1 = (KW - 1) * (NS + 1) + KS: N1 = NGCN(K1)
              NS1 = 2: If KS = 1 Then NS1 = 1
              NW1 = 2: If KW = 1 Then NW1 = 1
              NS2 = NSD(KS) + 1
              If KS < NS Then
                 If IDBLK(KSW + 1) > 0 Then NS2 = NSD(KS)
              End If
              NW2 = NWD(KW) + 1
              If KW < NW Then
                 If IDBLK(KSW + NS) > 0 Then NW2 = NWD(KW)
              End If
              For I = NW1 To NW2
                 IN1 = N1 + (I - 1) * NNS
                 For J = NS1 To NS2
                    IJ = IN1 + J - 1
                    NNAR(IJ) = 0
                 Next J
              Next I
              ICT = 0
              If NS2 = NSD(KS) Or NW2 = NWD(KW) Then ICT = 1
              If KS = NS Or KW = NW Then ICT = 1
              If ICT = 0 Then
                 If IDBLK(KSW + NS + 1) > 0 Then NNAR(IJ) = -1
              End If
           End If
        Next KS
     Next KW
     '--------  Node Identification for Side Merging ------
     If NSJ > 0 Then
        For I = 1 To NSJ
           I1 = MERG(I, 1): I2 = MERG(I, 2)
           Call SideDiv(I1, I2, IDIV)
           IA1 = NGCN(I1): IA2 = NGCN(I2)
           IASTP = (IA2 - IA1) / IDIV
           I1 = MERG(I, 3): I2 = MERG(I, 4)
           Call SideDiv(I1, I2, IDIV)
           IB1 = NGCN(I1): IB2 = NGCN(I2)
           IBSTP = (IB2 - IB1) / IDIV
           IAA = IA1 - IASTP
           For IBB = IB1 To IB2 Step IBSTP
              IAA = IAA + IASTP
              If IBB = IAA Then NNAR(IAA) = -1 Else NNAR(IBB) = IAA
           Next IBB
        Next I
     End If
     '----------  Final Node Numbers in the Array  --------
     NODE = 0
     For I = 1 To NNT
        If NNAR(I) > 0 Then
           II = NNAR(I): NNAR(I) = NNAR(II)
        ElseIf NNAR(I) < 0 Then
           NODE = NODE + 1: NNAR(I) = NODE
        End If
     Next I
End Sub
Private Sub SideDiv(I1, I2, IDIV)
     '===========  Number of Divisions  for Side I1,I2  ===========
     IMIN = I1: IMAX = I2
     If IMIN > I2 Then
        IMIN = I2
        IMAX = I1
     End If
     If (IMAX - IMIN) = 1 Then
        IDIV = NGCN(IMAX) - NGCN(IMIN)
     Else
        IDIV = (NGCN(IMAX) - NGCN(IMIN)) / NNS
     End If
End Sub
Private Sub CoordConnect()
     '------------  Nodal Coordinates  ---------------
     NN = NODE: NELM = 0
     ReDim X(NN, 2), XP(8, 2), NOC(2 * NNT, NEN), MAT(2 * NNT)
     For KW = 1 To NW
        For KS = 1 To NS
        KSW = NS * (KW - 1) + KS
        If IDBLK(KSW) <> 0 Then
           '---------  Extraction of Block Data  ----------
           NODW = NGCN(KSW + KW - 1) - NNS - 1
           For JW = 1 To NWD(KW) + 1
              ETA = -1 + 2 * (JW - 1) / NWD(KW)
              NODW = NODW + NNS: NODS = NODW
              For JS = 1 To NSD(KS) + 1
                 XI = -1 + 2 * (JS - 1) / NSD(KS)
                 NODS = NODS + 1: NODE = NNAR(NODS)
                 Call BlockXY(KW, KSW)
                 Call Shape(XI, ETA)
                 For J = 1 To 2
                    C1 = 0
                    For I = 1 To 8
                       C1 = C1 + SH(I) * XP(I, J)
                    Next I
                    X(NODE, J) = C1
                 Next J
                 '-----------------  Connectivity  ----------------
                 If JS <> NSD(KS) + 1 And JW <> NWD(KW) + 1 Then
                    N1 = NODE: N2 = NNAR(NODS + 1)
                    N4 = NNAR(NODS + NNS): N3 = NNAR(NODS + NNS + 1)
                    NELM = NELM + 1
                    If NEN = 3 Then
                       '------------- Triangular Elements ------------
                       NOC(NELM, 1) = N1: NOC(NELM, 2) = N2
                       NOC(NELM, 3) = N3: MAT(NELM) = IDBLK(KSW)
                       NELM = NELM + 1: NOC(NELM, 1) = N3: NOC(NELM, 2) = N4
                       NOC(NELM, 3) = N1: MAT(NELM) = IDBLK(KSW)
                    Else
                       '------------- Quadrilateral Elements ----------
                       NOC(NELM, 1) = N1: NOC(NELM, 2) = N2: MAT(NELM) = IDBLK(KSW)
                       NOC(NELM, 3) = N3: NOC(NELM, 4) = N4
                    End If
                 End If
              Next JS
           Next JW
        End If
     Next KS
   Next KW
     NE = NELM
     If NEN = 3 Then
     '--------- Readjustment for Triangle Connectivity ----------
        NE2 = NE / 2
        For I = 1 To NE2
           I1 = 2 * I - 1: N1 = NOC(I1, 1): N2 = NOC(I1, 2)
           N3 = NOC(I1, 3): N4 = NOC(2 * I, 2)
           X13 = X(N1, 1) - X(N3, 1): Y13 = X(N1, 2) - X(N3, 2)
           X24 = X(N2, 1) - X(N4, 1): Y24 = X(N2, 2) - X(N4, 2)
           If (X13 * X13 + Y13 * Y13) > 1.1 * (X24 * X24 + Y24 * Y24) Then
              NOC(I1, 3) = N4: NOC(2 * I, 3) = N2
           End If
        Next I
     End If
End Sub
Private Sub BlockXY(KW, KSW)
     '======  Coordinates of 8-Nodes of the Block  ======
     N1 = KSW + KW - 1
     XP(1, 1) = XB(N1, 1): XP(1, 2) = XB(N1, 2)
     XP(3, 1) = XB(N1 + 1, 1): XP(3, 2) = XB(N1 + 1, 2)
     XP(5, 1) = XB(N1 + NS + 2, 1): XP(5, 2) = XB(N1 + NS + 2, 2)
     XP(7, 1) = XB(N1 + NS + 1, 1): XP(7, 2) = XB(N1 + NS + 1, 2)
     XP(2, 1) = SR(KSW, 1): XP(2, 2) = SR(KSW, 2)
     XP(6, 1) = SR(KSW + NS, 1): XP(6, 2) = SR(KSW + NS, 2)
     XP(8, 1) = WR(N1, 1): XP(8, 2) = WR(N1, 2)
     XP(4, 1) = WR(N1 + 1, 1): XP(4, 2) = WR(N1 + 1, 2)
End Sub
Private Sub Shape(XI, ETA)
     '==============  Shape Functions  ================
     SH(1) = -(1 - XI) * (1 - ETA) * (1 + XI + ETA) / 4
     SH(2) = (1 - XI * XI) * (1 - ETA) / 2
     SH(3) = -(1 + XI) * (1 - ETA) * (1 - XI + ETA) / 4
     SH(4) = (1 - ETA * ETA) * (1 + XI) / 2
     SH(5) = -(1 + XI) * (1 + ETA) * (1 - XI - ETA) / 4
     SH(6) = (1 - XI * XI) * (1 + ETA) / 2
     SH(7) = -(1 - XI) * (1 + ETA) * (1 + XI - ETA) / 4
     SH(8) = (1 - ETA * ETA) * (1 - XI) / 2
End Sub
Private Sub Output()
     '===== Print Displacements, Stresses, and Reactions
     cdlDialog.Filter = "All Files (*.*)|*.*|<Output file> FileName.out|*.out"
     cdlDialog.FileName = ""
     cdlDialog.ShowSave
     File2 = cdlDialog.FileName
     Open File2 For Output As #2
     Print #2, "Program MESHGEN - CHANDRUPATLA & BELEGUNDU"
     Print #2, Title
     NDIM = 2: NDN = 2
     Print #2, "NN  NE  NM  NDIM  NEN  NDN"
     Print #2, NN; NE; NM; NDIM; NEN; NDN
     Print #2, "ND   NL   NMPC"
     Print #2, ND; NL; NMPC
     Print #2, "Node#   X    Y"
     For I = 1 To NN
        Print #2, I;
        For J = 1 To NDIM
           Print #2, X(I, J);
        Next J
        Print #2,
     Next I
     Print #2, "Elem#  Node1  Node2  Node3";
     If NEN = 3 Then Print #2, "  Material#"
     If NEN = 4 Then Print #2, "  Node4  Material#"
     For I = 1 To NE
        Print #2, I;
        For J = 1 To NEN
           Print #2, NOC(I, J);
        Next J
        Print #2, MAT(I)
     Next I
     Close #2
     picBox.Print "Data has been stored in the file "; File2
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

