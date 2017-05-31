// CmdView.cpp : implementation file
//

#include "stdafx.h"
#include "MainFrm.h"
#include "HomeWork2View.h"
#include "HomeWork2.h"
#include "CmdView.h"
#include "math.h"
#include "limits.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CCmdView

IMPLEMENT_DYNCREATE(CCmdView, CFormView)

CCmdView::CCmdView()
	: CFormView(CCmdView::IDD)
{
	//{{AFX_DATA_INIT(CCmdView)
	//}}AFX_DATA_INIT

	SAM=FALSE;
	
}

CCmdView::~CCmdView()
{
	

}

void CCmdView::DoDataExchange(CDataExchange* pDX)
{
	CFormView::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CCmdView)
	DDX_Control(pDX, IDC_PROGRESS1, m_Progress);
	//}}AFX_DATA_MAP
}


BEGIN_MESSAGE_MAP(CCmdView, CFormView)
	//{{AFX_MSG_MAP(CCmdView)
	ON_BN_CLICKED(IDC_SAMPLE, OnSample)
	ON_BN_CLICKED(IDC_CAL, OnCal)
	ON_BN_CLICKED(IDC_RADIO1, OnRadio1)
	ON_BN_CLICKED(IDC_RADIO2, OnRadio2)
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CCmdView diagnostics

#ifdef _DEBUG
void CCmdView::AssertValid() const
{
	CFormView::AssertValid();
}

void CCmdView::Dump(CDumpContext& dc) const
{
	CFormView::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CCmdView message handlers

CView* CCmdView::GetRenderView()
{
	CHomeWork2App *pApp = (CHomeWork2App *)AfxGetApp();
	CMainFrame *pMainFrame = (CMainFrame *)pApp->m_pMainWnd;
	pMainFrame->GetActiveFrame();
	CView *pView = (CView *)pMainFrame->m_wndSplitter.GetPane(0,1);

	return pView;
}

 void CCmdView::OnRadio1() 
 {
 	// TODO: Add your control notification handler code here
 	CHomeWork2View *pView = (CHomeWork2View *)GetRenderView();
 	pView->SelectModel=0;
 	pView->Invalidate();
 }

 void CCmdView::OnRadio2() 
 {
 	// TODO: Add your control notification handler code here
 	CHomeWork2View *pView = (CHomeWork2View *)GetRenderView();
 	pView->SelectModel=1;
 	pView->Invalidate();
 }



void CCmdView::OnInitialUpdate() 
{
	CFormView::OnInitialUpdate();
	
	// TODO: Add your specialized code here and/or call the base class
	((CButton*)GetDlgItem(IDC_RADIO1))->SetCheck(1);
}



void CCmdView::OnSample() 
{
	// TODO: Add your control notification handler code here

	//이미지와 라이트 매트릭스 크기 설정
	Img1.Resize(240,320);
	Img2.Resize(240,320);
	Img3.Resize(240,320);
	Img4.Resize(240,320);
	Img5.Resize(240,320);
	Img6.Resize(240,320);

	Lgt1.Resize(1,3);
	Lgt2.Resize(1,3);
	Lgt3.Resize(1,3);
	Lgt4.Resize(1,3);
	Lgt5.Resize(1,3);
	Lgt6.Resize(1,3);

	
		
//구
	CFileIO A("IMG3_1.bmp");
	CFileIO B("IMG3_2.bmp");
	CFileIO C("IMG3_3.bmp");
	CFileIO D("IMG3_4.bmp");
	CFileIO E("IMG3_5.bmp");
	CFileIO F("IMG3_6.bmp");
/*
	//주전자 
    CFileIO A("IMG1.bmp");
	CFileIO B("IMG2.bmp");
	CFileIO C("IMG3.bmp");
	CFileIO D("IMG4.bmp");
	CFileIO E("IMG5.bmp");
	CFileIO F("IMG6.bmp");

    //원뿔
	CFileIO A("IMG2_1.bmp");
	CFileIO B("IMG2_2.bmp");
	CFileIO C("IMG2_3.bmp");
	CFileIO D("IMG2_4.bmp");
	CFileIO E("IMG2_5.bmp");
	CFileIO F("IMG2_6.bmp");

    //구
	CFileIO A("IMG3_1.bmp");
	CFileIO B("IMG3_2.bmp");
	CFileIO C("IMG3_3.bmp");
	CFileIO D("IMG3_4.bmp");
	CFileIO E("IMG3_5.bmp");
	CFileIO F("IMG3_6.bmp");

    //교수님 모델 
	CFileIO A("IMG4_1.bmp");
	CFileIO B("IMG4_2.bmp");
	CFileIO C("IMG4_3.bmp");
	CFileIO D("IMG4_4.bmp");
	CFileIO E("IMG4_5.bmp");
	CFileIO F("IMG4_6.bmp");

	*/
	int i,j;

	for(i=0; i<240; ++i)
	{
		for(j=0; j<320; ++j)
		{
			//밝기값 매트릭스에 넣기
			Img1(i+1,j+1) = A.orgImg[i*320+j];
			Img2(i+1,j+1) = B.orgImg[i*320+j];
			Img3(i+1,j+1) = C.orgImg[i*320+j];
			Img4(i+1,j+1) = D.orgImg[i*320+j];
			Img5(i+1,j+1) = E.orgImg[i*320+j];
			Img6(i+1,j+1) = F.orgImg[i*320+j];

		}
	}
	
	
	//광원 위치 읽기 
	double L1,L2,L3;
	ifstream fp;
	fp.open("Light.txt");

	fp >> L1 >> L2 >> L3;
	Lgt1(1,1) = L1; Lgt1(1,2) = L2; Lgt1(1,3) = L3;
	fp >> L1 >> L2 >> L3;
	Lgt2(1,1) = L1; Lgt2(1,2) = L2; Lgt2(1,3) = L3;
	fp >> L1 >> L2 >> L3;
	Lgt3(1,1) = L1; Lgt3(1,2) = L2; Lgt3(1,3) = L3;
	fp >> L1 >> L2 >> L3;
	Lgt4(1,1) = L1; Lgt4(1,2) = L2; Lgt4(1,3) = L3;
	fp >> L1 >> L2 >> L3;
	Lgt5(1,1) = L1; Lgt5(1,2) = L2; Lgt5(1,3) = L3;
	fp >> L1 >> L2 >> L3;
	Lgt6(1,1) = L1; Lgt6(1,2) = L2; Lgt6(1,3) = L3;

	fp.close();

	SAM=TRUE;

	/*
	ofstream out;
	out.open("Sample.txt");
	for(i=1; i<=240; ++i)
	{
		for(j=1; j<=320; ++j)
		{
			out << S(i,j) << " ";
		}
		out << endl;
	}
	out.close();
	*/

}



void CCmdView::OnCal() 
{
	// TODO: Add your control notification handler code here
	
	//if loading 이 되었을때?

	if(SAM)
	{

	
	int i,j;

	//프로그래서바 
	m_Progress.SetRange(0,100);
	int iPos=0;

	//법선 벡터 Nx,Ny,Nz와 깊이 정보 Z를 구한다.
	//광원을 단위 벡터로 만든다. 
	Lgt1 = Lgt1 / NormF(Lgt1);
	Lgt2 = Lgt2 / NormF(Lgt2);
	Lgt3 = Lgt3 / NormF(Lgt3);
	Lgt4 = Lgt4 / NormF(Lgt4);
	Lgt5 = Lgt5 / NormF(Lgt5);
	Lgt6 = Lgt6 / NormF(Lgt6);
	

	//6 by 3의 하나의 S 메트릭스를 만든다. 
	S.Resize(6,3);
	S(1,1)=Lgt1(1,1); S(1,2) = Lgt1(1,2); S(1,3) = Lgt1(1,3);
	S(2,1)=Lgt2(1,1); S(2,2) = Lgt2(1,2); S(2,3) = Lgt2(1,3);
	S(3,1)=Lgt3(1,1); S(3,2) = Lgt3(1,2); S(3,3) = Lgt3(1,3);
	S(4,1)=Lgt4(1,1); S(4,2) = Lgt4(1,2); S(4,3) = Lgt4(1,3);
	S(5,1)=Lgt5(1,1); S(5,2) = Lgt5(1,2); S(5,3) = Lgt5(1,3);
	S(6,1)=Lgt6(1,1); S(6,2) = Lgt6(1,2); S(6,3) = Lgt6(1,3);

	//노멀 벡터 
	Matrix<float> Nx;
	//노멀 벡터 
	Matrix<float> Ny;
	//노멀 벡터 
	Matrix<float> Nz;
	//깊이 벡터 
	Matrix<float> Z;
	double ZMin = LONG_MAX;
	double ZMax = LONG_MIN;
	double xMax = LONG_MIN;
	double yMax = LONG_MIN;
	double zMax = LONG_MIN;

	//E 벡터 모델 밝기 
	Matrix<float> E(6,1); //밝기 벡터 
	Matrix<float> B; 
	Matrix<float> sB; 
	Matrix<float> p(240,320);
	Matrix<float> q(240,320);
	
	Nx.Resize(240,320);
	Ny.Resize(240,320);
	Nz.Resize(240,320);
	Z.Resize(240,320);

	//계산을 위한 임시 변수들	
	Matrix<float> LN(240,320);
	Matrix<float> tM(1,3);
	double tV;

	
	//b 구하는 공식 
	//b = (inv(S'*S))*S'*E;
	//미리 E 연산하기 전까지 구해놓음	
	sB= Pinv(Transpose(S)*S)*Transpose(S);

	
	for(i=0; i<240; i++)
	{
		for(j=0; j<320; j++)
		{
			//모델의 밝기를 6by1 로 넣는다.
			E(1,1) = Img1(i+1,j+1); E(2,1) = Img2(i+1,j+1); E(3,1) = Img3(i+1,j+1);
			E(4,1) = Img4(i+1,j+1); E(5,1) = Img5(i+1,j+1); E(6,1) = Img6(i+1,j+1);

			//psuedo code 를 푼다.
			//B = Inv(Transpose(S)*S)*Transpose(S)*E;
			//속도 개선을 위해 sB로 앞부분 구해놓음
			B = sB * E;
			//Albedo를 나누어 Normal 벡터로 만든다.
			//0인지 조사 
			tV = NormF(B);
			if(tV == 0)
			{
				Nx(i+1,j+1)=0; Ny(i+1,j+1)=0; Nz(i+1,j+1)=0;
			}else{
				B = B/tV;
				Nx(i+1,j+1)=B(1,1); Ny(i+1,j+1)=B(2,1); Nz(i+1,j+1)=B(3,1);
			}

			//Depth를 구한다. 
			tM(1,1) = Nx(i+1,j+1); tM(1,2) = Ny(i+1,j+1); tM(1,3) = Nz(i+1,j+1);
			//Z를 위한 pq를 구한다. 
			tV = NormF(tM);
			if(tV == 0)
			{
				tM(1,1) = 0; tM(1,2) = 0; tM(1,3) = 0;
			}else{
				tM = tM/tV;
			}

			p(i+1,j+1) = tM(1,1);
			q(i+1,j+1) = tM(1,2);

			//Max Min구하기
			if(xMax < Nx(i+1,j+1) )
				xMax = Nx(i+1,j+1);
			if(yMax < Ny(i+1,j+1) )
				yMax = Ny(i+1,j+1);
			if(yMax < Nz(i+1,j+1) )
				yMax = Nz(i+1,j+1);


			
			m_Progress.SetPos(iPos/1536);
			iPos++;	
		}
	}
	
	
	//Z를 구함 
	double pS,qS;
	int x;
	for(i=1; i<=240; i++)
	{
		for(j=1; j<=320; ++j)
		{
			pS=qS=0;
			//(0, V)를 누적 
			for(x=1; x<=i; ++x)
			{
				qS += q(x,1);
			}

			//(U,V)를 누적 
			for(x=1; x<=j; ++x)
			{
				pS += p(i,x);
			}
			Z(i,j) = pS+qS;

			//Max Min구하기
			if(ZMax < Z(i,j) )
				ZMax = Z(i,j);
			

			m_Progress.SetPos(iPos/1536);
			iPos++;	
		}		
	}

	CHomeWork2View *pView = (CHomeWork2View *)GetRenderView();
 	pView->Nx=  Nx;
	pView->Ny=  Ny;
	pView->Nz=  Nz;
	pView->Z=  Z;
	pView->b_RENDER = TRUE;
	pView->m_ZMax = ZMax;
	pView->m_xMax = xMax;
	pView->m_yMax = yMax;
	pView->m_zMax = zMax;
	pView->Invalidate();

	m_Progress.SetPos(0);

	}else{
		::AfxMessageBox("File is not Load!!");
	}
	/*
	//txt로 결과 확인 매트랩과 같았음 	
	ofstream out;
	out.open("S.txt");
	for(i=1; i<=6; ++i)
	{
		for(j=1; j<=3; ++j)
		{
			out << S(i,j) << " ";
		}
		out << endl;
	}
	out.close();

	out.open("Lgt.txt");
	for(j=1; j<=3; ++j)
	{
		out << Lgt1(1,j) << " ";
		out << Lgt2(1,j) << " ";
		out << Lgt3(1,j) << " ";
		out << Lgt4(1,j) << " ";
		out << Lgt5(1,j) << " ";
		out << Lgt6(1,j) << endl;

	}
	
	out.close();
	

	out.open("sB.txt");
	for(i=1; i<=3; ++i)
	{
		for(j=1; j<=6; ++j)
		{
			out << sB(i,j) << " ";
		}
		out << endl;
	}
	out.close();
	

	out.open("out.txt");

	for(i=0; i<240; ++i)
	{
		for(j=0; j<320; ++j)
		{
			out << Z(i+1,j+1) << " ";
		}
		out << endl;
	}
	*/
	
	
}

