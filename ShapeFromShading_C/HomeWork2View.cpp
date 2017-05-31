// HomeWork2View.cpp : implementation of the CHomeWork2View class
//

#include "stdafx.h"
#include "HomeWork2.h"

#include "math.h"
#include "HomeWork2Doc.h"
#include "HomeWork2View.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2View

IMPLEMENT_DYNCREATE(CHomeWork2View, CView)

BEGIN_MESSAGE_MAP(CHomeWork2View, CView)
	//{{AFX_MSG_MAP(CHomeWork2View)
	ON_WM_CREATE()
	ON_WM_SIZE()
	ON_WM_LBUTTONDOWN()
	ON_WM_LBUTTONUP()
	ON_WM_MOUSEMOVE()
	ON_WM_MOUSEWHEEL()
	ON_WM_RBUTTONDOWN()
	ON_WM_RBUTTONUP()
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2View construction/destruction

GLvoid *font_style2 = GLUT_BITMAP_HELVETICA_12;
void printw2 (float x, float y, float z, char* format, ...)
{
    va_list arg_list;
    char str[256];
	int i;
    
    va_start(arg_list, format);
    vsprintf(str, format, arg_list);
    va_end(arg_list);
    
    glRasterPos3f (x, y, z);

    for (i = 0; str[i] != '\0'; i++)
        glutBitmapCharacter(font_style2, str[i]);
}

CHomeWork2View::CHomeWork2View()
{
	// TODO: add construction code here
	b_RENDER=FALSE;
}

CHomeWork2View::~CHomeWork2View()
{
}

BOOL CHomeWork2View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2View drawing

void CHomeWork2View::OnDraw(CDC* pDC)
{
	CHomeWork2Doc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	// TODO: add draw code for native data here

	wglMakeCurrent(pDC->GetSafeHdc(),m_hGLContext);	//MFC에서 OpenGL을 뿌리기 위한 설정

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); //모드 
	glMatrixMode(GL_MODELVIEW); //모드 


	//Rendering 
	RenderDisplay(); //드로잉 
	// Double buffer
	SwapBuffers(pDC->GetSafeHdc()); //더블버퍼링 설정
	
}

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2View diagnostics

#ifdef _DEBUG
void CHomeWork2View::AssertValid() const
{
	CView::AssertValid();
}

void CHomeWork2View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CHomeWork2Doc* CHomeWork2View::GetDocument() // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CHomeWork2Doc)));
	return (CHomeWork2Doc*)m_pDocument;
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2View message handlers



int CHomeWork2View::OnCreate(LPCREATESTRUCT lpCreateStruct) 
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;
	
	// TODO: Add your specialized creation code here

	//////////////////////////////////////////////////////////////////////////
	//MFC에서 OpenGL을 쓰기 위한 설정 
	HWND hWnd = GetSafeHwnd();
	HDC hDC = ::GetDC(hWnd);	
	if(SetWindowPixelFormat(hDC)==FALSE)
		return 0;
	if(CreateViewGLContext(hDC)==FALSE)
		return 0;

	//////////////////////////////////////////////////////////////////////////
	// Default mode
	m_xScaling = 1.0f;
	m_yScaling = 1.0f;
	m_zScaling = 1.0f;
	
	
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA | GLUT_DEPTH);
	glClearColor(0.0, 0.0, 0.0, 0.0);		
	
	//초기값 
	InitValue();
	//광원 설정
	InitLight();
	
	return 0;
}

BOOL CHomeWork2View::SetWindowPixelFormat(HDC hDC)
{
	PIXELFORMATDESCRIPTOR pixelDesc;
	
	pixelDesc.nSize = sizeof(PIXELFORMATDESCRIPTOR);
	pixelDesc.nVersion = 1;
	
	pixelDesc.dwFlags = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL |
		PFD_DOUBLEBUFFER | PFD_STEREO_DONTCARE;
	
	pixelDesc.iPixelType = PFD_TYPE_RGBA;
	pixelDesc.cColorBits = 32;
	pixelDesc.cRedBits = 8;
	pixelDesc.cRedShift = 16;
	pixelDesc.cGreenBits = 8;
	pixelDesc.cGreenShift = 8;
	pixelDesc.cBlueBits = 8;
	pixelDesc.cBlueShift = 0;
	pixelDesc.cAlphaBits = 0;
	pixelDesc.cAlphaShift = 0;
	pixelDesc.cAccumBits = 64;
	pixelDesc.cAccumRedBits = 16;
	pixelDesc.cAccumGreenBits = 16;
	pixelDesc.cAccumBlueBits = 16;
	pixelDesc.cAccumAlphaBits = 0;
	pixelDesc.cDepthBits = 32;
	pixelDesc.cStencilBits = 8;
	pixelDesc.cAuxBuffers = 0;
	pixelDesc.iLayerType = PFD_MAIN_PLANE;
	pixelDesc.bReserved = 0;
	pixelDesc.dwLayerMask = 0;
	pixelDesc.dwVisibleMask = 0;
	pixelDesc.dwDamageMask = 0;

	m_GLPixelIndex = ChoosePixelFormat(hDC,&pixelDesc);
	if(m_GLPixelIndex == 0) // Choose default
	{
		m_GLPixelIndex = 1;
		if(DescribePixelFormat(hDC,m_GLPixelIndex,
			sizeof(PIXELFORMATDESCRIPTOR),&pixelDesc)==0)
			return FALSE;
	}
	
	if(!SetPixelFormat(hDC,m_GLPixelIndex,&pixelDesc))
		return FALSE;
	
	return TRUE;
}

BOOL CHomeWork2View::CreateViewGLContext(HDC hDC)
{
	m_hGLContext = wglCreateContext(hDC);
	
	if(m_hGLContext==NULL)
		return FALSE;
	
	if(wglMakeCurrent(hDC,m_hGLContext)==FALSE)
		return FALSE;	

	return TRUE;
}

void CHomeWork2View::OnSize(UINT nType, int cx, int cy) 
{
	CView::OnSize(nType, cx, cy);
	
	// TODO: Add your message handler code here
	CSize size(cx,cy);
	double aspect;
	aspect = (cy == 0) ? (double)size.cx : (double)size.cx/(double)size.cy;

	glViewport(0, 0, (GLsizei)size.cx, (GLsizei)size.cy); // 윈도우 사이즈 비율로 뿌림
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); //항등원 로딩 
	glOrtho(-1.0, 1.0, -1.0, 1.0, -10.0, 10.0); 
	//glOrtho (GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, 
	//GLdouble zNear, GLdouble zFar)
	//시점의 위치와 가시 거리의 범위를 설정
	/*
	glViewport(0,0,size.cx,size.cy);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,aspect,1,15.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glDrawBuffer(GL_BACK);
	glEnable(GL_DEPTH_TEST);
		*/

}

void CHomeWork2View::InitLight()
{

	
// Default mode
	glPolygonMode(GL_FRONT,GL_LINE);
	glPolygonMode(GL_BACK,GL_LINE);
	glShadeModel(GL_FLAT);
	glEnable(GL_NORMALIZE);	
	
	
	
	// Lights, material properties
	GLfloat	ambientProperties[]  = {0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat	diffuseProperties[]  = {0.8f, 0.8f, 0.8f, 1.0f};
	GLfloat	specularProperties[] = {1.0f, 1.0f, 1.0f, 1.0f};
	GLfloat	light_position[] = {1.0f, 10.0f, 1.0f, 0.0f};
	
	glClearDepth( 1.0 );	
	glLightfv( GL_LIGHT0, GL_AMBIENT, ambientProperties);
	glLightfv( GL_LIGHT0, GL_DIFFUSE, diffuseProperties);
	glLightfv( GL_LIGHT0, GL_SPECULAR, specularProperties);
	glLightfv( GL_LIGHT0, GL_POSITION, light_position);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 1.0);
	
	// Default : lighting
	glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHTING);
	glDisable(GL_LIGHTING);

	//glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);	

}

void CHomeWork2View::InitValue()
{
	m_diffuse =  0.0; //확산광 설정값 
	m_specular = 0.0; //반사광 설정값
	m_ambient = 0.0;  //주변광 설정값
	m_surface = 25.0; //표면 매끄럼 정도값
	m_Position = 0; // 광원 위치 0, 1이 있음

	//Left 마우스
	m_LeftButtonDown = FALSE;
	m_RightButtonDown = FALSE;

	m_xRotation = 0.0f;
	m_yRotation = 0.0f;
	m_zRotation = 0.0f;


	//Model Option
	SelectModel=0;
}

void CHomeWork2View::RenderDisplay()
{	
	//void APIENTRY gluLookAt 
	//( GLdouble eyex, GLdouble eyey, GLdouble eyez,
	//GLdouble centerx, GLdouble centery, GLdouble centerz, GLdouble upx, GLdouble upy, GLdouble upz)

	glDisable(GL_LIGHTING); //글자에 광원이 적용 되지 않도록 잠시 끈다.
	char buffer[10]="asdf";
	CString me;	
	me.Format("Lighting Ambient : %f",m_ambient); // 글자 만들기 
	printw2(-0.9,+0.9,0, (LPSTR)(LPCTSTR)me); //글자 뿌리기 
	me.Format("Lighting Diffuse: %f",m_diffuse); // 글자 만들기 
	printw2(-0.9,+0.85,0, (LPSTR)(LPCTSTR)me);//글자 뿌리기 
	me.Format("Lighting Specular : %f",m_specular);// 글자 만들기 
	printw2(-0.9,+0.8,0, (LPSTR)(LPCTSTR)me);//글자 뿌리기 
	me.Format("Material Shininess : %f",m_surface);// 글자 만들기 
	printw2(-0.9,+0.75,0, (LPSTR)(LPCTSTR)me);//글자 뿌리기 
	glEnable(GL_LIGHTING); //다시 광원을 켠다.
	
	glLoadIdentity();	
	glLoadIdentity();	
	//Rotation	
	glPushMatrix(); //셋팅 스택에 푸쉬
	glRotatef(m_xRotation, 1.0, 0.0, 0.0); //마우스로 돌린 로테이션 적용
	glRotatef(m_yRotation, 0.0, 1.0, 0.0);//마우스로 돌린 로테이션 적용
	glRotatef(m_zRotation, 0.0, 0.0, 1.0);//마우스로 돌린 로테이션 적용
	glScalef(m_xScaling,m_yScaling,m_zScaling);
//	gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0); //카메라 시점 설정 
	

	int a,b;

	int sft = 4;
	if(b_RENDER)
	{		
		if(SelectModel==0)
		{
			//법선 만들기 
			//깊이 만들기
			int x,y;
			x=0;
			
			for(int i=-120; i<120; i+=sft)
			{
				x+=sft;
				y=0;
				
				for(int j=-160; j<160; j+=sft)
				{
					y+=sft;

					glColor3f( -1*Z(x,y)/m_ZMax , Z(x,y)/m_ZMax ,0.1);
					glBegin(GL_LINES);					
					glVertex3f(i/120.0, j/160.0, 0);
					glVertex3f(i/120.0 + Nx(x,y)/5 , j/160.0 + Ny(x,y)/5 , Nz(x,y)/5);

					glEnd();
				}
				
			}

			
		}
		else if(SelectModel==1)
		{
			//깊이 만들기
			int x,y;
			x=0;
			
			for(int i=-120; i<120; i+=sft)
			{
				x+=sft;
				y=0;
				glBegin(GL_LINE_STRIP); 
				for(int j=-160; j<160; j+=sft)
				{
					y+=sft;
					glColor3f( -1*Z(x,y)/m_ZMax , Z(x,y)/m_ZMax ,0.1);
					glVertex3f(i/120.0, j/160.0, -1*Z(x,y)/m_ZMax*0.1);
				}
				glEnd();
			}

			x=y=0;
			for(int j=-160; j<160; j+=sft)
			{
				y+=sft;
				x=0;
				glBegin(GL_LINE_STRIP); 
				for(int i=-120; i<120; i+=sft)
				{
					
					x+=sft;
					glColor3f( -1*Z(x,y)/m_ZMax ,Z(x,y)/m_ZMax,0.1);
					glVertex3f(i/120.0, j/160.0, -1*Z(x,y)/m_ZMax*0.1);
				}
				glEnd();
			}
		}
	}else{

		float width = 1.5f;
		float step = 0.1f;
		
		float x = 0.0f;
		float y = 0.0f;
		float z = 0.0f;
		
		

	
		for(x = -width; x < width; x += step)
		{
			glBegin(GL_LINE_STRIP); 
			for(z = -width; z < width; z += step)
			{
				float r,g,b;
				g = 0.0f;
				
				y = (float)sin(3*x)*cos(3*z)/3.0f;
				r = y*3.0f;
				b = 1.0f-y;
				glColor3f(r,g,b);
				glVertex3d(x,y,z);
				
				y = (float)sin(3*(x+step))*cos(3*z)/3.0f;
				r = y*3.0f;
				b = 1.0f-y;
				glColor3f(r,g,b);
				glVertex3d(x+step,y,z);
			}
			glEnd();
		}
		
	}
		
	glPopMatrix(); //세팅을 다시 팝해서 없앰 
	glFlush(); //드로잉 재개 
}

void CHomeWork2View::OnLButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_LeftButtonDown = TRUE;
	m_LeftDownPos = point;
	CView::OnLButtonDown(nFlags, point);
}

void CHomeWork2View::OnLButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_LeftButtonDown = FALSE;
	CView::OnLButtonUp(nFlags, point);
}

void CHomeWork2View::OnMouseMove(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	if(m_LeftButtonDown)
	{	
		//	CMD_FormCommandView* pFormView ;		
		m_yRotation -= (float)(m_LeftDownPos.x - point.x)/3.0f;
		m_xRotation -= (float)(m_LeftDownPos.y - point.y)/3.0f;
		m_LeftDownPos = point;
		InvalidateRect(NULL,FALSE);
		
	}

	if(m_RightButtonDown)
	{	
		//	CMD_FormCommandView* pFormView ;		
		m_zRotation -= (float)(m_RightDownPos.x - point.x)/3.0f;		
		m_RightDownPos = point;
		InvalidateRect(NULL,FALSE);
		
	}

	CView::OnMouseMove(nFlags, point);
}






BOOL CHomeWork2View::OnMouseWheel(UINT nFlags, short zDelta, CPoint pt) 
{
	// TODO: Add your message handler code here and/or call default
	if(zDelta>0)
	{
		m_xScaling += 0.01;
		m_yScaling += 0.01;
		m_zScaling += 0.01;
		
	}else{
		if(m_xScaling-0.01 > 0){		
			m_xScaling -= 0.01;
			m_yScaling -= 0.01;
			m_zScaling -= 0.01;
		}		
	}
	InvalidateRect(NULL,FALSE);	

	return CView::OnMouseWheel(nFlags, zDelta, pt);
}

void CHomeWork2View::OnRButtonDown(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_RightButtonDown = TRUE;
	m_RightDownPos = point;
	CView::OnRButtonDown(nFlags, point);
}

void CHomeWork2View::OnRButtonUp(UINT nFlags, CPoint point) 
{
	// TODO: Add your message handler code here and/or call default
	m_RightButtonDown = FALSE;
	CView::OnRButtonUp(nFlags, point);
}
