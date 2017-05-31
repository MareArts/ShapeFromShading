// HomeWork2View.h : interface of the CHomeWork2View class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_HOMEWORK2VIEW_H__DEDD2375_65C0_48C4_BE61_4B7FBB5E6979__INCLUDED_)
#define AFX_HOMEWORK2VIEW_H__DEDD2375_65C0_48C4_BE61_4B7FBB5E6979__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class CHomeWork2Doc;

class CHomeWork2View : public CView
{
protected: // create from serialization only
	CHomeWork2View();
	DECLARE_DYNCREATE(CHomeWork2View)

// Attributes
public:
	CHomeWork2Doc* GetDocument();

	//노멀 벡터 
	Matrix<float> Nx;
	//노멀 벡터 
	Matrix<float> Ny;
	//노멀 벡터 
	Matrix<float> Nz;
	//깊이 벡터 
	Matrix<float> Z;

	double m_ZMax;
	
	double m_xMax;	
	double m_yMax;
	double m_zMax;
	

	BOOL b_RENDER;

	float m_diffuse ; //확산광 설정값 
	float m_specular; //반사광 설정값
	float m_ambient ; //주변광 설정값
	float m_surface ; //표면 매끄럼 정도값	
	int m_Position   ;  // 광원 위치 0, 1, 2, 3이 있음 
	

	float m_xScaling;
	float m_yScaling;
	float m_zScaling;

	// Mouse 
	//button click
	BOOL m_LeftButtonDown;
	BOOL m_RightButtonDown;
	CPoint m_LeftDownPos;
	CPoint m_RightDownPos;
	HCURSOR m_CursorRotation;
	// 회전
	float m_xRotation;
	float m_yRotation;
	float m_zRotation;

	//Option 1
	int SelectModel; //0 ,1, 2, 3


// Operations
public:
	//OpenGL SettingHGLRC m_hGLContext;
	int m_GLPixelIndex;
	HGLRC m_hGLContext;
// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CHomeWork2View)
	public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
	//}}AFX_VIRTUAL

// Implementation
public:
	void RenderDisplay();
	void InitValue();
	void InitLight();
	BOOL CreateViewGLContext(HDC hDC);
	BOOL SetWindowPixelFormat(HDC hDC);
	virtual ~CHomeWork2View();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CHomeWork2View)
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnSize(UINT nType, int cx, int cy);
	afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnLButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnMouseMove(UINT nFlags, CPoint point);
	afx_msg BOOL OnMouseWheel(UINT nFlags, short zDelta, CPoint pt);
	afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in HomeWork2View.cpp
inline CHomeWork2Doc* CHomeWork2View::GetDocument()
   { return (CHomeWork2Doc*)m_pDocument; }
#endif

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_HOMEWORK2VIEW_H__DEDD2375_65C0_48C4_BE61_4B7FBB5E6979__INCLUDED_)
