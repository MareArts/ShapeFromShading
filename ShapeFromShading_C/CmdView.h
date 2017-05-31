#if !defined(AFX_CMDVIEW_H__6DD4DDD7_4018_4024_A13D_043BA9E81BC0__INCLUDED_)
#define AFX_CMDVIEW_H__6DD4DDD7_4018_4024_A13D_043BA9E81BC0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// CmdView.h : header file
//

/////////////////////////////////////////////////////////////////////////////
// CCmdView form view
#include "FileIO.h"


#ifndef __AFXEXT_H__
#include <afxext.h>
#endif



class CCmdView : public CFormView
{
protected:
	CCmdView();           // protected constructor used by dynamic creation
	DECLARE_DYNCREATE(CCmdView)

// Form Data
public:
	//{{AFX_DATA(CCmdView)
	enum { IDD = CMD_VIEW };
	CProgressCtrl	m_Progress;
	//}}AFX_DATA


	BOOL SAM;
	//ÀÌ¹ÌÁö 
	Matrix<float> Img1;
	Matrix<float> Img2;
	Matrix<float> Img3;
	Matrix<float> Img4;
	Matrix<float> Img5;
	Matrix<float> Img6;

	//±¤¿ø
	Matrix<float> Lgt1;
	Matrix<float> Lgt2;
	Matrix<float> Lgt3;
	Matrix<float> Lgt4;
	Matrix<float> Lgt5;
	Matrix<float> Lgt6;

	
	//S º¤ÅÍ ±¤¿ø 
	Matrix<float> S;
	


// Attributes
public:

// Operations
public:
	CView* GetRenderView();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CCmdView)
	public:
	virtual void OnInitialUpdate();
	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	//}}AFX_VIRTUAL

// Implementation
protected:
	virtual ~CCmdView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

	// Generated message map functions
	//{{AFX_MSG(CCmdView)
	afx_msg void OnSample();
	afx_msg void OnCal();
	afx_msg void OnRadio1();
	afx_msg void OnRadio2();
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_CMDVIEW_H__6DD4DDD7_4018_4024_A13D_043BA9E81BC0__INCLUDED_)
