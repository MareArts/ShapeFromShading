// HomeWork2.h : main header file for the HOMEWORK2 application
//

#if !defined(AFX_HOMEWORK2_H__2532957C_87FF_4332_8929_921483ED3B1F__INCLUDED_)
#define AFX_HOMEWORK2_H__2532957C_87FF_4332_8929_921483ED3B1F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2App:
// See HomeWork2.cpp for the implementation of this class
//

class CHomeWork2App : public CWinApp
{
public:
	CHomeWork2App();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CHomeWork2App)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation
	//{{AFX_MSG(CHomeWork2App)
	afx_msg void OnAppAbout();
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_HOMEWORK2_H__2532957C_87FF_4332_8929_921483ED3B1F__INCLUDED_)
