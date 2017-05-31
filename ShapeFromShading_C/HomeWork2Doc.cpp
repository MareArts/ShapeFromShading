// HomeWork2Doc.cpp : implementation of the CHomeWork2Doc class
//

#include "stdafx.h"
#include "HomeWork2.h"

#include "HomeWork2Doc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2Doc

IMPLEMENT_DYNCREATE(CHomeWork2Doc, CDocument)

BEGIN_MESSAGE_MAP(CHomeWork2Doc, CDocument)
	//{{AFX_MSG_MAP(CHomeWork2Doc)
		// NOTE - the ClassWizard will add and remove mapping macros here.
		//    DO NOT EDIT what you see in these blocks of generated code!
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2Doc construction/destruction

CHomeWork2Doc::CHomeWork2Doc()
{
	// TODO: add one-time construction code here

}

CHomeWork2Doc::~CHomeWork2Doc()
{
}

BOOL CHomeWork2Doc::OnNewDocument()
{
	if (!CDocument::OnNewDocument())
		return FALSE;

	// TODO: add reinitialization code here
	// (SDI documents will reuse this document)

	return TRUE;
}



/////////////////////////////////////////////////////////////////////////////
// CHomeWork2Doc serialization

void CHomeWork2Doc::Serialize(CArchive& ar)
{
	if (ar.IsStoring())
	{
		// TODO: add storing code here
	}
	else
	{
		// TODO: add loading code here
	}
}

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2Doc diagnostics

#ifdef _DEBUG
void CHomeWork2Doc::AssertValid() const
{
	CDocument::AssertValid();
}

void CHomeWork2Doc::Dump(CDumpContext& dc) const
{
	CDocument::Dump(dc);
}
#endif //_DEBUG

/////////////////////////////////////////////////////////////////////////////
// CHomeWork2Doc commands
