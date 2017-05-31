// HomeWork2Doc.h : interface of the CHomeWork2Doc class
//
/////////////////////////////////////////////////////////////////////////////

#if !defined(AFX_HOMEWORK2DOC_H__49E69189_7C93_44A8_B9CC_F0397EC463DC__INCLUDED_)
#define AFX_HOMEWORK2DOC_H__49E69189_7C93_44A8_B9CC_F0397EC463DC__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class CHomeWork2Doc : public CDocument
{
protected: // create from serialization only
	CHomeWork2Doc();
	DECLARE_DYNCREATE(CHomeWork2Doc)

// Attributes
public:

// Operations
public:

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CHomeWork2Doc)
	public:
	virtual BOOL OnNewDocument();
	virtual void Serialize(CArchive& ar);
	//}}AFX_VIRTUAL

// Implementation
public:
	virtual ~CHomeWork2Doc();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	//{{AFX_MSG(CHomeWork2Doc)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_HOMEWORK2DOC_H__49E69189_7C93_44A8_B9CC_F0397EC463DC__INCLUDED_)
