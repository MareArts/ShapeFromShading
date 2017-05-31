// FileIO.h: interface for the CFileIO class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FILEIO_H__8A3F4899_7F07_4018_8A3F_237BBA09FA19__INCLUDED_)
#define AFX_FILEIO_H__8A3F4899_7F07_4018_8A3F_237BBA09FA19__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <stdio.h>
#include <MEMORY.H>
#include <windows.h>
#define WIDTHBYTES(bits) (((bits)+31)/32*4)
#define uchar unsigned char

class CFileIO  
{
public:
	CFileIO(CFileIO& a);
	CFileIO(char *inbuff);
	int GetWidth();
	int GetHeight();
	unsigned char* GetOrgImg();
	unsigned char* GetOutImg();
	void WriteBMP8(char *infile);
	void ReadBMP8(char *name);
	CFileIO();
	virtual ~CFileIO();
	//////////////////////////////////////////////////////////////////////////
	// '07 8/18 Ãß°¡ 
	void CreatBmpHeader(int width, int height);

	RGBQUAD* hRGB;
	BITMAPINFOHEADER hInfo;
	BITMAPFILEHEADER hf;
	unsigned char* orgImg;

private:
	
	unsigned char* outImg;
	int width;
	int height;
	
};

#endif // !defined(AFX_FILEIO_H__8A3F4899_7F07_4018_8A3F_237BBA09FA19__INCLUDED_)
