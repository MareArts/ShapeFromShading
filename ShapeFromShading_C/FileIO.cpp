// FileIO.cpp: implementation of the CFileIO class.
//
//////////////////////////////////////////////////////////////////////
#include "StdAfx.h"
#include "FileIO.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CFileIO::CFileIO()
{
	orgImg = NULL;
	outImg = NULL;
	hRGB   = NULL;
	
	height = width = 0;
}

CFileIO::~CFileIO()
{
	if(orgImg) delete []orgImg;
	if(outImg) delete []outImg;
	if(hRGB)   delete []hRGB;
}

void CFileIO::ReadBMP8(char *inbuff)
{
	printf("ISLAB(R), Made by D.J.Kang, 2005\n\n");
	FILE	*infile;
	if((infile = fopen(inbuff, "rb")) == NULL) 
	{ 
	 	printf("\nCan't open %s \n", inbuff);
		return;
	}

	// File, Image Headers
	fread(&hf,sizeof(BITMAPFILEHEADER),1,infile);
	if(hf.bfType!=0x4D42) return;

	fread(&hInfo,sizeof(BITMAPINFOHEADER),1,infile);
	if(hInfo.biBitCount !=8) 
	{ 
		printf("Bad File Format!!"); 
		return;
	}

	this->height = hInfo.biHeight;
	this->width  = hInfo.biWidth;

	// Pallette
	hRGB = new RGBQUAD [256];
	fread(hRGB,sizeof(RGBQUAD),256,infile);

	// Memory
	int rwsize = WIDTHBYTES( hInfo.biBitCount*hInfo.biWidth);
	unsigned char* lpImg = new unsigned char [hInfo.biHeight*rwsize];
	fread(lpImg,sizeof(char),hInfo.biHeight*rwsize,infile);
	fclose(infile);


	// raw Image
	orgImg = new unsigned char [hInfo.biHeight*hInfo.biWidth];
	outImg = new unsigned char [hInfo.biHeight*hInfo.biWidth];

	register int i,j;
	for(i=0; i<hInfo.biHeight; i++)
	{
		int index1 = (hInfo.biHeight-i-1)*rwsize;
		int index2 = i*hInfo.biWidth;
		for(j=0; j<hInfo.biWidth; j++)
		{
			orgImg[index2+j] = lpImg[index1+j];
		}
	}
	
	delete []lpImg; 
}

void CFileIO::WriteBMP8(char *outbuf)
{
	if(!outImg || !orgImg) return;

	FILE	*outfile;
	if((outfile = fopen(outbuf, "wb")) == NULL) 
	{ 
	 	printf("\nCan't open %s \n", outbuf);
		return;
	}

	// File, Image Headers
	fwrite(&hf,sizeof(char),sizeof(BITMAPFILEHEADER),outfile);
	fwrite(&hInfo,sizeof(char),sizeof(BITMAPINFOHEADER),outfile);

	// Pallette
	fwrite(hRGB,sizeof(RGBQUAD),256,outfile);

	// Data Inversion
	register int i,j;
	unsigned char* lpImg = new unsigned char [hInfo.biSizeImage];
	int rwsize = WIDTHBYTES( hInfo.biBitCount*hInfo.biWidth);
	for(i=0; i<hInfo.biHeight; i++)
	{
		int index1 = (hInfo.biHeight-i-1)*rwsize;
		int index2 = i*hInfo.biWidth;
		for(j=0; j<hInfo.biWidth; j++)
		{
			lpImg[index1+j] = outImg[index2+j]; 
		}
	}

	// Write data
	fwrite(lpImg,sizeof(char),hInfo.biSizeImage,outfile);

	fclose(outfile);

	delete []lpImg;
}

unsigned char* CFileIO::GetOutImg()
{
	return outImg;
}

unsigned char* CFileIO::GetOrgImg()
{
	return orgImg;
}	

int CFileIO::GetHeight()
{
	return height;
}

int CFileIO::GetWidth()
{
	return width;
}

CFileIO::CFileIO(char *inbuff)
{
	orgImg = NULL;
	outImg = NULL;
	hRGB   = NULL;
	height = width = 0;

	ReadBMP8(inbuff);
}

void CFileIO::CreatBmpHeader(int width, int height)
{
	// BITMAPFILEHEADER setting
	hf.bfOffBits = 1078;
	hf.bfReserved1 = 0;
	hf.bfReserved2 =0;
	hf.bfSize = hf.bfOffBits + width*height;
	hf.bfType = 19778;

	//BITMAPINFOHEADER setting
	hInfo.biSize = 40;
	hInfo.biWidth = width;
	hInfo.biHeight = height;
	hInfo.biPlanes = 1;
	hInfo.biBitCount = 8;
	hInfo.biCompression = 0;
	hInfo.biSizeImage = width * height;
	hInfo.biXPelsPerMeter = 0;
	hInfo.biYPelsPerMeter = 0;
	hInfo.biClrUsed = 256;
	hInfo.biClrImportant = 256;


	hRGB = new RGBQUAD [256];
	for ( int i = 0 ; i < 256; ++i)
	{
		hRGB[i].rgbBlue = hRGB[i].rgbGreen = hRGB[i].rgbRed = i;
		hRGB[i].rgbReserved = 0;
	}
	
	int rwsize = WIDTHBYTES( hInfo.biBitCount*hInfo.biWidth);
	unsigned char* lpImg = new unsigned char [hInfo.biHeight*rwsize];
	orgImg = new unsigned char [hInfo.biHeight*hInfo.biWidth];
	outImg = new unsigned char [hInfo.biHeight*hInfo.biWidth];

	memset( orgImg, 255, sizeof(char)* hInfo.biHeight * hInfo.biWidth );
	memset( outImg, 255, sizeof(char)* hInfo.biHeight * hInfo.biWidth );	
	
}

CFileIO::CFileIO(CFileIO &a)
{
	this->height = a.GetHeight();
	this->width = a.GetWidth();
//	this->orgImg = a.GetOrgImg();
//	this->outImg = a.GetOutImg();

	this->hf = a.hf;
	this->hInfo = a.hInfo;
//	this->hRGB = a.hRGB;
	
}
