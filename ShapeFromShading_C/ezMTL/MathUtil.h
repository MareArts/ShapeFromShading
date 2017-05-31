// -*- c++ -*-
////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2001 Oh-Wook Kwon, all rights reserved. ohwook@yahoo.com
//
//                          Easy Matrix Template Library
// 
// This Easy Matrix Template Library is provided "as is" without any express 
// or implied warranty of any kind with respect to this software. 
// In particular the authors shall not be liable for any direct, 
// indirect, special, incidental or consequential damages arising 
// in any way from use of the software.
// 
// Everyone is granted permission to copy, modify and redistribute this
// Easy Matrix Template Library, provided:
//  1.  All copies contain this copyright notice.
//  2.  All modified copies shall carry a notice stating who
//      made the last modification and the date of such modification.
//  3.  No charge is made for this software or works derived from it.  
//      This clause shall not be construed as constraining other software
//      distributed on the same medium as this software, nor is a
//      distribution fee considered a charge.
//
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Filename: MathUtil.h
////////////////////////////////////////////////////////////////////////////
#ifndef _MATH_UTIL_H_
#define _MATH_UTIL_H_

#include "./MatrixConfig.h"

const double LZERO=(-1.0E10);   // ~log(0)
const double minDouble=10*mtl_numeric_limits<double>::min();


inline int Mod(int x,int y) {
	return x%y;
}


inline double Log(double x){ // Safe log(x) to prevent underflow
	return (x<minDouble)?LZERO:log(x);
}


inline double Gammaln(double xx) {
	double x,tmp,ser;
	static double cof[6]={76.18009173,-86.50532033,24.01409822,
		-1.231739516,0.120858003e-2,-0.536382e-5};
	int j;
	x=xx-1.0;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.0;
	for (j=0;j<=5;j++) {
		x += 1.0;
		ser += cof[j]/x;
	}
	return -tmp+log(2.50662827465*ser);
}


inline double Gamma(double x){
	return exp(Gammaln(x));
}


inline double Digamma(double x){
	return (Gammaln(x*1.001)-Gammaln(x*0.999)) / (x*0.002);
}


template <typename T> inline T Ipow(T x, int y){
	if(y<0) { x=T(1)/x; y=-y; }
	T z=T(1);
	for(int i=1;i<=y;i++) z*=x;
	return z;
}


#endif
