// -*- c++ -*-
///////////////////////////////////////////////////////////////////////////////
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
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Filename: MatrixConfig.h
///////////////////////////////////////////////////////////////////////////////
#ifndef _MATRIX_CONFIG_H_
#define _MATRIX_CONFIG_H_

///////////////////////////////////////////////////////////////////////////////
// User-configurable settings
///////////////////////////////////////////////////////////////////////////////
//
// 1. Define USE_LOCAL_RAND to use local random number generators.
#define USE_LOCAL_RAND
//
// 2. Define ALLOW_ACCESS_BY_BRACKET to allow accessing array elements by a[].
#define ALLOW_ACCESS_BY_BRACKET
//
// 3. Define number formats when matrices and vectors are printed.
#define mtl_format_int "%10d"
#define mtl_format_float "%10.5g"
#define mtl_format_double "%10.5g"
#define mtl_format_complex "%-20s"
//
// 4. Define DISABLE_COMPLEX to disable complex linear algebra.
//#define DISABLE_COMPLEX
//
// 5. Define USE_NRC_CODE to use the "Numerical Recipes in C" compatible codes.
//    Define the file name which you want to provide if you have your own.
//    If the file name is not defined, the default library is used.
//    The "Numerical Recipes in C" compatible codes are not provided 
//    in this package due to copyright.
//    I do not recommend to use the code in "Numerical Recipes in C"
//    if matrices are large (m,n > 10) or SVD of non-square matrices are needed.
//    I found that the code yields numerical errors for non-square matrices.
//
//#define USE_NRC_CODE
//#define NRC_FILE_NAME "../../ezNRC-1.0/ezNRC/nrc_Matrix.h"
//
///////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <float.h>


#pragma warning(disable: 4275 4700 4786 4800)
#include <algorithm>
#include <functional>
#include <bitset>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <valarray>
#include <vector>
#include <complex>
using namespace std;


inline bool IsComplex(const int x) { return false; }
inline bool IsComplex(const float x) { return false; }
inline bool IsComplex(const double x) { return false; }
template <typename T> inline bool IsComplex(const complex<T> x) { return true; }


inline int conj(const int x) { return x; }
inline float conj(const float x) { return x; }
inline double conj(const double x) { return x; }

inline int real(const int x) { return x; }
inline float real(const float x) { return x; }
inline double real(const double x) { return x; }

inline int imag(const int x) { return 0; }
inline float imag(const float x) { return 0; }
inline double imag(const double x) { return 0; }

inline bool mtl_iszero(const int x) { return (x==0); }
inline bool mtl_iszero(const float x) { return (x==0); }
inline bool mtl_iszero(const double x) { return (x==0); }
template <typename T> inline bool mtl_iszero(const complex<T> x) { return (real(x)==0 && imag(x)==0); }

inline bool mtl_iscomplex(const int x) { return false; }
inline bool mtl_iscomplex(const float x) { return false; }
inline bool mtl_iscomplex(const double x) { return false; }
template <typename T> inline bool mtl_iscomplex(const complex<T> x) { return true; }

inline bool mtl_isdouble(const int x) { return false; }
inline bool mtl_isdouble(const float x) { return false; }
inline bool mtl_isdouble(const double x) { return true; }
inline bool mtl_isdouble(const complex<int> x) { return false; }
inline bool mtl_isdouble(const complex<float> x) { return false; }
inline bool mtl_isdouble(const complex<double> x) { return true; }

inline int Abs(const int x) { return abs(x); }
inline float Abs(const float x) { return (float)fabs(x); }
inline double Abs(const double x) { return fabs(x); }
template <typename T> inline double Abs(const complex<T> x) {
	double re=fabs(x.real());
	double im=fabs(x.imag());
	if (re > im){ // re != 0
		double w=im/re;
		return (re*sqrt(1.0+w*w));
	}
	else if(im == 0){ // re=im=0
		return 0;
	}
	else{ // im >= re, re != 0, im != 0
		double w=re/im;
		return (im*sqrt(1.0+w*w));
	}
}
 
// If x=0, it returns 1.
// For real input, it returns +/- 1.
// For complex input, it returns a unit vector to the direction of x, x/|x|.
template <typename T> inline T mtl_sign(const T x) {
	double absx=Abs(x);
	if(absx==0) return T(1);
	else return (x/T(Abs(x)));
}


template <typename T> inline complex<T> pow(const complex<T> z, double p) {
	double mag=pow(Abs(z),p);
	double theta=atan2(imag(z),real(z));
	return complex<T>(mag*cos(p*theta),mag*sin(p*theta));
}


///////////////////////////////////////////////////////////////////////////////
// mtl_numeric_limits
//   Because the Visual C++ 6.0 has a bug in numeric_limits of <limits>
//   (numeric_limits<complex<double> >::epsilon() = 0),
//   I use my own mtl_numeric_limits template.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
class mtl_numeric_limits {
public:
	static T max();
	static T min();
	static T epsilon();
};


inline int mtl_numeric_limits<int>::max() { return INT_MAX; }
inline int mtl_numeric_limits<int>::min() { return INT_MIN; }
inline double mtl_numeric_limits<double>::max() { return DBL_MAX; }
inline double mtl_numeric_limits<double>::min() { return DBL_MIN; }
inline double mtl_numeric_limits<double>::epsilon() { return DBL_EPSILON; }
inline float mtl_numeric_limits<float>::max() { return FLT_MAX; }
inline float mtl_numeric_limits<float>::min() { return FLT_MIN; }
inline float mtl_numeric_limits<float>::epsilon() { return FLT_EPSILON; }
inline complex<double> mtl_numeric_limits<complex<double> >::max() { return complex<double>(DBL_MAX,DBL_MAX); }
inline complex<double> mtl_numeric_limits<complex<double> >::min() { return complex<double>(DBL_MIN,DBL_MIN); }
inline complex<double> mtl_numeric_limits<complex<double> >::epsilon() { return complex<double>(DBL_EPSILON,DBL_EPSILON); }
inline complex<float> mtl_numeric_limits<complex<float> >::max() { return complex<float>(FLT_MAX,FLT_MAX); }
inline complex<float> mtl_numeric_limits<complex<float> >::min() { return complex<float>(FLT_MIN,FLT_MIN); }
inline complex<float> mtl_numeric_limits<complex<float> >::epsilon() { return complex<float>(FLT_EPSILON,FLT_EPSILON); }


#pragma warning(default: 4800)
#pragma warning(disable: 4244)


#endif
