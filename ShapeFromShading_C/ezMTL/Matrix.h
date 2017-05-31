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
// Filename: Matrix.h
//
// Version: 1.0 (November 5, 2001)
//
// Programmer: Oh-Wook Kwon (ohwook@yahoo.com)
//
// Installed Operating Systems/Compilers:
//     Windows 32/Visual C++ 6.0
//     FreeBSD/GNU C++ 2.95.2
//
// Other Information:
//     Refer to README.txt in this directory.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _MATRIX_H_
#define _MATRIX_H_

// Indicate that ezMTL is used.
#ifndef __EZMTL__
#define __EZMTL__
#endif

#include "./Vector.h"

template <typename T> class Vector;
template <typename T> class Matrix;

// Template functions with default arguments
#ifdef __GNUC__
template <typename T> Vector<T> Sum(const Matrix<T>& A,int idx=1);
template <typename T> Vector<T> Mean(const Matrix<T>& A,int idx=1);
template <typename T> double Norm(const Matrix<T>& A, int p=2);
template <typename T> double NormF(const Matrix<T>& A);
template <typename T> Vector<T> Var(const Matrix<T>& A,int k=0);
template <typename T> double Cond(const Matrix<T>& A,int p=2);
template <typename T> Matrix<T> Inv(const Matrix<T>& A,int* retCode=0);
template <typename T> Matrix<T> PinvRecipe(const Matrix<T>& A,int* retCode=0);
template <typename T> Matrix<T> PinvMeschach(const Matrix<T>& A, int* retCode=0);
template <typename T> Matrix<T> Pinv(const Matrix<T>& A,int* retCode=0);
#else
template <typename T> Vector<T> Sum(const Matrix<T>& A,int idx);
template <typename T> Vector<T> Mean(const Matrix<T>& A,int idx);
template <typename T> double Norm(const Matrix<T>& A, int p);
template <typename T> double NormF(const Matrix<T>& A);
template <typename T> Vector<T> Var(const Matrix<T>& A,int k);
template <typename T> double Cond(const Matrix<T>& A,int p);
template <typename T> Matrix<T> Inv(const Matrix<T>& A,int* retCode);
template <typename T> Matrix<T> PinvRecipe(const Matrix<T>& A,int* retCode);
template <typename T> Matrix<T> PinvMeschach(const Matrix<T>& A, int* retCode);
template <typename T> Matrix<T> Pinv(const Matrix<T>& A,int* retCode);
#endif

// friend template functions
template <typename T> Matrix<T> mtl_plus(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_minus(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_plus(const Matrix<T>& A, const double b);
template <typename T> Matrix<T> mtl_minus(const Matrix<T>& A, const double b);
template <typename T> Matrix<T> mtl_times(const Matrix<T>& A, const double b);
template <typename T> Matrix<T> mtl_divide(const Matrix<T>& A, const double b);
template <typename T> Matrix<T> mtl_plus(const double b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_minus(const double b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_times(const double b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_divide(const double b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_times(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_mtimes(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Vector<T> mtl_mvtimes(const Matrix<T>& M, const Vector<T>& V);
template <typename T> Vector<T> mtl_vmtimes(const Vector<T>& V, const Matrix<T>& M);
template <typename T> Matrix<T> mtl_mrdivide(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_rdivide(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_mldivide(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_ldivide(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> mtl_power(const Matrix<T>& A, double x);
template <typename T> Matrix<T> mtl_power(const Matrix<T>& A, int x);
template <typename T> Matrix<T> mtl_mpower(const Matrix<T>& A, double x);
template <typename T> Matrix<T> mtl_mpower(const Matrix<T>& A, int x);
template <typename T> ostream& mtl_ostream(ostream& s, const Matrix<T>& A);
template <typename T> istream& mtl_istream(istream& s, Matrix<T>& A);

#ifndef DISABLE_COMPLEX
template <typename T> Matrix<T> mtl_plus(const Matrix<T>& A, const complex<T> b);
template <typename T> Matrix<T> mtl_minus(const Matrix<T>& A, const complex<T> b);
template <typename T> Matrix<T> mtl_times(const Matrix<T>& A, const complex<T> b);
template <typename T> Matrix<T> mtl_divide(const Matrix<T>& A, const complex<T> b);
template <typename T> Matrix<T> mtl_plus(const complex<T> b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_minus(const complex<T> b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_times(const complex<T> b, const Matrix<T>& A);
template <typename T> Matrix<T> mtl_divide(const complex<T> b, const Matrix<T>& A);
#endif

template <typename T> Vector<int> Size(const Matrix<T>& A);
template <typename T> int Size(const Matrix<T>& A,int d);
template <typename T> int NumRows(const Matrix<T>& M);
template <typename T> int NumCols(const Matrix<T>& M);
template <typename T> T Multiply3(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y);
template <typename T> T MultiplyXtAY(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y);
template <typename T> Matrix<T> ArrayMultiply(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> ArrayDivide(const Matrix<T>& A, const Matrix<T>& B);
template <typename T> Matrix<T> Abs(const Matrix<T>& A);
template <typename T> Matrix<T> Sign(const Matrix<T>& A);
template <typename T> Matrix<T> Pow(const Matrix<T>& A, double x);
template <typename T> Matrix<T> Pow(const Matrix<T>& A, int x);
template <typename T> Matrix<T> Powm(const Matrix<T>& A, double x);
template <typename T> Matrix<T> Powm(const Matrix<T>& A, int x);
template <typename T> Matrix<T> Exp(const Matrix<T>& A);
template <typename T> Matrix<T> Expm(const Matrix<T>& A);
template <typename T> Matrix<T> Log(const Matrix<T>& A);
template <typename T> Matrix<T> Logm(const Matrix<T>& A);
template <typename T> Matrix<T> Sqrt(const Matrix<T>& A);
template <typename T> Matrix<T> Sqrtm(const Matrix<T>& A);
template <typename T> Matrix<int> Int(const Matrix<T>& A);
template <typename T> Matrix<float> Float(const Matrix<T>& A);
template <typename T> Matrix<double> Double(const Matrix<T>& A);

template <typename T> int Eig(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int EigS(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int EigJacobi(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int EigHHolder(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int EigReal(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int EigHermitian(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int EigComplex(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
template <typename T> int SimDiag(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& E, Matrix<T>& D);
template <typename T> Matrix<T> Transpose(const Matrix<T>& A);
template <typename T> Matrix<T> ArrayTranspose(const Matrix<T>& A);
template <typename T> Matrix<T> Conj(const Matrix<T>& A);
template <typename T> Matrix<T> Cov(const Matrix<T>& A);
template <typename T> int Svd(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
template <typename T> int SvdRecipe(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
template <typename T> int SvdMeschach(const Matrix<T>& A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
template <typename T> int SvdSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> int Schur(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
template <typename T> int SchurReal(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
template <typename T> int SchurComplex(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
template <typename T> int Hess(const Matrix<T>& A, Matrix<T>& H, Matrix<T>& P);
template <typename T> Matrix<T> Orth(const Matrix<T>& A);
template <typename T> Matrix<T> Null(const Matrix<T>& A);
template <typename T> int Lu(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U, Matrix<T>& P);
template <typename T> int Lu(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U);
template <typename T> int LuSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> int Chol(const Matrix<T>& A, Matrix<T>& L);
template <typename T> int CholSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> int Ldl(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& D);
template <typename T> int LdlSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> T Det(const Matrix<T>& A);
template <typename T> double LogAbsDet(const Matrix<T>& A);
template <typename T> int Qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
template <typename T> int QrSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> int Qrcp(const Matrix<T>& A, Matrix<T>& QRCP, Vector<T>& diag, Vector<int>& pivot);
template <typename T> int QrcpSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> int Bkp(const Matrix<T>& A, Matrix<T>& BKP, Vector<int>& pivot, Vector<int>& blocks);
template <typename T> int BkpSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> double CondLU(const Matrix<T>& A);
template <typename T> double CondQR(const Matrix<T>& A,int p);
template <typename T> Vector<T> Diag(const Matrix<T>& A);
template <typename T> T Trace(const Matrix<T>& A);
template <typename T> Vector<T> Moment(const Matrix<T>& A, int order);
template <typename T> Vector<T> Std(const Matrix<T>& A);
template <typename T> Vector<T> Skewness(const Matrix<T>& A);
template <typename T> Vector<T> Kurtosis(const Matrix<T>& A);
template <typename T> Matrix<T> Reshape(const Matrix<T>& A);
template <typename T> int Solve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
template <typename T> int Rank(const Matrix<T>& A);


// ------------------------------------------------------------------------
// Matrix declaration
// ------------------------------------------------------------------------
template <typename T>
class Matrix {
private:
	// I chose to implement a matrix by a Vector of Vector because it is
	// much simpler than other methods. For sparse matrices, the method
	// wastes memory.
	int m;
	int n;
	Vector<Vector<T> > mat;
public:
	Matrix(int m_=0,int n_=0,const T& v=T());
	Matrix(int m_, int n_, const char* v_);
	Matrix(const Vector<T>& v, int m_,int n_);
	//Matrix(const Matrix<T>& A);
	virtual ~Matrix();
	//Matrix<T>&  operator=(const Matrix<T>& A);
	Matrix<T>&  operator=(const T& a);
#ifdef ALLOW_ACCESS_BY_BRACKET
	// The [] operator does not check if the arguments are in the valid range.
	// I recommend to use it only in implementation of the Matrix class.
	Vector<T>&			operator[](int i);
	const Vector<T>&	operator[](int i) const;
#endif
	Vector<T>&			operator()(int i);
	const Vector<T>&	operator()(int i) const;
	T&			operator()(int i,int j);
	const T&	operator()(int i,int j) const;
	Matrix<T>   operator+();
	Matrix<T>   operator-();
	Matrix<T>&  operator+=(const Matrix<T>& B);
	Matrix<T>&  operator-=(const Matrix<T>& B);
	Matrix<T>&  operator+=(const double b);
	Matrix<T>&  operator-=(const double b);
	Matrix<T>&  operator*=(const double b);
	Matrix<T>&  operator/=(const double b);
#ifndef DISABLE_COMPLEX
	Matrix<T>&  operator+=(const complex<T> b);
	Matrix<T>&  operator-=(const complex<T> b);
	Matrix<T>&  operator*=(const complex<T> b);
	Matrix<T>&  operator/=(const complex<T> b);
#endif
	/////////////////////////////
	// friend operator overloading
	/////////////////////////////
	// GCC-2.95-2 did not allow friend operator overloading definition outside class declaration.
	friend ostream& operator<<(ostream& s, const Matrix<T>& A) {
		return mtl_ostream(s,A);
	}
	friend ostream& operator<<(ostream& s, Matrix<T>& A) {
		return mtl_ostream(s,A);
	}
	friend istream& operator>>(istream& s, Matrix<T>& A) {
		return mtl_istream(s,A);
	}
	friend Matrix<T> operator+(const Matrix<T>& A, const Matrix<T>& B) {
		return mtl_plus(A,B);
	}
	friend Matrix<T> operator-(const Matrix<T>& A, const Matrix<T>& B) {
		return mtl_minus(A,B);
	}
	friend Matrix<T> operator+(const Matrix<T>& A, const double b) {
		return mtl_plus(A,b);
	}
	friend Matrix<T> operator-(const Matrix<T>& A, const double b) {
		return mtl_minus(A,b);
	}
	friend Matrix<T> operator*(const Matrix<T>& A, const double b) {
		return mtl_times(A,b);
	}
	friend Matrix<T> operator/(const Matrix<T>& A, const double b) {
		return mtl_divide(A,b);
	}
	friend Matrix<T> operator+(const double b, const Matrix<T>& A) {
		return mtl_plus(b,A);
	}
	friend Matrix<T> operator-(const double b, const Matrix<T>& A) {
		return mtl_minus(b,A);
	}
	friend Matrix<T> operator*(const double b, const Matrix<T>& A) {
		return mtl_times(b,A);
	}
	friend Matrix<T> operator/(const double b, const Matrix<T>& A) {
		return mtl_divide(b,A);
	}
	// Multiply Matrix * Matrix
	friend Matrix<T> operator*(const Matrix<T>& A, const Matrix<T>& B) {
		return mtl_mtimes(A,B);
	}
	// Multiply Matrix*Vector'
	friend Vector<T> operator*(const Matrix<T>& M, const Vector<T>& V) {
		return mtl_mvtimes(M,V);
	}
	// Multiply Vector*Matrix
	friend Vector<T> operator*(const Vector<T>& V, const Matrix<T>& M) {
		return mtl_vmtimes(V,M);
	}

#ifndef DISABLE_COMPLEX
	friend Matrix<T> operator+(const Matrix<T>& A, const complex<T> b) {
		return mtl_plus(A,b);
	}
	friend Matrix<T> operator-(const Matrix<T>& A, const complex<T> b) {
		return mtl_minus(A,b);
	}
	friend Matrix<T> operator*(const Matrix<T>& A, const complex<T> b) {
		return mtl_times(A,b);
	}
	friend Matrix<T> operator/(const Matrix<T>& A, const complex<T> b) {
		return mtl_divide(A,b);
	}
	friend Matrix<T> operator+(const complex<T> b, const Matrix<T>& A) {
		return mtl_plus(b,A);
	}
	friend Matrix<T> operator-(const complex<T> b, const Matrix<T>& A) {
		return mtl_minus(b,A);
	}
	friend Matrix<T> operator*(const complex<T> b, const Matrix<T>& A) {
		return mtl_times(b,A);
	}
	friend Matrix<T> operator/(const complex<T> b, const Matrix<T>& A) {
		return mtl_divide(b,A);
	}
#endif

	/////////////////////////////
	// friend functions
	/////////////////////////////
#ifdef __GNUC__
	friend Matrix<T> mtl_plus<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_minus<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_plus<T>(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_minus<T>(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_times<T>(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_divide<T>(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_plus<T>(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_minus<T>(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_times<T>(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_divide<T>(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_times<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_mtimes<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Vector<T> mtl_mvtimes<T>(const Matrix<T>& M, const Vector<T>& V);
	friend Vector<T> mtl_vmtimes<T>(const Vector<T>& V, const Matrix<T>& M);
	friend Matrix<T> mtl_mrdivide<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_rdivide<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_mldivide<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_ldivide<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_power<T>(const Matrix<T>& A, double x);
	friend Matrix<T> mtl_power<T>(const Matrix<T>& A, int x);
	friend Matrix<T> mtl_mpower<T>(const Matrix<T>& A, double x);
	friend Matrix<T> mtl_mpower<T>(const Matrix<T>& A, int x);
	friend ostream& mtl_ostream<T>(ostream& s, const Matrix<T>& A);
	friend istream& mtl_istream<T>(istream& s, Matrix<T>& A);

#ifndef DISABLE_COMPLEX
	friend Matrix<T> mtl_plus<T>(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_minus<T>(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_times<T>(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_divide<T>(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_plus<T>(const complex<T> b, const Matrix<T>& A);
	friend Matrix<T> mtl_minus<T>(const complex<T> b, const Matrix<T>& A);
	friend Matrix<T> mtl_times<T>(const complex<T> b, const Matrix<T>& A);
	friend Matrix<T> mtl_divide<T>(const complex<T> b, const Matrix<T>& A);
#endif

	friend Vector<int> Size<T>(const Matrix<T>& A);
	friend int Size<T>(const Matrix<T>& A,int d);
	friend int NumRows<T>(const Matrix<T>& M);
	friend int NumCols<T>(const Matrix<T>& M);
	friend T Multiply3<T>(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y);
	friend T MultiplyXtAY<T>(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y);
 	friend Matrix<T> ArrayMultiply<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> ArrayDivide<T>(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> Abs<T>(const Matrix<T>& A);
	friend Matrix<T> Sign<T>(const Matrix<T>& A);
	friend Matrix<T> Pow<T>(const Matrix<T>& A, double x);
	friend Matrix<T> Pow<T>(const Matrix<T>& A, int x);
	friend Matrix<T> Powm<T>(const Matrix<T>& A, double x);
	friend Matrix<T> Powm<T>(const Matrix<T>& A, int x);
	friend Matrix<T> Exp<T>(const Matrix<T>& A);
	friend Matrix<T> Expm<T>(const Matrix<T>& A);
	friend Matrix<T> Log<T>(const Matrix<T>& A);
	friend Matrix<T> Logm<T>(const Matrix<T>& A);
	friend Matrix<T> Sqrt<T>(const Matrix<T>& A);
	friend Matrix<T> Sqrtm<T>(const Matrix<T>& A);
	friend Matrix<int> Int<T>(const Matrix<T>& A);
	friend Matrix<float> Float<T>(const Matrix<T>& A);
	friend Matrix<double> Double<T>(const Matrix<T>& A);
	friend Matrix<T> Transpose<T>(const Matrix<T>& A);
	friend Matrix<T> ArrayTranspose<T>(const Matrix<T>& A);
	friend Matrix<T> Conj<T>(const Matrix<T>& A);
	friend Matrix<T> Cov<T>(const Matrix<T>& A);
	friend Matrix<T> Orth<T>(const Matrix<T>& A);
	friend Matrix<T> Null<T>(const Matrix<T>& A);

	friend int Eig<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigS<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigJacobi<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigHHolder<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigReal<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigHermitian<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigComplex<T>(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int SimDiag<T>(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& E, Matrix<T>& D);
	friend int Svd<T>(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
	friend int SvdRecipe<T>(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
	friend int SvdMeschach<T>(const Matrix<T>& A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
	friend int SvdSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Schur<T>(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
	friend int SchurReal<T>(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
	friend int SchurComplex<T>(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
	friend int Hess<T>(const Matrix<T>& A, Matrix<T>& H, Matrix<T>& P);
	friend int Lu<T>(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U, Matrix<T>& P);
	friend int Lu<T>(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U);
	friend int LuSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Chol<T>(const Matrix<T>& A, Matrix<T>& L);
	friend int CholSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Ldl<T>(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& D);
	friend int LdlSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Qr<T>(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
	friend int QrSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Qrcp<T>(const Matrix<T>& A, Matrix<T>& QRCP, Vector<T>& diag, Vector<int>& pivot);
	friend int QrcpSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Bkp<T>(const Matrix<T>& A, Matrix<T>& BKP, Vector<int>& pivot, Vector<int>& blocks);
	friend int BkpSolve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend double CondLU<T>(const Matrix<T>& A);
	friend double CondQR<T>(const Matrix<T>& A,int p);
	friend T Det<T>(const Matrix<T>& A);
	friend double LogAbsDet<T>(const Matrix<T>& A);
	friend Vector<T> Diag<T>(const Matrix<T>& A);
	friend T Trace<T>(const Matrix<T>& A);
	friend Vector<T> Moment<T>(const Matrix<T>& A, int order);
	friend Vector<T> Std<T>(const Matrix<T>& A);
	friend Vector<T> Skewness<T>(const Matrix<T>& A);
	friend Vector<T> Kurtosis<T>(const Matrix<T>& A);
	friend Matrix<T> Reshape<T>(const Matrix<T>& A);
	friend int Solve<T>(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Rank<T>(const Matrix<T>& A);

#else
	friend Matrix<T> mtl_plus(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_minus(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_plus(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_minus(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_times(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_divide(const Matrix<T>& A, const double b);
	friend Matrix<T> mtl_plus(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_minus(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_times(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_divide(const double b, const Matrix<T>& A);
	friend Matrix<T> mtl_times(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_mtimes(const Matrix<T>& A, const Matrix<T>& B);
	friend Vector<T> mtl_mvtimes(const Matrix<T>& M, const Vector<T>& V);
	friend Vector<T> mtl_vmtimes(const Vector<T>& V, const Matrix<T>& M);
	friend Matrix<T> mtl_mrdivide(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_rdivide(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_mldivide(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_ldivide(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> mtl_power(const Matrix<T>& A, double x);
	friend Matrix<T> mtl_power(const Matrix<T>& A, int x);
	friend Matrix<T> mtl_mpower(const Matrix<T>& A, double x);
	friend Matrix<T> mtl_mpower(const Matrix<T>& A, int x);
	friend ostream& mtl_ostream(ostream& s, const Matrix<T>& A);
	friend istream& mtl_istream(istream& s, Matrix<T>& A);

#ifndef DISABLE_COMPLEX
	friend Matrix<T> mtl_plus(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_minus(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_times(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_divide(const Matrix<T>& A, const complex<T> b);
	friend Matrix<T> mtl_plus(const complex<T> b, const Matrix<T>& A);
	friend Matrix<T> mtl_minus(const complex<T> b, const Matrix<T>& A);
	friend Matrix<T> mtl_times(const complex<T> b, const Matrix<T>& A);
	friend Matrix<T> mtl_divide(const complex<T> b, const Matrix<T>& A);
#endif

	friend Vector<int> Size(const Matrix<T>& A);
	friend int Size(const Matrix<T>& A,int d);
	friend int NumRows(const Matrix<T>& M);
	friend int NumCols(const Matrix<T>& M);
	friend T Multiply3(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y);
	friend T MultiplyXtAY(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y);
	friend Matrix<T> ArrayMultiply(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> ArrayDivide(const Matrix<T>& A, const Matrix<T>& B);
	friend Matrix<T> Abs(const Matrix<T>& A);
	friend Matrix<T> Sign(const Matrix<T>& A);
	friend Matrix<T> Pow(const Matrix<T>& A, double x);
	friend Matrix<T> Pow(const Matrix<T>& A, int x);
	friend Matrix<T> Powm(const Matrix<T>& A, double x);
	friend Matrix<T> Powm(const Matrix<T>& A, int x);
	friend Matrix<T> Exp(const Matrix<T>& A);
	friend Matrix<T> Expm(const Matrix<T>& A);
	friend Matrix<T> Log(const Matrix<T>& A);
	friend Matrix<T> Logm(const Matrix<T>& A);
	friend Matrix<T> Sqrt(const Matrix<T>& A);
	friend Matrix<T> Sqrtm(const Matrix<T>& A);
	friend Matrix<int> Int(const Matrix<T>& A);
	friend Matrix<float> Float(const Matrix<T>& A);
	friend Matrix<double> Double(const Matrix<T>& A);
	friend Matrix<T> Transpose(const Matrix<T>& A);
	friend Matrix<T> ArrayTranspose(const Matrix<T>& A);
	friend Matrix<T> Conj(const Matrix<T>& A);
	friend Matrix<T> Cov(const Matrix<T>& A);
	friend Matrix<T> Orth(const Matrix<T>& A);
	friend Matrix<T> Null(const Matrix<T>& A);

	friend int Eig(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigS(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigJacobi(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigHHolder(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigReal(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigHermitian(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int EigComplex(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal);
	friend int SimDiag(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& E, Matrix<T>& D);
	friend int Svd(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
	friend int SvdRecipe(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
	friend int SvdMeschach(const Matrix<T>& A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V);
	friend int SvdSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Schur(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
	friend int SchurReal(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
	friend int SchurComplex(const Matrix<T>& A, Matrix<T>& t, Matrix<T>& Q);
	friend int Lu(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U, Matrix<T>& P);
	friend int Lu(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U);
	friend int LuSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Chol(const Matrix<T>& A, Matrix<T>& L);
	friend int CholSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Ldl(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& D);
	friend int LdlSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend T Det(const Matrix<T>& A);
	friend double LogAbsDet(const Matrix<T>& A);
	friend int Qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R);
	friend int QrSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Qrcp(const Matrix<T>& A, Matrix<T>& QRCP, Vector<T>& diag, Vector<int>& pivot);
	friend int QrcpSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Bkp(const Matrix<T>& A, Matrix<T>& BKP, Vector<int>& pivot, Vector<int>& blocks);
	friend int BkpSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend double CondLU(const Matrix<T>& A);
	friend double CondQR(const Matrix<T>& A,int p);
	friend Vector<T> Diag(const Matrix<T>& A);
	friend T Trace(const Matrix<T>& A);
	friend Vector<T> Moment(const Matrix<T>& A, int order);
	friend Vector<T> Std(const Matrix<T>& A);
	friend Vector<T> Skewness(const Matrix<T>& A);
	friend Vector<T> Kurtosis(const Matrix<T>& A);
	friend Matrix<T> Reshape(const Matrix<T>& A);
	friend int Solve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x);
	friend int Rank(const Matrix<T>& A);

#endif

	/////////////////////////////
	// member functions
	/////////////////////////////
	Matrix<T> t() const;
	Matrix<T> ct() const;
	Matrix<T> i() const;
	int p(string& name) const;
	int NormRows(int i=0);
	int NormCols(int j=0);
	Vector<T> Sum_(int idx=1) const;
	int Print(ostream& s, const char *name=0);
	int Print(ostream& s, string& name);
	int Read(istream& s, string& name);
	int Read(istream& s, const char *name=0);
	int Resize(int new_m, int new_n);
	Vector<T> OutputAsRow();
	Vector<T> OutputAsCol();
	Vector<T> Row(int m1) const;
	Matrix<T> Rows(int m1,int m2) const;
	Vector<T> Col(int n1) const;
	Matrix<T> Cols(int n1,int n2) const;
	Matrix<T> SubMatrix(int m1,int m2,int n1,int n2) const;
	T Scalar(const Matrix<T>& A) const;
	int SetRow(int rowX,const Vector<T>& A);
	int SetCol(int colX,const Vector<T>& A);
	int SetIdentity();
	int Save(const string& fname);
	int SaveMatlab(const string& fname);
	int Input(const Vector<T>& V);
	int Input(const Matrix<T>& A);

	bool IsSquare() const;
	bool IsSymmetric(double tolerance=1e-6) const;
	bool IsZero(double tolerance=1e-6) const;
	bool IsEye(double tolerance=1e-6) const;
	bool IsNormal(double tolerance=1e-6) const;
	bool IsOrthogonal(double tolerance=1e-6) const;
	bool IsDiagonal(double tolerance=1e-6) const;
	bool IsPD(double tolerance=1e-6) const;
	bool IsND(double tolerance=1e-6) const;
	bool IsLowerTri(double tolerance=1e-6) const;
	bool IsUpperTri(double tolerance=1e-6) const;
	bool IsTri(double tolerance=1e-6) const;
	bool IsUnitary(double tolerance=1e-6) const;
	bool IsHermitian(double tolerance=1e-6) const;
	bool IsReal() const;
	bool IsComplex() const;
	bool IsDouble() const;

	void FromInt(const Matrix<int>& A);
	void FromFloat(const Matrix<float>& A);
	void FromDouble(const Matrix<double>& A);
	Matrix<int> ToInt() const;
	Matrix<float> ToFloat() const;
	Matrix<double> ToDouble() const;
#ifndef DISABLE_COMPLEX
	void FromCDouble(const Matrix<complex<double> >& A);
	Matrix<complex<double> > ToCDouble() const;
#endif
	Matrix<double> Real() const;
	Matrix<double> Imag() const;

public:
	/////////////////////////////
	// exported static functions
	/////////////////////////////
	static Matrix<T> Randn(int m,int n);
	static Matrix<T> Eye(int m_,int n_);
	static Matrix<T> Zeros(int m_,int n_);
	static Matrix<T> CRandn(int m,int n);

	/////////////////////////////
	// local static functions
	/////////////////////////////
	static double normnxn(const Matrix<T>& A);
	static int UpperBackSubstitute(const Matrix<T>& A,Vector<T>& b);
	static int LowerBackSubstitute(const Matrix<T>& A,Vector<T>& b);
	static int cleanNearZero(Matrix<T>& A,double tolerance=1e-6);
	static int svdClean(Vector<T>& s, double tolerance=1e-6);
	static int pinv(Matrix<T>& U, Matrix<T>& S, Matrix<T>& V, Matrix<T>& pinvA);
	static int sortEigen(Matrix<T>& EVec,Matrix<T>& EVal);
	static int sortEigenCompact(Matrix<T>& EVec,Matrix<T>& EVal);
	static int Compact2Complex(Matrix<T>& EVec, Matrix<T>& EVal, Matrix<complex<T> >& CVec, Matrix<complex<T> >& CVal);
	static int Complex2Compact(Matrix<complex<T> >& CVec, Matrix<complex<T> >& CVal, Matrix<T>& EVec, Matrix<T>& EVal);
	static double pythagoras(const T a, const T b);
	static double pythagoras3(const T a, const T b, const T c);
	static int jacobiTransform(Matrix<T>& A, Vector<T>& d, Matrix<T>& v,int *rotN);
	static void jacobiRotate(Matrix<T>& a,int i,int j,int k,int l,double s,double tau);
	static Matrix<T> ExpmReal(const Matrix<T>& A);
	static Matrix<T> ExpmComplex(const Matrix<T>& A);
	static Matrix<T> LogmReal(const Matrix<T>& A);
	static Matrix<T> LogmComplex(const Matrix<T>& A);

	/////////////////////////////
	// From Meschach library
	/////////////////////////////
	// mes_eig.h: Sensitive to precision. Use double precision.
	static T mes_in_prod(const Vector<T>& a,const Vector<T>& b);
	static T mes_in_prod_offset(const Vector<T>& a,const Vector<T>& b,int i0);
	static double mes_norm_inf(const Vector<T>& x);
	static double mes_norm2(const Vector<T>& x);
	static double mes_norm2_offset(const Vector<T>& a,int i0);	
	static void	mes_get_row(const Matrix<T>& mat,int row,Vector<T>& vec);
	static void	mes_get_col(const Matrix<T>& mat,int col,Vector<T>& vec);
	static void	mes_set_row(Matrix<T>& mat,int row,Vector<T>& vec);
	static void	mes_set_col(Matrix<T>& mat,int col,Vector<T>& vec);	
	static void	mes_copy_offset(const Vector<T>& in,Vector<T>& out,int i0);
	static int mes_trieig(Vector<T>& a,Vector<T>& b,Matrix<T>& Q,int maxIter=100);
	static int mes_symmeig(const Matrix<T>& A,Matrix<T>& Q,Vector<T>& out);
	static int mes_hhldr3cols(Matrix<T>& A, int k, int j0, double beta, double nu1, double nu2, double nu3);
	static int mes_hhldr3rows(Matrix<T>& A, int k, int i0, double beta, double nu1, double nu2, double nu3);
	static int mes_schur(Matrix<T>& A,Matrix<T>& Q,int maxIter=100);
	static int mes_schur_evals(const Matrix<T>& t,Vector<T>& real_pt,Vector<T>& imag_pt);
	static int mes_schur_vecs(const Matrix<T>& t,const Matrix<T>& Q,Matrix<T>& X_re,Matrix<T>& X_im);
	static int mes_Hfactor(Matrix<T>& A, Vector<T>& diag, Vector<T>& beta);
	static int mes_makeHQ(const Matrix<T>& H, Vector<T>& diag, Vector<T>& beta, Matrix<T>& Qout);
	static int mes_makeH(const Matrix<T>& H,Matrix<T>& Hout);
	static int mes_hhvec(Vector<T>& vec,int i0,double* beta,Vector<T>& out,T* newval);
	static int mes_hhtrvec(Vector<T>& hh,double beta,int i0,Vector<T>& in,Vector<T>& out);
	static int mes_hhtrrows(Matrix<T>& M,int i0,int j0,Vector<T>& hh,double beta);
	static int mes_hhtrcols(Matrix<T>& M,int i0,int j0,Vector<T>& hh,double beta);
	static int mes_rot_vec(Vector<T>& x,int i,int k,double c,T s,Vector<T>& out);
	static int mes_rot_rows(Matrix<T>& mat,int i,int k,double c,T s,Matrix<T>& out);
	static int mes_rot_cols(Matrix<T>& mat,int i,int k,double c,T s,Matrix<T>& out);
	static void	mes_hhldr3(double x, double y, double z, double *nu1, double *beta, double *newval);
	static void	mes_givens(T x,T y,double* c,T* s);
	// for complex matrices
	static int mes_zHfactor(Matrix<T>& A, Vector<T>& diag);
	static int mes_zHQunpack(Matrix<T>& HQ,Vector<T>& diag,Matrix<T>& Q,Matrix<T>& H);
	static int mes_zschur(Matrix<T>& A,Matrix<T>& Q,int maxIter=100);
	static int mes_zschur_evals(const Matrix<T>& t,Vector<T>& EVal);
	static int mes_zschur_vecs(const Matrix<T>& t,const Matrix<T>& Q,Matrix<T>& EVec);

	// mes_rand.h
	static double mes_mrand();
	static int mes_mrandlist(T *a, int len);
	static void mes_smrand(int seed);
	static void	mes_v_rand(Vector<T>& x);
	// mes_svd.h: Sensitive to precision. Use double precision.
	static int mes_fixsvd(Vector<T>& dd, Matrix<T>& U, Matrix<T>& V);
	static int mes_bisvd(Vector<T>& dd, Vector<T>& ff, Matrix<T>& U, Matrix<T>& V,int maxIter=100);
	static int mes_svd(const Matrix<T>& A, Matrix<T>& U, Matrix<T>& V, Vector<T>& dd);
	static int mes_bifactor(Matrix<T>& A, Matrix<T>& U, Matrix<T>& V);
	static int mes_rotsvd(Vector<T>& dd, Matrix<T>& U, Matrix<T>& V);
	// mes_solve.h
	static int mes_CHfactor(Matrix<T>& A);
	static int mes_CHsolve(const Matrix<T>& A,const Vector<T>& b,Vector<T>& x);
	static int mes_LDLfactor(Matrix<T>& A);
	static int mes_LDLsolve(const Matrix<T>& LDL,const Vector<T>& b,Vector<T>& x);
	static int mes_Usolve(const Matrix<T>& matrix,const Vector<T>& b,Vector<T>& out,double diag);
	static int mes_Lsolve(const Matrix<T>& matrix,const Vector<T>& b,Vector<T>& out,double diag);
	static int mes_UTsolve(const Matrix<T>& U,const Vector<T>& b,Vector<T>& out,double diag);
	static int mes_Dsolve(const Matrix<T>& A,const Vector<T>& b,Vector<T>& x);
	static int mes_LTsolve(const Matrix<T>& L,const Vector<T>& b,Vector<T>& out,double diag);
	// mes_factor.h
	static int mes_LUfactor(Matrix<T>& A,Vector<int>& pivot);
	static int mes_LUsolve(Matrix<T>& A,Vector<int>& pivot,const Vector<T>& b,Vector<T>& x);
	static int mes_LUTsolve(Matrix<T>& LU,Vector<int>& pivot,const Vector<T>& b,Vector<T>& x);
	static int mes_inverse(const Matrix<T>& A,Matrix<T>& out);
	static double mes_LUcondest(Matrix<T>& LU,Vector<int>& pivot);
	static int mes_QRfactor(Matrix<T>& A,Vector<T>& diag);
	static int mes_QRCPfactor(Matrix<T>& A,Vector<T>& diag,Vector<int>& px);
	static int mes_Qsolve_(Matrix<T>& QR,Vector<T>& diag,const Vector<T>& b,Vector<T>& x,Vector<T>& tmp);
	static int mes_makeQ(Matrix<T>& QR,Vector<T>& diag,Matrix<T>& Qout);
	static int mes_makeR(Matrix<T>& QR,Matrix<T>& Rout);
	static int mes_QRsolve(Matrix<T>& QR,Vector<T>& diag,const Vector<T>& b,Vector<T>& x);
	static int mes_QRCPsolve(Matrix<T>& QR,Vector<T>& diag,Vector<int>& pivot,const Vector<T>& b,Vector<T>& x);
	static int mes_Umlt(const Matrix<T>& U,const Vector<T>& x,Vector<T>& out);
	static int mes_UTmlt(const Matrix<T>& U,const Vector<T>& x,Vector<T>& out);
	static int mes_QRTsolve(const Matrix<T>& A,Vector<T>& diag,Vector<T>& c,Vector<T>& sc);
	static double mes_QRcondest(const Matrix<T>& QR,int p=2);
	static void mes_interchange(Matrix<T>& A,int i,int j);
	static int mes_BKPfactor(Matrix<T>& A,Vector<int>& pivot,Vector<int>& blocks);
	static int mes_BKPsolve(const Matrix<T>& A,Vector<int>& pivot,Vector<int>& block,const Vector<T>& b,Vector<T>& x);
	// mes_perm.h
	static int mes_px_inv(Vector<int>& px,Vector<int>& out);
	static int mes_px_mlt(Vector<int>& px1,Vector<int>& px2,Vector<int>& out);
	static int mes_px_vec(Vector<int>& px,Vector<T>& vector,Vector<T>& out);
	static int mes_pxinv_vec(Vector<int>& px,Vector<T>& x,Vector<T>& out);
	static int mes_px_transp(Vector<int>& px,int i1,int i2);
	static int mes_myqsort(int *a,int num);
	static int mes_px_sign(Vector<int>& px);
	static int mes_px_cols(Vector<int>& px,Matrix<T>& A,Matrix<T>& out);
	static int mes_px_rows(Vector<int>& px,Matrix<T>& A,Matrix<T>& out);
	// mes_math.h
	static int mes_exp(Matrix<T>& A,double eps,Matrix<T>& out,int* q_out,int* j_out);

	/////////////////////////////
	// From Numerical Recipes in C
	/////////////////////////////
	static int nrc_choldcmp(Matrix<T>& A, Vector<T>& d);
	static int nrc_cholsolv(const Matrix<T>& L, const Vector<T>& d, const Vector<T>& b, Vector<T>& x);
	static int nrc_ludcmp(Matrix<T>& A,Vector<int>& index,double* d,Vector<int>& index2);
	static int nrc_lusolv(const Matrix<T>& LU,const Vector<int>& index,const Vector<T>& b,Vector<T>& x);
	static int nrc_qrdcmp(Matrix<T>& A,Vector<T>& c,Vector<T>& d);
	static int nrc_qrsolv(const Matrix<T>& QR,const Vector<T>& c,const Vector<T>& d,const Vector<T>& b, Vector<T>& x);
	static int nrc_svdcmp(Matrix<T>& A, Vector<T>& S, Matrix<T>& V,int maxIter=100);
	static int nrc_svsolv(Matrix<T>& U, Vector<T>& S, Matrix<T>& V,const Vector<T>& b,Vector<T>& x);
	static int nrc_tred2(Matrix<T>& A,Vector<T>& d,Vector<T>& e,int compEVec);
	static int nrc_tqli(Vector<T>& d,Vector<T>& e,Matrix<T>& Z,int compEVec);

};

typedef Matrix<float> FMatrix;
typedef Matrix<double> DMatrix;
typedef Matrix<FComplex> CFMatrix;
typedef Matrix<DComplex> CDMatrix;


// ------------------------------------------------------------------------
// Matrix implementation
// ------------------------------------------------------------------------
#ifndef local_max
#define local_max(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef local_min
#define local_min(a, b)  (((a) < (b)) ? (a) : (b)) 
#endif


/////////////////////////////
// constructor and destructor
/////////////////////////////
template <typename T> Matrix<T>::Matrix(int m_, int n_, const T& v_) : m(m_), n(n_) {
	int i;
	mat.Resize(m);
	for(i=1;i<=m;i++) { mat[i].SetType(ROW_VECTOR); mat[i].Resize(n); }
	for(i=1;i<=m;i++) for(int j=1;j<=n;j++)  mat[i][j]=v_;
}


template <typename T> Matrix<T>::Matrix(int m_, int n_, const char* v_) : m(m_), n(n_) {
	int i,k;
	mat.Resize(m);
	for(i=1;i<=m;i++) { mat[i].SetType(ROW_VECTOR); mat[i].Resize(n); }
	Vector<string> num=mtl_split(v_);
	for(i=1,k=1;i<=m;i++) for(int j=1;j<=n;j++) {
		if(k<=num.Size()) { mat[i][j]=T(atof(num[k].c_str())); k++; }
		else mat[i][j]=T();
	}
}


template <typename T> Matrix<T>::Matrix(const Vector<T>& v_, int m_, int n_) : m(m_), n(n_) {
	int i,k;
	mat.Resize(m);
	for(i=1;i<=m;i++) { mat[i].SetType(ROW_VECTOR); mat[i].Resize(n); }
	for(i=1,k=1;i<=m;i++) for(int j=1;j<=n;j++) {
		if(k<=v_.Size()) { mat[i][j]=v_[k]; k++; }
		else mat[i][j]=T();
	}
}


//template <typename T> inline Matrix<T>::Matrix(const Matrix<T>& A) : mat(A.mat),m(A.m),n(A.n) {
//}


template <typename T> inline Matrix<T>::~Matrix() {
}


/////////////////////////////
// operator overloading
/////////////////////////////
//template <typename T> inline Matrix<T>& Matrix<T>::operator=(const Matrix<T>& A) {
//	if(this == &A) return *this;
//	if ( A.n != n || A.m != m) Resize(A.m, A.n);
//	mat = A.mat;
//	return *this;
//}


template <typename T> Matrix<T>& Matrix<T>::operator=(const T& a) {
	for(int i=1; i<=m; ++i) for(int j=1;j<=n; ++j) mat[i][j] = a;
	return *this;
}


#ifdef ALLOW_ACCESS_BY_BRACKET
template <typename T> inline Vector<T>& Matrix<T>::operator[](int i) { 
	return mat[i]; 
}


template <typename T> inline const Vector<T>& Matrix<T>::operator[](int i) const { 
	return mat[i]; 
}
#endif


template <typename T> inline Vector<T>& Matrix<T>::operator()(int i) {
	return mat(i); 
}


template <typename T> inline const Vector<T>& Matrix<T>::operator()(int i) const {
	return mat(i); 
}


template <typename T> inline T& Matrix<T>::operator()(int i, int j) {
	return (mat(i))(j); 
}


template <typename T> inline const T& Matrix<T>::operator()(int i, int j) const {
	return (mat(i))(j); 
}


template <typename T> inline Matrix<T> Matrix<T>::operator+() {
	return (*this);
}


template <typename T> Matrix<T> Matrix<T>::operator-() {
	Matrix<T> C=(*this);
	C*=(-1);
	return C;
}


template <typename T> Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& B) {
	assert(m == B.m && n == B.n);
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]+=B.mat[i][j];
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& B) {
	assert(m == B.m && n == B.n);
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]-=B.mat[i][j];
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator+=(const double b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]+=b;
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator-=(const double b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]-=b;
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator*=(const double b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]*=b;
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator/=(const double b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]/=b;
	return (*this);
}


#ifndef DISABLE_COMPLEX
template <typename T> Matrix<T>& Matrix<T>::operator+=(const complex<T> b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]+=b;
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator-=(const complex<T> b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]-=b;
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator*=(const complex<T> b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]*=b;
	return (*this);
}


template <typename T> Matrix<T>& Matrix<T>::operator/=(const complex<T> b) {
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++) mat[i][j]/=b;
	return (*this);
}
#endif


/////////////////////////////
// friend operator overloading
/////////////////////////////
template <typename T> inline Matrix<T> mtl_plus(const Matrix<T>& A, const Matrix<T>& B) {
	Matrix<T> C=A; C+=B; return C;
}


template <typename T> inline Matrix<T> mtl_minus(const Matrix<T>& A, const Matrix<T>& B) {
	Matrix<T> C=A; C-=B; return C;
}


template <typename T> inline Matrix<T> mtl_plus(const Matrix<T>& A, const double b) {
	Matrix<T> C=A; C+=b; return C;
}


template <typename T> inline Matrix<T> mtl_minus(const Matrix<T>& A, const double b) {
	Matrix<T> C=A; C-=b; return C;
}


template <typename T> inline Matrix<T> mtl_times(const Matrix<T>& A, const double b) {
	Matrix<T> C=A; C*=b; return C;
}


template <typename T> inline Matrix<T> mtl_divide(const Matrix<T>& A, const double b) {
	Matrix<T> C=A; C/=b; return C;
}


template <typename T> inline Matrix<T> mtl_plus(const double b, const Matrix<T>& A) {
	Matrix<T> C=A; C+=b; return C;
}


template <typename T> inline Matrix<T> mtl_minus(const double b, const Matrix<T>& A) {
	Matrix<T> C=-A; C+=b; return C;
}


template <typename T> inline Matrix<T> mtl_times(const double b, const Matrix<T>& A) {
	Matrix<T> C=A; C*=b; return C;
}


template <typename T> Matrix<T> mtl_divide(const double b, const Matrix<T>& A) {
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=C.m;i++) for(int j=1;C.n;j++)
		C.mat[i][j]=b/A.mat[i][j];
	return C;
}


// Multiply Matrix * Matrix
template <typename T> Matrix<T> mtl_mtimes(const Matrix<T>& A, const Matrix<T>& B) {
	assert(A.n == B.m);
	Matrix<T> C(A.m, B.n);
	for(int i=1;i<=A.m;i++) {
		for(int j=1;j<=B.n;j++){
			T sum=T();
			for(int k=1;k<=A.n;k++)
				sum+=A.mat[i][k]*B.mat[k][j];
			C.mat[i][j]=sum;
		}
	}
	return C;
}


#ifndef DISABLE_COMPLEX
template <typename T> inline Matrix<T> mtl_plus(const Matrix<T>& A, const complex<T> b) {
	Matrix<T> C=A; C+=b; return C;
}


template <typename T> inline Matrix<T> mtl_minus(const Matrix<T>& A, const complex<T> b) {
	Matrix<T> C=A; C-=b; return C;
}


template <typename T> inline Matrix<T> mtl_times(const Matrix<T>& A, const complex<T> b) {
	Matrix<T> C=A; C*=b; return C;
}


template <typename T> inline Matrix<T> mtl_divide(const Matrix<T>& A, const complex<T> b) {
	Matrix<T> C=A; C/=b; return C;
}


template <typename T> inline Matrix<T> mtl_plus(const complex<T> b, const Matrix<T>& A) {
	Matrix<T> C=A; C+=b; return C;
}


template <typename T> inline Matrix<T> mtl_minus(const complex<T> b, const Matrix<T>& A) {
	Matrix<T> C=-A; C+=b; return C;
}


template <typename T> inline Matrix<T> mtl_times(const complex<T> b, const Matrix<T>& A) {
	Matrix<T> C=A; C*=b; return C;
}
#endif


// Multiply Matrix*Vector'
template <typename T> Vector<T> mtl_mvtimes(const Matrix<T>& M, const Vector<T>& V) {
	assert(V.Type() == COL_VECTOR);
	assert(V.Size() == M.n); 
	Vector<T> C(M.m,COL_VECTOR);
	for(int i=1; i<=M.m; ++i) {
		T sum = 0;
		for(int k=1; k<=M.n; ++k) sum += M.mat[i][k] * V[k];
		C[i] = sum;
	}
	return C;
}


// Multiply Vector*Matrix
template <typename T> Vector<T> mtl_vmtimes(const Vector<T>& V, const Matrix<T>& M) {
	assert(V.Type() == ROW_VECTOR);
	assert(V.Size() == M.m); 
	Vector<T> C(M.n,ROW_VECTOR);
	for(int j=1; j<=M.n; ++j) {
		T sum = 0;
		for(int k=1; k<=M.m; ++k) sum += V[k] * M.mat[k][j];
		C[j] = sum;
	}
	return C;
}


// Array Multiply
template <typename T> Matrix<T>  mtl_times(const Matrix<T>& A, const Matrix<T>& B){
	assert(A.m==B.m && A.n==B.n);
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C.mat[i][j]=A.mat[i][j]*B.mat[i][j];
	return C;
}


// Array right Division A(i,j)/B(i,j)
template <typename T> Matrix<T>  mtl_rdivide(const Matrix<T>& A, const Matrix<T>& B){
	assert(A.m==B.m && A.n==B.n);
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C.mat[i][j]=A.mat[i][j]/B.mat[i][j];
	return C;
}


// Array left Division B(i,j)/A(i,j)
template <typename T> Matrix<T>  mtl_ldivide(const Matrix<T>& A, const Matrix<T>& B){
	assert(A.m==B.m && A.n==B.n);
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C.mat[i][j]=B.mat[i][j]/A.mat[i][j];
	return C;
}


// Matrix right Division A\B = pinv(A)*b
template <typename T> Matrix<T>  mtl_mldivide(const Matrix<T>& A, const Matrix<T>& B){
	Matrix<T> C(A.n,B.n);
	Vector<T> b(B.m),x(A.n);
	for(int j=1;j<=B.n;j++){
		b=B.Col(j);
		int isOK=Solve(A,b,x);
		C.SetCol(j,x);
	}
	return C;
}


// Matrix left Division A/B = A*inv(B). More precisely, A/B = (B'\A')'
template <typename T> Matrix<T>  mtl_mrdivide(const Matrix<T>& A, const Matrix<T>& B){
	return Transpose(mtl_mldivide(Transpose(B),Transpose(A)));
}


template <typename T> Matrix<T> mtl_power(const Matrix<T>& A, double x){
	Matrix<T> C(A.m,A.n);
	if(x==0) return Matrix<T>::Eye(A.m,A.n);
	else{
		if(Abs(x-(int)x)<real(mtl_numeric_limits<T>::epsilon())){
			return mtl_power(A,(int)x);
		}
		for(int i=1;i<=A.m;i++){
			for(int j=1;j<=A.n;j++) {
				if(A.IsReal() && real(A[i][j])<0){
					cerr << "ERROR: The A matrix has negative elements.\n";
					cerr << "       Use a complex matrix.\n";
					return Matrix<T>();
				}
				C[i][j]=pow(A[i][j],x);
			}	
		}
		return C;
	}
}


template <typename T> Matrix<T> mtl_power(const Matrix<T>& A, int x){
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C[i][j]=Ipow(A[i][j],x);
	return C;
}


template <typename T> Matrix<T> mtl_mpower(const Matrix<T>& A, double x){
	assert(A.m==A.n);
	if(x==0){
		return Matrix<T>::Eye(A.m,A.n);
	}
	else if(x>0){
		if(Abs(x-(int)x)<real(mtl_numeric_limits<T>::epsilon())){
			return mtl_mpower(A,(int)x);
		}
		Matrix<T> EVec, EVal;
		int isOK=Eig(A,EVec,EVal);
		Matrix<T> C(A.m,A.n); C=T();
		for(int i=1;i<=A.m;i++){
			if(A.IsReal() && real(EVal[i][i])<0){
				cerr << "ERROR: The A matrix has negative eigenvalues.\n";
				cerr << "       Use a complex matrix.\n";
				return Matrix<T>();
			}
			C[i][i]=pow(EVal[i][i],x);
		}
		C=EVec*C*Inv(EVec);
		return C;
	}
	else{
		return Powm(Inv(A),-x);
	}
}


template <typename T> Matrix<T> mtl_mpower(const Matrix<T>& A, int x){
	assert(A.m==A.n);
	if(x==0){
		return Matrix<T>::Eye(A.m,A.n);
	}
	else if(x>0){
		Matrix<T> C(A.m,A.n);
		C.SetIdentity();
		for(int i=1;i<=x;i++)
			C=C*A;
		return C;
	}
	else{
		int retCode;
		Matrix<T> C=Inv(A,&retCode);
		if(retCode==0) {
			cerr << "ERROR: Singular matrix.\n";
			assert(0);
			Matrix<T> Inf(A.m,A.n,mtl_numeric_limits<T>::max());
			return Inf;
		}
		return Powm(C,-x);
	}
}


template <typename T> ostream& mtl_ostream(ostream& s, const Matrix<T>& A) {
	for(int i=1; i<=A.m; i++) {
		s << "        ";
		for(int j=1; j<=A.n; j++) {
			s << Num2str(A[i][j]) << " ";
		}
		s << endl;
	}
	return s;
}


template <typename T> istream& mtl_istream(istream& s, Matrix<T>& A) {
	for (int i=1; i<=A.m; i++){
		for (int j=1; j<=A.n; j++)
			s >> A[i][j];
	}
	return s;
}



/////////////////////////////
// member functions
/////////////////////////////
template <typename T> int Matrix<T>::Resize(int new_m, int new_n) {
	if(m==new_m && n==new_n) return 1;
	// must preserve data
	int i;
	for(i=1;i<=local_min(m,new_m);i++) { mat[i].Resize(new_n); mat[i].SetType(ROW_VECTOR); }
	mat.Resize(new_m);
	for(i=m+1;i<=new_m;i++) { mat[i].Resize(new_n); mat[i].SetType(ROW_VECTOR); }
	m=new_m; n=new_n;
	return 1;
}


template <typename T> Vector<T> Matrix<T>::Sum_(int idx) const{
	assert(idx==1 || idx==2);
	Vector<T> C;
	int i,j;
	if(idx==1){
		C.SetType(ROW_VECTOR);
		C.Resize(n);
		for(j=1;j<=n;j++){
			T sum=0; for(i=1;i<=m;i++) sum+=mat[i][j]; C[j]=sum;
		}
	}
	else if(idx==2){
		C.SetType(COL_VECTOR);
		C.Resize(m);
		for(i=1;i<=m;i++){
			T sum=0; for(j=1;j<=n;j++) sum+=mat[i][j]; C[i]=sum;
		}
	}
	return C;
}


template <typename T> Vector<T> Matrix<T>::OutputAsRow(){
	Vector<T> C(m*n,ROW_VECTOR);
	int k=1;
	for(int i=1;i<=m;i++) for(int j=1;j<=n;j++)
		C[k++]=mat[i][j];
	return C;
}


template <typename T> Vector<T> Matrix<T>::OutputAsCol(){ 
	Vector<T> C(m*n,COL_VECTOR);
	int k=1;
	for(int j=1;j<=n;j++) for(int i=1;i<=m;i++) 
		C[k++]=mat[i][j];
	return C;
}


template <typename T> Vector<T> Matrix<T>::Row(int rowX) const {
	assert(rowX>=1 && rowX<=m);
	return mat[rowX];
}


template <typename T> Matrix<T> Matrix<T>::Rows(int m1,int m2) const {
	assert(m1>=1 && m1<=m && m2>=1 && m2<=m && m2>=m1);
	Matrix<T> C(m2-m1+1,n);
	for(int i=1;i<=m2-m1+1;i++) for(int j=1;j<=n;j++) 
		C[i][j]=mat[i+m1-1][j];
	return C;
}


template <typename T> Vector<T> Matrix<T>::Col(int colX) const {
	assert(colX>=1 && colX<=n);
	Vector<T> C(m,COL_VECTOR);
	for(int i=1;i<=m;i++) C[i]=mat[i][colX];
	return C;
}


template <typename T> Matrix<T> Matrix<T>::Cols(int n1,int n2) const {
	assert(n1>=1 && n1<=n && n2>=1 && n2<=n && n2>=n1);
	Matrix<T> C(m,n2-n1+1);
	for(int j=1;j<=n2-n1+1;j++) for(int i=1;i<=m;i++)
		C[i][j]=mat[i][j+n1-1];
	return C;
}


template <typename T> int Matrix<T>::SetRow(int rowX,const Vector<T>& A){
	assert(n==A.Size());
	mat(rowX)=A;
	return 1;
}


template <typename T> int Matrix<T>::SetCol(int colX,const Vector<T>& A){
	assert(m==A.Size());
	for(int i=1;i<=m;i++) mat[i][colX]=A[i];
	return 1;
}


template <typename T> 	int Matrix<T>::Print(ostream& s, const char *name) {
	if(name != 0 && name[0] != 0) s << name << " ";
	s << 2 << " " << m << " " << n << endl;
	s << (*this) << endl;
	return 1;
}


template <typename T> 	int Matrix<T>::Print(ostream& s, string& name) {
	return Print(s,name.c_str());
}


template <typename T> 	int Matrix<T>::Read(istream& s, const char *name) {
	string dummy;
	s >> dummy;
	if(name != 0 && name[0] != 0){
		assert(name==dummy);
	}
	int itmp, mm, nn;
	s >> itmp >> mm >> nn;
	assert(itmp==2);
	Resize(mm,nn);
	s >> (*this);
	return 1;
}


template <typename T> 	int Matrix<T>::Read(istream& s, string& name) {
	return Read(s,name.c_str());
}


template <typename T> int Matrix<T>::SetIdentity() {
	(*this)=T();
	for (int i=1; i<=local_min(m,n); i++) mat[i][i] = 1;
	return 1;
}


template <typename T> Matrix<T> Matrix<T>::SubMatrix(int m1,int m2, int n1, int n2) const{
	Matrix<T> C(m2-m1+1,n2-n1+1);
	for(int i=1;i<=m2-m1+1;i++){
		for(int j=1;j<=n2-n1+1;j++){
			C.mat[i][j]=mat[m1-1+i][n1-1+j];
		}
	}
	return C;
}


template <typename T> int Matrix<T>::NormCols(int j)
{
	if(j==0) for(int k=1;k<=n;k++) NormCols(k);	
	else{
		int i;
		double sqSum=0;
		for(i=1;i<=m;i++) sqSum+=real(mat[i][j]*conj(mat[i][j]));
		if(sqSum>0){
			double normWj=sqrt(sqSum);
			for(i=1;i<=m;i++) mat[i][j]/=normWj;
		}
	}
	return 1;
}


template <typename T> int Matrix<T>::NormRows(int i)
{
	if(i==0) for(int k=1;k<=n;k++) NormRows(k);	
	else{
		int j;
		double sqSum=0;
		for(j=1;j<=n;j++) sqSum+=real(mat[i][j]*conj(mat[i][j]));
		if(sqSum>0){
			double normWi=sqrt(sqSum);
			for(j=1;j<=n;j++) mat[i][j]/=normWi;
		}
	}
}


template <typename T> int Matrix<T>::Save(const string& fname){
	FILE* fp=fopen(fname.c_str(),"wb");
	if(!fp) return 0;
	fprintf(fp,"FMAT %d %d\n",m,n);
	for(int i=1;i<=m;i++){
		for(int j=1;j<=n;j++)
			fprintf(fp,"%f ",mat[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return m*n;
}


template <typename T> int Matrix<T>::SaveMatlab(const string& fname){
	FILE* fp=fopen(fname.c_str(),"wb");
	if(!fp) return 0;
	for(int i=1;i<=m;i++){
		for(int j=1;j<=n;j++)
			fprintf(fp,"%f ",mat[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
	return m*n;
}


template <typename T> int Matrix<T>::Input(const Vector<T>& V){
	int k=1;
	for(int i=1;i<=m;i++){
		for(int j=1;j<=n;j++)
			mat[i][j]=V[k++];
	}
	return m*n;
}


template <typename T> int Matrix<T>::Input(const Matrix<T>& A){
	int mm=local_min(m,A.m);
	int nn=local_min(n,A.n);
	for(int i=1;i<=mm;i++){
		for(int j=1;j<=nn;j++)
			mat[i][j]=A.mat[i][j];
	}
	return mm*nn;
}



/////////////////////////////
// basic friend functions
/////////////////////////////
template <typename T> inline Vector<int> Size(const Matrix<T>& A) {
	Vector<int> C(2,ROW_VECTOR);
	C(1)=A.m; C(2)=A.n;
	return C;
}


template <typename T> inline int Size(const Matrix<T>& A,int d) {
	assert(d==1 || d==2);
	return ((d==1) ? A.m : A.n);
}


template <typename T> inline int NumRows(const Matrix<T>& M) {
	return M.m;
}


template <typename T> inline int NumCols(const Matrix<T>& M) {
	return M.n;
}


template <typename T> inline double Scalar(const Matrix<T>& A) {
	assert(A.m==1 && A.n==1);
	return A[1][1];
}


// To remove memory allocation
template <typename T> T Multiply3(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y){
	assert(x.Type() == ROW_VECTOR);
	assert(y.Type() == COL_VECTOR);
	assert(x.Size()==B.m && B.n==y.Size());
	T sum=T();
	for(int j=1;j<=B.n;j++){
		T s=T();
		for(int i=1;i<=B.m;i++) s+=x[i]*B.mat[i][j];
		sum += s*y[j];
	}
	return sum;
}


// To remove memory allocation
template <typename T> T MultiplyXtAY(const Vector<T>& x, const Matrix<T>& B, const Vector<T>& y){
	assert(x.Type() == COL_VECTOR);
	assert(y.Type() == COL_VECTOR);
	assert(x.Size()==B.m && B.n==y.Size());
	T sum=T();
	for(int j=1;j<=B.n;j++){
		T s=T();
		for(int i=1;i<=B.m;i++) s+=x[i]*B.mat[i][j];
		sum += s*y[j];
	}
	return sum;
}


// Transpose: return transposed matrix of given matrix.
template <typename T> Matrix<T> Matrix<T>::t() const {
	int i,j;
	Matrix<T> B(n, m);
	for (i=1; i<=m; i++) for (j=1; j<=n; j++){
			B.mat[j][i]=(mat[i][j]);
	}
	return B;
}


// Conjugate Transpose: return conjugate transposed matrix of given matrix.
template <typename T> Matrix<T> Matrix<T>::ct() const {
	int i,j;
	Matrix<T> B(n, m);
	for (i=1; i<=m; i++) for (j=1; j<=n; j++){
			B.mat[j][i]=conj(mat[i][j]);
	}
	return B;
}


template <typename T> inline Matrix<T> Matrix<T>::i() const {
	return Inv(*this);
}


template <typename T> inline int Matrix<T>::p(string& name) const {
	return (*this).Print(cout,name);
}


template <typename T> inline Matrix<T> Transpose(const Matrix<T>& A) {
	return A.ct();
}


template <typename T> inline Matrix<T> ArrayTranspose(const Matrix<T>& A) {
	return A.ct();
}


// Conjugate: return conjugate transposed matrix of given matrix.
template <typename T> Matrix<T> Conj(const Matrix<T>& A) {
	int i,j;
	Matrix<T> B(A.m, A.n);
	for (i=1; i<=m; i++) for (j=1; j<=n; j++){
			B.mat[i][j]=conj(mat[i][j]);
	}
	return B;
}


// return the real matrix of a complex matrix.
template <typename T> Matrix<double> Matrix<T>::Real() const {
	int i,j;
	Matrix<double> B(m, n);
	for (i=1; i<=m; i++) for (j=1; j<=n; j++){
			B[i][j]=real((*this)[i][j]);
	}
	return B;
}


// return the imaginary matrix of a complex matrix.
template <typename T> Matrix<double> Matrix<T>::Imag() const {
	int i,j;
	Matrix<double> B(m, n);
	for (i=1; i<=m; i++) for (j=1; j<=n; j++){
			B[i][j]=imag((*this)[i][j]);
	}
	return B;
}


template <typename T> Vector<T> Diag(const Matrix<T>& A){
	Vector<T> C(A.m,COL_VECTOR);
	for(int i=1;i<=A.m;i++) C[i]=A.mat[i][i];
	return C;
}


template <typename T> inline Vector<T> Sum(const Matrix<T>& A,int idx=1){
	return A.Sum_(idx);
}


template <typename T> Vector<T> Mean(const Matrix<T>& A,int idx=1){
	assert(idx==1 || idx==2);
	int sampleN=(idx==1)?NumRows(A):NumCols(A);
	return Sum(A,idx)/T(sampleN);
}


template <typename T> Matrix<T> Cov(const Matrix<T>& A){
	Vector<T> mean=Mean(A);
	return (Transpose(A)*A)/(double)(A.m)-OuterProd(mean,mean);
}


// Trace
template <typename T> T Trace(const Matrix<T>& A){
	assert(A.m==A.n);
	T sum=T();
	for(int i=1;i<=A.m;i++) sum+=A[i][i];
	return T(sum);
}


// Array Multiply: A(i,j)*B(i,j)
template <typename T> Matrix<T>  ArrayMultiply(const Matrix<T>& A, const Matrix<T>& B){
	return mtl_times(A,B);
}


// Array Division: A(i,j)/B(i,j)
template <typename T> Matrix<T>  ArrayDivide(const Matrix<T>& A, const Matrix<T>& B){
	return mtl_rdivide(A,B);
}


template <typename T> Matrix<T> Abs(const Matrix<T>& A){
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C[i][j]=T(Abs(A[i][j]));
	return C;
}


template <typename T> Matrix<T> Sign(const Matrix<T>& A, double x){
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C[i][j]=T( ( (A[i][j]>0) ? 1 : ((A[i][j]<0) ? (-1) : 0) ) );
	return C;
}


template <typename T> Matrix<T> Pow(const Matrix<T>& A, double x){
	return mtl_power(A,x);
}


template <typename T> Matrix<T> Pow(const Matrix<T>& A, int x){
	return mtl_power(A,x);
}


template <typename T> Matrix<T> Powm(const Matrix<T>& A, double x){
	return mtl_mpower(A,x);
}


template <typename T> Matrix<T> Powm(const Matrix<T>& A, int x){
	return mtl_mpower(A,x);
}


template <typename T> Matrix<T> Exp(const Matrix<T>& A){
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C[i][j]=exp(A[i][j]);
	return C;
}


template <typename T> Matrix<T> Log(const Matrix<T>& A){
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C[i][j]=log(A[i][j]);
	return C;
}


template <typename T> Matrix<T> Sqrt(const Matrix<T>& A){
	Matrix<T> C(A.m,A.n);
	for(int i=1;i<=A.m;i++)
		for(int j=1;j<=A.n;j++) C[i][j]=sqrt(A[i][j]);
	return C;
}


template <typename T> void Matrix<T>::FromInt(const Matrix<int>& A){
	Resize(NumRows(A),NumCols(A));
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) (*this)[i][j]=T(A[i][j]);
}


template <typename T> void Matrix<T>::FromFloat(const Matrix<float>& A){
	Resize(NumRows(A),NumCols(A));
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) (*this)[i][j]=T(A[i][j]);
}


template <typename T> void Matrix<T>::FromDouble(const Matrix<double>& A){
	Resize(NumRows(A),NumCols(A));
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) (*this)[i][j]=T(A[i][j]);
}


#ifndef DISABLE_COMPLEX
inline void Matrix<complex<float> >::FromCDouble(const Matrix<complex<double> >& A){
	Resize(NumRows(A),NumCols(A));
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) (*this)[i][j]=complex<float>(real(A[i][j]),imag(A[i][j]));
}
#endif


template <typename T> Matrix<int> Matrix<T>::ToInt() const {
	Matrix<int> C(m,n);
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) C[i][j]=(int)(*this)[i][j];
	return C;
}


template <typename T> Matrix<float> Matrix<T>::ToFloat() const {
	Matrix<float> C(m,n);
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) C[i][j]=(float)(*this)[i][j];
	return C;
}


template <typename T> Matrix<double> Matrix<T>::ToDouble() const {
	Matrix<double> C(m,n);
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) C[i][j]=(double)(*this)[i][j];
	return C;
}


#ifndef DISABLE_COMPLEX
template <typename T> Matrix<complex<double> > Matrix<T>::ToCDouble() const {
	Matrix<complex<double> > C(m,n);
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++) C[i][j]=complex<double>(real((*this)[i][j]), imag((*this)[i][j]));
	return C;
}
#endif


template <typename T> Matrix<int> Int(const Matrix<T>& A) {
	Matrix<int> C(NumRows(A),NumCols(A));
	for(int i=1;i<=NumRows(A);i++)
		for(int j=1;j<=NumCols(A);j++) C[i][j]=(int)A[i][j];
	return C;
}


template <typename T> Matrix<float> Float(const Matrix<T>& A) {
	Matrix<float> C(NumRows(A),NumCols(A));
	for(int i=1;i<=NumRows(A);i++)
		for(int j=1;j<=NumCols(A);j++) C[i][j]=(float)A[i][j];
	return C;
}


template <typename T> Matrix<double> Double(const Matrix<T>& A) {
	Matrix<double> C(NumRows(A),NumCols(A));
	for(int i=1;i<=NumRows(A);i++)
		for(int j=1;j<=NumCols(A);j++) C[i][j]=(double)A[i][j];
	return C;
}


template <typename T> Vector<T> Moment(const Matrix<T>& A, int order){
	Matrix<T> C=A;
	Vector<T> mean=Mean(C,1);
	for(int i=1;i<=C.m;i++){
		C[i]-=mean;
		C[i]=Pow(C[i],order);
	}
	Vector<T> mom=Mean(C,1);
	return mom;
}


template <typename T> Vector<T> Var(const Matrix<T>& A,int k=0){
	Matrix<T> C=A;
	Vector<T> mean=Mean(C,1);
	for(int i=1;i<=NumRows(C);i++){
		C[i]-=mean;
		C[i]=Pow(C[i],2);
	}
	if(k==1) return Sum(C,1)/T(NumRows(C));
	else return Sum(C,1)/T(NumRows(C)-1);
}


template <typename T> inline Vector<T> Std(const Matrix<T>& A){
	return Pow(Var(A),0.5);
}


template <typename T> Vector<T> Skewness(const Matrix<T>& A){
	Vector<T> mom3=Moment(A,3);
	Vector<T> std3=Pow(Moment(A,2),1.5);
	return ArrayDivide(mom3,std3);
}


template <typename T> Vector<T> Kurtosis(const Matrix<T>& A){
	Vector<T> mom4=Moment(A,4);
	Vector<T> std4=Pow(Moment(A,2),2);
	return ArrayDivide(mom4,std4);
}


template <typename T> Matrix<T> Reshape(const Matrix<T>& A, int new_m, int new_n){
	assert(new_m*new_n == NumRows(A)*NumCols(A));
	Matrix<T> C(new_m,new_n);
	int i,j;
	int rowN=NumRows(A);
	for(j=1;j<=new_n;j++){
		for(i=1;i<=new_m;i++){
			int k=(j-1)*new_m+i-1;
			int colX=k/rowN+1;
			int rowX=k%rowN+1;
			C[i][j]=A[rowX][colX];
		}
	}
	return C;
}


template <typename T> double Norm(const Matrix<T>& A, int p=2) {
	int m=NumRows(A);
	int n=NumCols(A);
	int i,j;
	if(p==1){
		double maxsumabs=0;
		for(j=1;j<=n;j++){
			double sumabs=0;
			for(i=1;i<=m;i++) sumabs+=Abs(A[i][j]);
			if(maxsumabs<sumabs) maxsumabs=sumabs;
		}
		return maxsumabs;
	}
	else if(p==2){
		if(m==n){ return Matrix<T>::normnxn(A); }
		int l=local_min(m,n);
		Matrix<T> B(l,l);
		Matrix<T> trA=Transpose(A);
		if(m<n)
			B=A*trA;
		else
			B=trA*A;
		return sqrt(Matrix<T>::normnxn(B));
	}
	else if(p==mtl_numeric_limits<int>::max()){
		double maxsumabs=0;
		for(i=1;i<=m;i++) {
			double sumabs=0;
			for(j=1;j<=n;j++) sumabs+=Abs(A[i][j]);
			if(maxsumabs<sumabs) maxsumabs=sumabs;
		}
		return maxsumabs;
	}
	else{
		cerr << "ERROR: Invalid p-norm\n";
		assert(0);
		return 0;
	}
}


template <typename T> double NormF(const Matrix<T>& A) {
	return sqrt(real(Sum(Diag(Transpose(A)*A))));
}


template <typename T> double Matrix<T>::normnxn(const Matrix<T>& A) {
	int m=A.m; int n=A.n;
	Matrix<T> U(m,m), S(m,n), V(n,n);
	int a1=Svd(A,U,S,V);
	if(a1==0) {
		cerr << "Svd failed\n";
		assert(0);
		return 0;
	}
	double maxW=0;
	for(int i=1;i<=local_min(m,n);i++){
		if(maxW<Abs(S[i][i])) maxW=Abs(S[i][i]);
	}
	return maxW;
}


// m==n
template <typename T> inline bool Matrix<T>::IsSquare() const {
	return (m==n);
}


// A(i,j)=conj(A(j,i))
template <typename T> bool Matrix<T>::IsSymmetric(double tolerance) const {
	int i,j;
	double err;
	T c1,c2;
	if(m!=n) return false;	
	for (i=1; i<=m-1; i++) {
		for (j=i+1; j<=m; j++) {
			c1 = mat[i][j];
			c2 = conj(mat[j][i]);
			if (!mtl_iszero(c1+c2)) {
				err = Abs(T(0.5)*(c1-c2)/(c1+c2));
				if (err>tolerance) {
					return false;
				} 
			}
		}
	}
	return true;
}


// A(i,j)=conj(A(j,i)) for complex
template <typename T> bool Matrix<T>::IsHermitian(double tolerance) const {
	return IsSymmetric(tolerance);
}


// A(i,i)=1, A(i,j)=0, i neq j
template <typename T> bool Matrix<T>::IsEye(double tolerance) const {
	int i,j;
	if(m!=n) return false;
	for (int i=1; i<=m; i++) {
		for (int j=1; j<=n; j++) {
			if(i==j) { if(Abs(mat[i][i]-1)>tolerance) return false; }
			else if(Abs(mat[i][j])>tolerance) return false;
		}
	}
	return true;
}


// A(i,j)=0
template <typename T> bool Matrix<T>::IsZero(double tolerance) const {
	for (int i=1; i<=m; i++){
		for (int j=1; j<=n; j++){
			if(Abs(mat[i][j])>tolerance) return false;
		}
	}
	return true;
}


// A(i,j)=0, i neq j
template <typename T> bool Matrix<T>::IsDiagonal(double tolerance) const {
	if(m!=n) return false;
	for (int i=1; i<=m; i++) for (int j=1; j<=n; j++) {
		if(i!=j && Abs(mat[i][j])>tolerance) return false;
	}
	return true;
}


// A(i,j)=0, i<j
template <typename T> bool Matrix<T>::IsLowerTri(double tolerance) const {
	if(m!=n) return false;
	for (int i=1; i<=m; i++) for (int j=i+1; j<=n; j++) {
		if(Abs(mat[i][j])>tolerance) return false;
	}
	return true;
}


// A(i,j)=0, i>j
template <typename T> bool Matrix<T>::IsUpperTri(double tolerance) const {
	if(m!=n) return false;
	for (int i=1; i<=m; i++) for (int j=1; j<i; j++) {
		if(Abs(mat[i][j])>tolerance) return false;
	}
	return true;
}


// IsLowerTri or IsUpperTri
template <typename T> bool Matrix<T>::IsTri(double tolerance) const {
	return (IsUpperTri(tolerance) || IsLowerTri(tolerance));
}


// A'*A=I
template <typename T> bool Matrix<T>::IsOrthogonal(double tolerance) const {
	if(m!=n) return false;
	Matrix<T> C=Transpose(*this)*(*this);
	return C.IsEye();
}


// A'*A=I for complex
template <typename T> bool Matrix<T>::IsUnitary(double tolerance) const {
	return C.IsOrthogonal(tolerance);
}


// real matrix
template <typename T> bool Matrix<T>::IsReal() const {
	T x;
	return (mtl_iscomplex(x)==false);
}


// complex matrix
template <typename T> bool Matrix<T>::IsComplex() const {
	T x;
	return mtl_iscomplex(x);
}


// complex matrix
template <typename T> bool Matrix<T>::IsDouble() const {
	T x;
	return mtl_isdouble(x);
}


// A'*A=A*A'
template <typename T> bool Matrix<T>::IsNormal(double tolerance) const {
	if(m!=n) return false;
	Matrix<T> C=Transpose(*this)*(*this)-(*this)*Transpose(*this);
	return C.IsZero();
}


// Positive definite
template <typename T> bool Matrix<T>::IsPD(double tolerance) const {
	if(IsSymmetric()==false) return false;
	Matrix<T> L;
	int isOK=Chol(*this,L);
	return ((isOK!=0) ? true : false);
}


// Negative definite
template <typename T> bool Matrix<T>::IsND(double tolerance) const {
	Matrix<T> C=-(*this);
	return C.IsPD();
}


// cleaning up near-zero elements to prevent underflow in SVD.
template <typename T> int Matrix<T>::cleanNearZero(Matrix<T>& A,double tolerance) {
	int  i,j;
	int m=NumRows(A); int n=NumCols(A);
	int nearzero=0;
	double d=0;
	for (i=1; i<=m; i++) for(j=1;j<=n;j++) if(Abs(A[i][j])>d) d=Abs(A[i][j]);
	if(d==0) return m*n;
	for (i=1; i<=m; i++){
		for (j=1; j<=n; j++){
			if (Abs(A[i][j]/T(d))<tolerance) { A[i][j] = 0; nearzero++; }
		}
	}
	return nearzero;
}


// The input is matrix 'EVal' with eigenvalues as diagonal elements
// and matrix 'EVec' with the eigenvectors as columns.
// The eigenvectors are sorted in the decreasing order of absolute eigenvalues.  
template <typename T> int Matrix<T>::sortEigen(Matrix<T>& EVec, Matrix<T>& EVal) {
	if(EVal.IsReal()==true && EVal.IsDiagonal()==false){
		return sortEigenCompact(EVec,EVal);
	}
	int m=NumRows(EVal);
	int n=NumCols(EVal);
	int j,col;
	for (col=1; col<=n-1; col++){
		double pmax=-1;
		int kmax=col;
		for (j=col+1; j<=n; j++){
			double p=Abs(EVal[j][j]);
			if (p > pmax) { pmax=p; kmax=j; }
		}
		if(kmax != col){
			T tmp = EVal[kmax][kmax];
			EVal[kmax][kmax] = EVal[col][col];
			EVal[col][col] = tmp;
			for(int i=1;i<=m;i++){
				tmp=EVec[i][col];
				EVec[i][col]=EVec[i][kmax];
				EVec[i][kmax]=tmp;
			}
		}
	}
	return 1;
}


// The input is matrix 'EVal' with eigenvalues as diagonal or 2x2 diagonal elements
// Eval and EVec are of real types.
template <typename T> int Matrix<T>::sortEigenCompact(Matrix<T>& EVec, Matrix<T>& EVal) {
	Matrix<complex<T> > V,D;
	Compact2Complex(EVec,EVal,V,D);
	Matrix<complex<T> >::sortEigen(V,D);
	Complex2Compact(V,D,EVec,EVal);
	return 1;
}


inline int Matrix<complex<float> >::sortEigenCompact(Matrix<complex<float> >& EVec, Matrix<complex<float> >& EVal) {
	cerr << "ERROR: Must not arrive here.\n";
	assert(0);
	return 0;
}


inline int Matrix<complex<double> >::sortEigenCompact(Matrix<complex<double> >& EVec, Matrix<complex<double> >& EVal) {
	cerr << "ERROR: Must not arrive here.\n";
	assert(0);
	return 0;
}


template <typename T> int Matrix<T>::Compact2Complex(Matrix<T>& EVec, Matrix<T>& EVal, Matrix<complex<T> >& CVec, Matrix<complex<T> >& CVal) {
	int m=NumRows(EVal);
	int n=NumCols(EVal);
	CVal.Resize(m,n); CVal=complex<T>();
	CVec.Resize(m,n);
	int i,j;
	for(j=1;j<=NumCols(EVal);j++){
		if(j+1>NumCols(EVal) || EVal[j+1][j] == 0){ // real
			CVal[j][j]=EVal[j][j];
			for(i=1;i<=m;i++){
				CVec[i][j]=complex<T>(EVec[i][j],0);
			}
		}
		else{ // complex
			CVal[j][j]=complex<T>(EVal[j][j],EVal[j+1][j]);
			for(i=1;i<=m;i++){
				CVec[i][j]=complex<T>(EVec[i][j],EVec[i][j+1]);
			}
			j++;
			CVal[j][j]=complex<T>(EVal[j][j],EVal[j-1][j]);
			for(i=1;i<=m;i++){
				CVec[i][j]=complex<T>(EVec[i][j-1],-EVec[i][j]);
			}
		}
	}
	return 1;
}


template <typename T> int Matrix<T>::Complex2Compact(Matrix<complex<T> >& CVec, Matrix<complex<T> >& CVal, Matrix<T>& EVec, Matrix<T>& EVal) {
	int m=NumRows(CVal);
	int n=NumCols(CVal);
	EVal.Resize(m,n); EVal=T();
	EVec.Resize(m,n);
	int i,j;
	for(j=1;j<=NumCols(CVal);j++){
		if(j+1>NumCols(CVal) || imag(CVal[j][j]) == 0){ // real
			EVal[j][j]=real(CVal[j][j]);
			for(i=1;i<=m;i++){
				EVec[i][j]=real(CVec[i][j]);
			}
		}
		else{ // complex
			EVal[j][j]=real(CVal[j][j]);
			EVal[j+1][j]=imag(CVal[j][j]);
			for(i=1;i<=m;i++){
				EVec[i][j]=real(CVec[i][j]);
			}
			j++;
			EVal[j][j]=real(CVal[j][j]);
			EVal[j-1][j]=imag(CVal[j][j]);
			if(EVal[j-1][j]+EVal[j][j-1] != 0){
				cerr << "ERROR: Input matrix is not in the compact form.\n";
				assert(0);
			}
			for(i=1;i<=m;i++){
				EVec[i][j]=imag(CVec[i][j-1]);
			}
		}
	}
	return 1;
}


// Computes sqrt(|a|^2 + |b|^2) without underflow or overflow
template <typename T> double Matrix<T>::pythagoras(const T a, const T b) {
	double x=Abs(a);
	double y=Abs(b);
	if (x > y){ // x != 0
		double w=y/x;
		return (x*sqrt(1.0+w*w));
	}
	else if(y == 0){ // x=y=0
		return 0;
	}
	else{ // y >= x, x != 0, y != 0
		double w=x/y;
		return (y*sqrt(1.0+w*w));
	}
}


// Computes sqrt(|a|^2 + |b|^2 + |c|^2) without underflow or overflow
template <typename T> double Matrix<T>::pythagoras3(const T a, const T b, const T c) {
	double x=Abs(a);
	double y=Abs(b);
	double z=Abs(c);
	double maxx=local_max(local_max(x,y),z);
	if(maxx==0)
		return 0;
	else
		return (maxx*sqrt((x/maxx)*(x/maxx)+(y/maxx)*(y/maxx)+(z/maxx)*(z/maxx)));
}


//Compute Jacobi rotation
template <typename T> void Matrix<T>::jacobiRotate(Matrix<T>& a,int i,int j,int k,int l,double s,double tau) {
	T g=a[i][j];
	T h=a[k][l];
	a[i][j]=g-s*(h+g*tau);
	a[k][l]=h+s*(g-h*tau);
}


// Jacobi-Transformation
template <typename T> int Matrix<T>::jacobiTransform(Matrix<T>& a, Vector<T>& d, Matrix<T>& v,int *rotN){
	int maxIter=100;
	int isOK=1;
	int j,iq,ip,i;
	double thresh,theta,tau,t,sm,s,h,g,c;
	int n=NumCols(a);
	Vector<T> b(n);
	Vector<T> z(n);
	for(ip=1;ip<=n;ip++){
		for(iq=1;iq<=n;iq++) v[ip][iq]=0;
		v[ip][ip]=1;
	}
	for(ip=1;ip<=n;ip++){
		b[ip]=d[ip]=a[ip][ip];
		z[ip]=0;
	}
	*rotN=0;
	for(i=1;i<=maxIter;i++){
		sm=0;
		for(ip=1;ip<=n-1;ip++){
			for(iq=ip+1;iq<=n;iq++)
				sm += Abs(a[ip][iq]);
		}
		if(sm == 0.0){
			isOK=1;
			return isOK;
		}
		if(i<4)
			thresh=0.2*sm/(n*n);
		else
			thresh=0.0;
		for(ip=1;ip<=n-1;ip++){
			for(iq=ip+1;iq<=n;iq++){
				g=100.0*Abs(a[ip][iq]);
				if(i>4 && (float)(Abs(d[ip])+g) == (float)Abs(d[ip])
					&& (float)(Abs(d[iq])+g) == (float)Abs(d[iq]))
					a[ip][iq]=0;
				else if(Abs(a[ip][iq])>thresh){
					h=d[iq]-d[ip];
					if((float)(Abs(h)+g) == (float)Abs(h))
						t=(a[ip][iq])/h;
					else{
						theta=0.5*h/(a[ip][iq]);
						t=1.0/(Abs(theta)+sqrt(1.0+theta*theta));
						if(theta<0) t = -t;
					}
					c=1.0/sqrt(1+t*t);
					s=t*c;
					tau=s/(1.0+c);
					h=t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq]=0;
					for(j=1;j<=ip-1;j++){
						Matrix<T>::jacobiRotate(a,j,ip,j,iq,s,tau);
					}
					for(j=ip+1;j<=iq-1;j++){
						Matrix<T>::jacobiRotate(a,ip,j,j,iq,s,tau);
					}
					for(j=iq+1;j<=n;j++){
						Matrix<T>::jacobiRotate(a,ip,j,iq,j,s,tau);
					}
					for(j=1;j<=n;j++){
						Matrix<T>::jacobiRotate(v,j,ip,j,iq,s,tau);
					}
					(*rotN)++;
				}
			}
		}
		for(ip=1;ip<=n;ip++){
			b[ip] += z[ip];
			d[ip]=b[ip];
			z[ip]=0;
		}
	}
	isOK=0; // did not converge in maxIter.
	return isOK;
}


// Eigen analysis by the Jacobi method
// Note: The input matrix A must be real symmetric.
template <typename T> int EigJacobi(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal) {
	assert(NumRows(A)==NumCols(A));
	double sym_tolerance=1e-3;
	int dim=NumRows(A);
	assert(A.IsReal() && A.IsSymmetric(sym_tolerance));
	Matrix<T> C=A;
	Vector<T> d(dim);
	EVec.Resize(dim,dim);
	int rotN;
	int isOK=Matrix<T>::jacobiTransform(C,d,EVec,&rotN);
	if(isOK==0){
		cerr<<"Eig: Not converged in EigJacobi()\n";
		return isOK;
	}
	
	int i;
	EVal.Resize(dim,dim); EVal=T();
	for(i=1;i<=dim;i++){
		EVal[i][i]=d[i];
	}
	return isOK;
}


// Eigen analysis by the Householder transformation and QL algorithm
// Note: The input matrix A must be real symmetric.
// The combination of tred2 and tqli is the most efficient technique to find all
// eigenvalues and eigenvectors of a real symmetric matrix.
template <typename T> int EigHHolder(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal) {
	assert(NumRows(A)==NumCols(A));
	int dim=NumRows(A);
	double sym_tolerance=1e-3;
	assert(A.IsReal() && A.IsSymmetric(sym_tolerance));
	int isOK;
#ifdef USE_NRC_CODE
	EVec=A;
	Vector<T> d(dim),e(dim);
	isOK=Matrix<T>::nrc_tred2(EVec,d,e,1);
	if(isOK==0){
		cerr<<"Eig: Not converged in EigHHolder()\n";
		return isOK;
	}
	isOK=Matrix<T>::nrc_tqli(d,e,EVec,1);
	if(isOK==0){
		cerr<<"Eig: Not converged in EigHHolder()\n";
		return isOK;
	}
#else
	Vector<T> d(dim);
	EVec.Resize(dim,dim);
	isOK=Matrix<T>::mes_symmeig(A,EVec,d);
#endif
	int i;
	EVal.Resize(dim,dim); EVal=T();
	for(i=1;i<=dim;i++){
		EVal[i][i]=d[i];
	}
	return isOK;
}


inline int EigHHolder(const Matrix<float>& A, Matrix<float>& EVec, Matrix<float>& EVal) {
	Matrix<double> AA=A.ToDouble();
	Matrix<double> EEVec=EVec.ToDouble();
	Matrix<double> EEVal=EVal.ToDouble();
	int isOK=EigHHolder(AA,EEVec,EEVal);
	EVec.FromDouble(EEVec);
	EVal.FromDouble(EEVal);
	return isOK;
}


template <typename T> int EigSymmetric(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal) {
	int dim=NumRows(A);
	if(A.IsDouble()==false){
			Matrix<double> AA=A.ToDouble();
			Matrix<double> EEVec;
			Matrix<double> EEVal;
			int isOK=EigSymmetric(AA,EEVec,EEVal);
			EVec.FromDouble(EEVec);
			EVal.FromDouble(EEVal);
			return isOK;
	}
	if(dim < 10)
		return EigJacobi(A,EVec,EVal);
	else
		return EigHHolder(A,EVec,EVal);
}


#ifndef DISABLE_COMPLEX
inline int EigSymmetric(const Matrix<complex<double> >& A, Matrix<complex<double> >& EVec, Matrix<complex<double> >& EVal) {
	cerr << "ERROR Must not arrive here!\n";
	assert(0);
	return 0;
}


inline int EigSymmetric(const Matrix<complex<float> >& A, Matrix<complex<float> >& EVec, Matrix<complex<float> >& EVal) {
	cerr << "ERROR Must not arrive here!\n";
	assert(0);
	return 0;
}
#endif


// the main eigen routine
template <typename T> int Eig(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal) {
	assert(NumRows(A)==NumCols(A));
	double sym_tolerance=1e-3;
	int dim=NumRows(A);
	if(A.IsReal() && A.IsSymmetric(sym_tolerance)){ // real symmetric	
		return EigSymmetric(A,EVec,EVal);
	}
	else if(A.IsComplex() && A.IsSymmetric(sym_tolerance)){ // complex symmetric (Hermitian)
		return EigHermitian(A,EVec,EVal);
	}
	else if(A.IsReal()){ // real general
		return EigReal(A,EVec,EVal);
	}
	else{ // complex general
		return EigComplex(A,EVec,EVal);
	}
}


// the sorted eigen routine
template <typename T> int EigS(const Matrix<T>& A, Matrix<T>& EVec, Matrix<T>& EVal) {
	int isOK=Eig(A,EVec,EVal);
	Matrix<T>::sortEigen(EVec,EVal);
	return isOK;
}


///////////////////////////////////////////////////////////////////////////////
// Simultaneous diagonalization
// Given two symmetric matrices A and B,
// find a transformation matrix E and a diagonal matrix D such that
//   E'*A*E=I
//   E'*B*E=D
// Returns 0 if input matrices are singular or nonsymmetric.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int SimDiag(const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& E, Matrix<T>& D) {
	int dim=NumRows(A);
	double thresh=1e-6;
	double sym_tolerance=1e-3;
	double zero_tolerance=1e-6;
	// test dimensions and symmetry
	int retCode;
	retCode = A.IsSymmetric(sym_tolerance);
	assert(retCode);
	retCode = B.IsSymmetric(sym_tolerance);
	assert(retCode);
	
	Matrix<T> W(dim,dim),P1(dim,dim);
	// eigenvalues of W in W, eigenvector in P1
	retCode = Eig(A,P1,W);
	assert(retCode); if(retCode==0) return 0;
	retCode = P1.NormCols();
	assert(retCode); if(retCode==0) return 0;
	
	// ------ whitening -----
	for (int i=1; i<=dim; i++) {
		double d = W[i][i];
		if (d > 0.0) W[i][i]= 1.0/sqrt(d);
		else if (d < 0.0) W[i][i]= -1.0/sqrt(-d);
	}
	W = P1 * W;	

	// ------ Here, A is transformed to I and B is transformed K.
	Matrix<T> K = Transpose(W) * B * W;
	int zeroN = Matrix<T>::cleanNearZero(K,zero_tolerance);
	// ----- new eigenvalues in K, eigenvector in P2 -----
	Matrix<T> P2(dim,dim);
	D.Resize(dim,dim);
	retCode = Eig(K,P2,D);
	assert(retCode); if(retCode==0) return 0;
	retCode = P2.NormCols();
	assert(retCode); if(retCode==0) return 0;

	E.Resize(dim,dim);
	E = W * P2; // E=EVec2*Whiten*EVec1
	retCode = Matrix<T>::sortEigen(E,D);
	assert(retCode); if(retCode==0) return 0;
	return retCode;
}


// A component of the diag vector 's' of the SVD of 'a' was smaller than 1.0e-6 and consequently set to zero.
template <typename T> int Matrix<T>::svdClean(Vector<T>& s, double tolerance) {
	int i;
	int n = s.Size();
	double smax=0;
	for(i=1;i<=n;i++)
		if(smax<Abs(s[i])) smax=Abs(s[i]);
    double smin = Abs(smax * tolerance);
    for (i=1; i<=n; i++) if ( Abs(s[i]) < smin) { s[i] = 0.0;}
	return 1;
}


template <typename T> int SvdRecipe(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V) {	
	int i,j;
	int m=A.m, n=A.n;
	int dim=local_max(m,n);

	// Append zeros to deal with non-square matrices
	U.Resize(dim,dim);
	Vector<T> d(dim);
	V.Resize(dim,dim);
	U.Input(A);
	double thresh=1e-6;
	int zeroN=Matrix<T>::cleanNearZero(U,thresh); //To eliminate underflow in SVD.

	// Call SVD main routine
	int a1=Matrix<T>::nrc_svdcmp(U, d, V);
	if(a1==0) {
		U=0; U.Input(A);
		zeroN=Matrix<T>::cleanNearZero(U,10*thresh); //To eliminate underflow in SVD.
		int maxIter=500;
		a1=Matrix<T>::nrc_svdcmp(U, d, V, maxIter);
		if(a1==0){
			U.Resize(m,m); U=0;
			S.Resize(m,n); S=0;
			V.Resize(n,n); V=0;
			cerr << "ERROR: No convergence in " << maxIter << " svdcmp iterations" << endl;
			return a1;
		}
	}
	Matrix<T>::svdClean(d);

	// Bookkeep the results
	Vector<T> Stmp;
	Vector<int> order;
	VecSort(d, Stmp, order, -1);
	Matrix<T> Utmp(dim,dim), Vtmp(dim,dim);
	for(i=1;i<=dim;i++) for(j=1;j<=dim;j++) {
		Utmp[i][j]=U[i][order[j]];
		Vtmp[i][j]=V[i][order[j]];
	}
	Utmp.Resize(m,m);U.Resize(m,m);
	Vtmp.Resize(n,n);V.Resize(n,n);
	U=Utmp;V=Vtmp;
	S.Resize(m,n); S=0;
	for(i=1;i<=local_min(m,n);i++) S[i][i]=Stmp[i];
	return a1;
}


template <typename T> int SvdSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	double thresh=1e-6;
	assert(A.m>=A.n);
	int m=A.m, n=A.n;
	int dim=local_max(m,n);
#ifdef USE_NRC_CODE
	Matrix<T> U(dim,dim),V(dim,dim);
	Vector<T> d(dim);
	U=0; U.Input(A);
	int zeroN=Matrix<T>::cleanNearZero(U,thresh); //To eliminate underflow in SVD.
	x.Resize(n);
	int a1=Matrix<T>::nrc_svdcmp(U, d, V);
	Matrix<T>::svdClean(d);
	return Matrix<T>::nrc_svsolv(U,d,V,b,x);
#else
	return QrSolve(A,b,x);
#endif
}

///////////////////////////////////////////////////////////////////////////////
// the main SVD routine
// Given a matrix A(mxn), this routines computes its singular value
// decomposition, A = U*S*V'. The matrix U replaces A on output. The
// diagonal matrix of singular values S(nxn) is output as (1xn). The
// matrix V (not the Transpose V') is output as (nxn).
// Note: Diagonal elements are sorted in the non-increasing order.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Svd(const Matrix<T> A, Matrix<T>& U, Matrix<T>& S, Matrix<T>& V) {	
#ifdef USE_NRC_CODE
	return SvdRecipe(A,U,S,V);
#else
	return SvdMeschach(A,U,S,V);
#endif
}


#ifndef DISABLE_COMPLEX
inline int SvdComplex(const Matrix<complex<float> > A, Matrix<complex<float> >& U, Matrix<complex<float> >& S, Matrix<complex<float> >& V) {	
	static Matrix<complex<double> > AA,UU,SS,VV;
	AA=A.ToCDouble();
	int isOK=Svd(A,U,S,V);
	U.FromCDouble(UU);
	S.FromCDouble(SS);
	V.FromCDouble(VV);
	return isOK;
}
#endif


// Orthogonal space
template <typename T> Matrix<T> Orth(const Matrix<T>& A) {
	Matrix<T> U,S,V;
	int a1=Svd(A,U,S,V);
	int r=0;
	for(int i=1;i<=NumRows(S);i++) if(S[i][i] != 0) r++;
	return U.SubMatrix(1,NumRows(U),1,r);
}


// Null space
template <typename T> Matrix<T> Null(const Matrix<T>& A) {
	Matrix<T> U,S,V;
	int a1=Svd(A,U,S,V);
	int r=0;
	for(int i=1;i<=NumRows(S);i++) if(S[i][i] != 0) r++;
	return V.SubMatrix(1,NumRows(V),NumCols(V)-r,NumCols(V));
}


// Matrix inversion by using SVD
// Using the singular value decomposition this routines returns B the
// inverse of a square matrix A or at least the best choice in a
// least-squares sense. Matrix A won't be changed unless B=A.
// A = U *  S   *  V'     result of the SVD
// B = V * S^-1 *  U'     inverse (note: S diagonal, U V orthonormal)
template <typename T> Matrix<T> PinvRecipe(const Matrix<T>& A, int* retCode=0) {
	Matrix<T> U,S,V,Result;	
	*retCode=1;
	SvdRecipe(A, U, S, V);
	int isOK=Matrix<T>::pinv(U,S,V,Result);
	if(retCode) *retCode=isOK;
	return Result;
}


template <typename T> Matrix<T> Pinv(const Matrix<T>& A, int* retCode=0) {
#ifdef USE_NRC_CODE
	return PinvRecipe(A,retCode);
#else
	return PinvMeschach(A,retCode);
#endif
}


template <typename T> int Matrix<T>::UpperBackSubstitute(const Matrix<T>& A,Vector<T>& b){
	assert(A.IsUpperTri());
	int i,j;
	T sum;
	b[A.n] /= A[A.n][A.n];
	for(i=A.m-1;i>=1;i--){
		for(sum=0,j=i+1;j<=A.n;j++) sum += A[i][j]*b[j];
		b[i]=(b[i]-sum)/A[i][i];
	}
	return 1;
}


template <typename T> int Matrix<T>::LowerBackSubstitute(const Matrix<T>& A,Vector<T>& b){
	assert(A.IsLowerTri());
	int i,j;
	T sum;
	b[1] /= A[1][1];
	for(i=2;i<=A.m;i++){
		for(sum=0,j=1;j<i;j++) sum += A[i][j]*b[j];
		b[i]=(b[i]-sum)/A[i][i];
	}
	return 1;
}


// ------------------------------------------------------------------------
// Lu
//    LU Decomposition
//	A=L*U;
//	where A is an NxN square matrix and L is a lower triangular matrix 
//	and U is an upper triangular matrix.
// Det
//    Determinant
// LogAbsDet
//	  Log of absolute determinant
// ------------------------------------------------------------------------
///////////////////////
// LU decomposition
// L*U=P*A;
///////////////////////
template <typename T> int Lu(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U, Matrix<T>& P) {
	assert(NumRows(A)==NumCols(A));
	int n=NumRows(A);
	int i,j,isOK;
#ifdef USE_NRC_CODE
	Vector<int> indx(n),indx2(n);
	Matrix<T> B=A;
	double d=0;
	isOK=Matrix<T>::nrc_ludcmp(B,indx,&d,indx2);
	if(isOK==0) return isOK;
	L.Resize(n,n); L=0;
	U.Resize(n,n); U=0;
	for(i=1;i<=n;i++){
		for(j=i;j<=n;j++){
			U[i][j]=B[i][j];
		}
		for(j=1;j<i;j++){
			L[i][j]=B[i][j];
		}
		L[i][i]=1;
	}
	P.Resize(n,n); P=0;
	for(i=1;i<=n;i++) P[i][indx2[i]]=1;
	return isOK;
#else
	Vector<int> pivot(n);
	Matrix<T> B=A;
	isOK=Matrix<T>::mes_LUfactor(B,pivot);
	L.Resize(n,n); L=T();
	U.Resize(n,n); U=T();
	for(i=1;i<=n;i++){
		for(j=i;j<=n;j++){
			U[i][j]=B[i][j];
		}
		for(j=1;j<i;j++){
			L[i][j]=B[i][j];
		}
		L[i][i]=1;
	}
	P.Resize(n,n); P=T();
	for(i=1;i<=n;i++) P[i][pivot[i]]=1;
	return isOK;
#endif
	
}


///////////////////////
// LU decomposition
// A=L*U;
///////////////////////
template <typename T> int Lu(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U) {
#ifdef USE_NRC_CODE
	Matrix<T> P;
	int isOK=Lu(A,L,U,P);
	if(isOK==0) return isOK;
	L=Transpose(P)*L;
#else
	Matrix<T> P;
	int isOK=Lu(A,L,U,P);
	if(isOK==0) return isOK;
	L=Transpose(P)*L;
#endif
	return isOK;
}


///////////////////////
// LU solution
// Solves A*x=b
///////////////////////
template <typename T> int LuSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	assert(A.m==A.n);
	int isOK;
	Matrix<T> L=A;
#ifdef USE_NRC_CODE
	Vector<int> indx(A.n),indx2(A.n);
	double d=0;
	isOK=Matrix<T>::nrc_ludcmp(L,indx,&d,indx2);
	if(isOK==0){
		cerr << "Input matrix A is singular.\n";
		return isOK;
	}
	isOK=Matrix<T>::nrc_lusolv(L,indx,b,x);
	return isOK;
#else
	Vector<int> pivot(A.n);
	isOK=Matrix<T>::mes_LUfactor(L,pivot);
	isOK=Matrix<T>::mes_LUsolve(L,pivot,b,x);
	return isOK;
#endif
}


///////////////////////
// Determinant
///////////////////////
template <typename T> T Det(const Matrix<T>& A) {
	assert(NumRows(A)==NumCols(A));
	int n=NumRows(A);
#ifdef USE_NRC_CODE
	Vector<int> indx(n),indx2(n);
	double d;
	Matrix<T> B=A;
	int a1=Matrix<T>::nrc_ludcmp(B,indx,&d,indx2);
	if(a1==0){
		cerr << "ERROR: Singular matrix in Det()\n";
		assert(0);
		return a1;
	}
#else
	Vector<int> pivot(n);
	Matrix<T> B=A;
	int a1=Matrix<T>::mes_LUfactor(B,pivot);
	if(a1==0){
		cerr << "ERROR: Singular matrix in Det()\n";
		assert(0);
		return a1;
	}
	int d = Matrix<T>::mes_px_sign(pivot);
#endif
	T det=T(d);
	for(int i=1;i<=n;i++) {
		det*=(B[i][i]);
	}
	return det;
}


///////////////////////
// Log Absolute Determinant
///////////////////////
template <typename T> double LogAbsDet(const Matrix<T>& A) {
	assert(NumRows(A)==NumCols(A));
	int n=NumRows(A);
#ifdef USE_NRC_CODE
	Vector<int> indx(n),indx2(n);
	double d;
	Matrix<T> B=A;
	int a1=Matrix<T>::nrc_ludcmp(B,indx,&d,indx2);
	if(a1==0){
		cerr << "ERROR: Singular matrix in Det()\n";
		assert(0);
		return 0;
	}
#else
	Vector<int> pivot(n);
	Matrix<T> B=A;
	int a1=Matrix<T>::mes_LUfactor(B,pivot);
	if(a1==0){
		cerr << "ERROR: Singular matrix in Det()\n";
		assert(0);
		return 0;
	}
#endif
	double logdet=0;
	for(int i=1;i<=n;i++) {
		logdet+=Log(Abs(B[i][i]));
	}
	return logdet;
}


///////////////////////
// Inverse
///////////////////////
template <typename T> Matrix<T> Inv(const Matrix<T>& A, int* retCode=0) {
	int m=NumRows(A);
	int n=NumCols(A);
	assert(m==n);
#ifdef USE_NRC_CODE
	Vector<int> indx(n),indx2(n);
	double d=0;
	Matrix<T> B=A;
	if(retCode) *retCode=1;
	Matrix<T> Y(m,n);
	int a1=Matrix<T>::nrc_ludcmp(B,indx,&d,indx2);
	if(a1==0) {
		if(retCode) *retCode=a1;
		return Y;
	}
	int i,j;
	Vector<T> col(n);
	Vector<T> x(n);
	for(j=1;j<=n;j++){
		for(i=1;i<=n;i++) col[i]=0;
		col[j]=1;
		Matrix<T>::nrc_lusolv(B,indx,col,x);
		for(i=1;i<=n;i++) Y[i][j]=x[i];
	}
	return Y;
#else
	Matrix<T> Y(m,n);
	int isOK=Matrix<T>::mes_inverse(A,Y);
	if(retCode) *retCode=isOK;
	return Y;
#endif
}


///////////////////////////////////////////////////////////////
// Condition number
///////////////////////////////////////////////////////////////
// For symmetric matrices
// k=EVal(max)/EVal(min)  
///////////////////////
template <typename T> double Cond(const Matrix<T>& A, int p=2) {
#ifdef USE_NRC_CODE
	int m=NumRows(A);
	int n=NumRows(A);
	assert(m==n);
	assert(p==2);
	if(A.IsSymmetric()){
		Matrix<T> EVec, EVal;
		if(EigS(A,EVec,EVal)){
			if(EVal[m][m]==0) return mtl_numeric_limits<double>::max();
			else return Abs(EVal[1][1]/EVal[m][m]);
		}
		else{
			cerr << "ERROR: Failed in computing eigenvalues.\n";
			return mtl_numeric_limits<double>::max();
		}
	}
	else{
		cerr << "ERROR: In the current implementation, condition number is computed for symmetric matrices.\n";
		assert(0);
		return HUGE_VAL;
	}
#else
	return CondQR(A,p);
#endif
}


#if !defined(USE_NRC_CODE)
///////////////////////
// LU Condition number
///////////////////////
template <typename T> double CondLU(const Matrix<T>& A) {
	Vector<int> pivot;
	Matrix<T> LU=A;
	Matrix<T>::mes_LUfactor(LU,pivot);
	return Matrix<T>::mes_LUcondest(LU,pivot);
}


///////////////////////
// QR Condition number
///////////////////////
template <typename T> double CondQR(const Matrix<T>& A, int p=2) {
	Matrix<T> QR=A;
	Vector<T> diag;
	Matrix<T>::mes_QRfactor(QR,diag);
	return Matrix<T>::mes_QRcondest(QR,p);
}
#endif


// ------------------------------------------------------------------------
// Cholesky decomposition
// 	A=L*L' where the input A is a positive-definite symmetric NxN matrix
// 	and L is a lower triangular matrix.
// ------------------------------------------------------------------------
///////////////////////
// Cholesky decomposition
// A=L*L';
///////////////////////
template <typename T> int Chol(const Matrix<T>& A, Matrix<T>& L) {
	assert(A.m==A.n);
	int i,j;
#ifdef USE_NRC_CODE
	Vector<T> d(A.m);
	L=A;
	int isOK=Matrix<T>::nrc_choldcmp(L,d);
	for(i=1;i<=A.m;i++){
		L[i][i]=d[i]; for(j=i+1;j<=A.n;j++) L[i][j]=0;
	}
	return isOK;
#else
	L=A;
	int isOK=Matrix<T>::mes_CHfactor(L);
	for(i=1;i<=A.m;i++){
		for(j=i+1;j<=A.n;j++) L[i][j]=0;
	}
	return isOK;
#endif
}


///////////////////////
// Cholesky solution
// Solves A*x=b
///////////////////////
template <typename T> int CholSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	assert(A.m==A.n);
	int isOK;
#ifdef USE_NRC_CODE
	Vector<T> d(A.m);
	Matrix<T> L=A;
	isOK=Matrix<T>::nrc_choldcmp(L,d);
	if(isOK==0){
		cerr << "Input matrix A is not positive-definite symmetric.\n";
		return isOK;
	}
	x.Resize(A.m);
	isOK=Matrix<T>::nrc_cholsolv(L,d,b,x);
	return isOK;
#else
	Matrix<T> L=A;
	isOK=Matrix<T>::mes_CHfactor(L);
	isOK=Matrix<T>::mes_CHsolve(L,b,x);
	return isOK;
#endif
}


///////////////////////
// LDL decomposition
// A=L*D*L';
///////////////////////
template <typename T> int Ldl(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& D) {
	assert(A.m==A.n);
	int i,j,isOK;
#ifdef USE_NRC_CODE
	isOK=Chol(A,L);
	D.Resize(A.m,A.n); D=0;
	for(i=1;i<=A.m;i++){
		T d=L(i,i);
		D(i,i)=d*d;
		L(i,i)=1;
		for(j=1;j<i;j++)
			L(i,j)/=d;
		for(j=i+1;j<=A.n;j++)
			L(i,j)=0;
	}
	return isOK;
#else
	L=A;
	isOK=Matrix<T>::mes_LDLfactor(L);
	D.Resize(A.m,A.n); D=0;
	for(i=1;i<=A.m;i++){
		D(i,i)=L(i,i);
		L(i,i)=1;
		for(j=i+1;j<=A.n;j++)
			L(i,j)=0;
	}
	return isOK;
#endif
}


///////////////////////
// LDL solution
// Solves A*x=b
///////////////////////
template <typename T> int LdlSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	assert(A.m==A.n);
	int isOK;
#ifdef USE_NRC_CODE	
	isOK=CholSolve(A,b,x);
	return isOK;
#else
	Matrix<T> LDL=A;
	LDL=A;
	isOK=Matrix<T>::mes_LDLfactor(LDL);
	isOK=Matrix<T>::mes_LDLsolve(LDL,b,x);
	return isOK;
#endif
}


// ------------------------------------------------------------------------
// QR decomposition
//	A=Q*R where the input A is an arbitrary MxN matrix
//	and Q is an orthogonal MxM matrix and R is an upper triangular MxN matrix.
//	Q'*Q=I
// ------------------------------------------------------------------------
///////////////////////
// QR decomposition
// A=Q*R;
///////////////////////
template <typename T> int Qr(const Matrix<T>& A, Matrix<T>& Q,Matrix<T>& R) {
	int dim=local_max(A.m,A.n);
	int isOK;
#ifdef USE_NRC_CODE
	Vector<T> c(dim),d(dim);
	R.Resize(dim,dim);
	R=0; R.Input(A);
	isOK=Matrix<T>::nrc_qrdcmp(R,c,d);
	Q.Resize(dim,dim); Q.SetIdentity();
	Matrix<T> I(dim,dim); I.SetIdentity();
	Vector<T> uj(dim);
	int i,j;
	for(j=1;j<=dim-1;j++){
		for(i=1;i<=j-1;i++) uj[i]=0;
		for(i=j;i<=dim;i++) uj[i]=R[i][j];
		Q=Q*(I-OuterProd(uj,uj)/c[j]);
	}
	if(A.m<dim) Q.Resize(A.m,A.m);
	for(i=1;i<=A.m;i++){
		for(j=1;j<i;j++) R[i][j]=0;
		R[i][i]=d[i];
	}
	R.Resize(A.m,A.n);
#else
	Vector<T> diag(dim);
	R=A;
	isOK=Matrix<T>::mes_QRfactor(R,diag);
	isOK=Matrix<T>::mes_makeQ(R,diag,Q);
	isOK=Matrix<T>::mes_makeR(R,R);
#endif
	return isOK;
}


///////////////////////
// QR solution
// Solves A*x=b
///////////////////////
template <typename T> int QrSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	assert(A.m==b.Size());
	int isOK;
	int dim=local_max(A.m,A.n);
	// Ax=b, QRx=b, Rx=Q'b
	if(A.m>=A.n){
		// In this case, this routine gives the least-square solution.
#ifdef USE_NRC_CODE
		Vector<T> c(dim),d(dim);
		Matrix<T> QR(dim,dim); QR=0; QR.Input(A);
		isOK=Matrix<T>::nrc_qrdcmp(QR,c,d);
		x.Resize(A.n);
		isOK=Matrix<T>::nrc_qrsolv(QR,c,d,b,x);
		return isOK;
#else
		Vector<T> diag(dim);
		Matrix<T> QR=A;
		isOK=Matrix<T>::mes_QRfactor(QR,diag);
		isOK=Matrix<T>::mes_QRsolve(QR,diag,b,x);
		return isOK;
#endif
	}
	else{
		// This implementation is just to produce the same answer as MATLAB.
		// This solution is not equal to Pinv(A)*b.
		// Please see the explanation in the MATLAB reference.
		Matrix<T> Q,R;
		int isOK=Qr(A,Q,R); // A=QR
		Vector<T> b2(A.m);
		int i,j;
		for(i=1;i<=A.m;i++){ // Q'*b
			T sum=0;
			for(j=1;j<=A.m;j++) sum+=conj(Q[j][i])*b[j];
			b2[i]=sum;
		}
		Matrix<T> A2=R.SubMatrix(1,A.m,A.n-A.m+1,A.n);
		Vector<T> x2;
		isOK=QrSolve(A2,b2,x2);
		x.Resize(A.n);
		for(i=1;i<=A.n-A.m;i++) x[i]=0;
		for(i=A.n-A.m+1,j=1;i<=A.n;i++,j++) x[i]=x2[j];
		return isOK;
	}
}


///////////////////////
// QRCP decomposition
// A=;
///////////////////////
template <typename T> int Qrcp(const Matrix<T>& A, Matrix<T>& QRCP, Vector<T>& diag, Vector<int>& pivot) {
	QRCP=A;
	int isOK=Matrix<T>::mes_QRCPfactor(QRCP,diag,pivot);
	return isOK;
}


///////////////////////
// QRCP solution
// Solves A*x=b
///////////////////////
template <typename T> int QrcpSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	assert(A.m==b.Size());
	Matrix<T> QR=A;
	Vector<T> diag;
	Vector<int> px;
	int isOK=Matrix<T>::mes_QRCPfactor(QR,diag,px);
	if(isOK==0) return isOK;
	isOK=Matrix<T>::mes_QRCPsolve(QR,diag,px,b,x);
	return isOK;
}


///////////////////////
// BKP decomposition
// A=;
///////////////////////
template <typename T> int Bkp(const Matrix<T>& A, Matrix<T>& BKP, Vector<int>& pivot, Vector<int>& blocks) {
	BKP=A;
	int isOK=Matrix<T>::mes_BKPfactor(BKP,pivot,blocks);
	return isOK;
}


///////////////////////
// BKP solution
// Solves A*x=b
///////////////////////
template <typename T> int BkpSolve(const Matrix<T>& A, const Vector<T>& b, Vector<T>& x) {
	assert(A.m==b.Size());
	Matrix<T> BKP=A;
	Vector<int> blocks;
	Vector<int> px;
	int isOK=Matrix<T>::mes_BKPfactor(BKP,px,blocks);
	if(isOK==0) return isOK;
	isOK=Matrix<T>::mes_BKPsolve(BKP,px,blocks,b,x);
	return isOK;
}


// ------------------------------------------------------------------------
// Solve
//	Solve linear equations for unknown x: Ax=b
//	If A is triangular matrix, backsubstitution is used.
//	If A is symmetric positive definite, the Cholesky factorization is used.
//	If A is square, the LU decomposition is used.
//	If A is not square, the QR decomposition is used.
// ------------------------------------------------------------------------
template <typename T> int Solve(const Matrix<T>& A,const Vector<T>& b,Vector<T>& x) {
	if(A.IsUpperTri()){
		x=b; return Matrix<T>::UpperBackSubstitute(A,x);
	}
	else if(A.IsLowerTri()){
		x=b; return Matrix<T>::LowerBackSubstitute(A,x);
	}
	else if(A.IsSymmetric() && A.IsPD()){
		return CholSolve(A,b,x); // A=L*L', x=L'\(L\b)
	}
	else if(A.IsSquare()){
		return LuSolve(A,b,x); // A=L*U, x=U\(L'\b)
	}
	else{
		return QrSolve(A,b,x);
	}
	return 0;
}


// ------------------------------------------------------------------------
//Rank
//	Number of independent rows or columns
// ------------------------------------------------------------------------
template <typename T> int Rank(const Matrix<T>& A) {
	double eps=mtl_numeric_limits<double>::epsilon();
	Matrix<T> U,S,V;
	int isOK=Svd(A,U,S,V);
	double tol = Max(Size(A))*real(S[1][1])*eps;
	int r=0;
	for(int i=1;i<=local_min(NumRows(S),NumCols(S));i++){
		if(real(S[i][i])>tol) r++;
	}
	return r;
}


// EXPORT->Svd: Calculate the decompostion of matrix A.
//   NOTE: on return that U and V hold U' and V' respectively!
template <typename T> int SvdMeschach(const Matrix<T>& AA, Matrix<T>& UU, Matrix<T>& SS, Matrix<T>& VV)
{
	int m=NumRows(AA), n=NumCols(AA);
	Matrix<T> BB=AA;
	static Vector<T> dd;
	double thresh=1e-6;
	int zeroN=Matrix<T>::cleanNearZero(BB,thresh); //To eliminate underflow in SVD.
	int isOK=Matrix<T>::mes_svd(BB,UU,VV,dd);
	Matrix<T>::svdClean(dd);
	UU=Transpose(UU);
    VV=Transpose(VV);	
	SS.Resize(m,n); SS=T();
	for(int i=1;i<=local_min(m,n);i++)
		SS[i][i]=T(dd[i]);
	return isOK;
}

inline int SvdMeschach(const Matrix<float>& A, Matrix<float>& U, Matrix<float>& S, Matrix<float>& V)
{
	static Matrix<double> BB,UU,SS,VV;
	BB=A.ToDouble();
	int isOK=SvdMeschach(BB,UU,SS,VV);
	U.FromDouble(UU);
	S.FromDouble(SS);
	V.FromDouble(VV);
	return isOK;
}


#ifndef DISABLE_COMPLEX
inline int SvdMeschach(const Matrix<complex<float> >& A, Matrix<complex<float> >& U, Matrix<complex<float> >& S, Matrix<complex<float> >& V)
{
	static Matrix<complex<double> > BB,UU,SS,VV;
	BB=A.ToCDouble();
	int isOK=SvdMeschach(BB,UU,SS,VV);
	U.FromCDouble(UU);
	S.FromCDouble(SS);
	V.FromCDouble(VV);
	return isOK;
}
#endif


// EXPORT->Pinv: Inverted Singular Value Decomposition (calls Svd)
//   and inverse of A is returned in Result
template <typename T> Matrix<T> PinvMeschach(const Matrix<T>& A, int *retCode=0)
{
   static Matrix<T> U,S,V,Result;
   SvdMeschach(A, U, S, V);
   int isOK=Matrix<T>::pinv(U,S,V,Result);
   if(retCode) *retCode=isOK;
   return Result;
}


// pinv(A)=V * S^(-1) * U^T
template <typename T> int Matrix<T>::pinv(Matrix<T>& U, Matrix<T>& S, Matrix<T>& V, Matrix<T>& pinvA)
{
   int i, j;
   int m = NumRows(U);
   int n = NumCols(V);
   int dim=local_min(m,n);
   static Matrix<T> tmp1;
   tmp1.Resize(n,m);

   // tmp1 = S^(-1)*U^T
   // tmp1 will be the product of the matrix U^T and the diagonal 
   // matrix of singular values stored in vector w. tmp1 is then
   // pre-multiplied by the matrix v to produce the inverse which is returned
   for (i = 1; i <= n; i++){
	   double d;
	   if(i<=dim) d=real(S[i][i]);
	   else d=0;
	   // Only non-zero elements of diag matrix w are used to compute the inverse.
	   if(d>0) d=1/d;
	   for (j = 1; j <= m; j++){
		   if(i<=m)
			   tmp1[i][j] = T(d)*conj(U[j][i]);
		   else
			   tmp1[i][j]=0;
	   }
   }
   pinvA=V*tmp1;
   return 1;
}


template <typename T> int Schur(const Matrix<T>& AA, Matrix<T>& TT, Matrix<T>& QQ)
{
	if(AA.IsReal()){
		return SchurReal(AA,TT,QQ);
	}
	else{
		return SchurComplex(AA,TT,QQ);
	}
}


template <typename T> int SchurReal(const Matrix<T>& AA, Matrix<T>& TT, Matrix<T>& QQ)
{
	TT=AA;
	int isOK=Matrix<double>::mes_schur(TT,QQ); // compute Schur form: A = Q*T*Q^T
	return isOK;
}



inline int Schur(const Matrix<float>& A, Matrix<float>& t, Matrix<float>& Q)
{
	static Matrix<double> AA,TT,QQ;
	int m=NumRows(A);
	int n=NumCols(A);
	AA=A.ToDouble();
	QQ.Resize(m,n);
	int isOK=Schur(AA,TT,QQ); // compute Schur form: A = Q*T*Q^T
	t.FromDouble(TT);
	Q.FromDouble(QQ);
	return isOK;
}


#ifndef DISABLE_COMPLEX
inline int SchurReal(const Matrix<complex<float> >& AA, Matrix<complex<float> >& TT, Matrix<complex<float> >& QQ)
{
	cerr << "ERROR Must not arrive here!\n";
	assert(0);
	return 0;
}


inline int SchurReal(const Matrix<complex<double> >& AA, Matrix<complex<double> >& TT, Matrix<complex<double> >& QQ)
{
	cerr << "ERROR Must not arrive here!\n";
	assert(0);
	return 0;
}
#endif


template <typename T> int SchurComplex(const Matrix<T>& AA, Matrix<T>& TT, Matrix<T>& QQ)
{
	TT = AA;
	QQ.Resize(NumRows(AA),NumCols(AA));
	int isOK=Matrix<T>::mes_zschur(TT,QQ); // compute Schur form: A = Q*T*Q^T
	return isOK;
}


#ifndef DISABLE_COMPLEX
inline int SchurComplex(const Matrix<complex<float> >& A, Matrix<complex<float> >& t, Matrix<complex<float> >& Q)
{
	static Matrix<complex<double> > AA,TT,QQ;
	AA=A.ToCDouble();
	int isOK=SchurComplex(AA,TT,QQ);
	t.FromCDouble(TT);
	Q.FromCDouble(QQ);
	return isOK;
}
#endif


///////////////////////////////////////////////////////////////////////////////
// Hess
//
// A Hessenberg matrix is zero below the first subdiagonal.
// If the matrix is symmetric or Hermitian, the form is tridiagonal.
// This matrix has the same eigenvalues as the original,
// but less computation is needed to reveal them.
//    A = P*H*P^T  and P^T*P = eye(size(A))
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Hess(const Matrix<T>& A, Matrix<T>& H, Matrix<T>& P)
{
	static Vector<double> diag, beta;
	int n = NumCols(A);
    diag.Resize(n);
    beta.Resize(n);
	H=A;
    int isOK=Matrix<double>::mes_Hfactor(H,diag,beta);
	if(isOK==0) return isOK;
	isOK=Matrix<double>::mes_makeHQ(H,diag,beta,P);
	if(isOK==0) return isOK;
    isOK=Matrix<double>::mes_makeH(H,H);
	if(isOK==0) return isOK;
	return isOK;
}


inline int Hess(const Matrix<float>& A, Matrix<float>& H, Matrix<float>& P)
{
	static Matrix<double> AA,HH,PP;
	AA=A.ToDouble();
	int isOK=Hess(AA,HH,PP);
	H.FromDouble(HH);
	P.FromDouble(PP);
	return isOK;
}


#ifndef DISABLE_COMPLEX
template <typename T> int HessComplex(const Matrix<T>& A, Matrix<T>& H, Matrix<T>& P)
{
	static Vector<complex<double> > diag;
	int n = NumCols(A);
    diag.Resize(n);
	H=A;
	int isOK=Matrix<complex<double> >::mes_zHfactor(H,diag);
	if(isOK==0) return isOK;
    isOK=Matrix<complex<double> >::mes_zHQunpack(H,diag,P,H);
	if(isOK==0) return isOK;
	return isOK;
}


inline int HessComplex(const Matrix<complex<float> >& A, Matrix<complex<float> >& H, Matrix<complex<float> >& P)
{
	Matrix<complex<double> > AA,HH,PP;
	AA=A.ToCDouble();
	int isOK=HessComplex(AA,HH,PP);
	H.FromCDouble(HH);
	P.FromCDouble(PP);
	return isOK;
}
#endif


template <typename T> int EigReal(const Matrix<T>& AA, Matrix<T>& EEVec, Matrix<T>& EEVal)
{
	static Matrix<T> TT,QQ;
	int isOK=2;
	int m=NumRows(AA);
	int n=NumCols(AA);
	QQ.Resize(m,n);
	TT=AA;
	
	// compute Schur form: A = Q*T*Q^T
	isOK=Schur(AA,TT,QQ);
	if(isOK==0) return isOK;
	
	// extract eigenvalues
	Vector<T>   evals_re(m);
	Vector<T>   evals_im(m);
	isOK=Matrix<T>::mes_schur_evals(TT,evals_re,evals_im);
	if(isOK==0) return isOK;
	
	// Q not needed for eigenvalues
	Matrix<T> EEVec_re(m,n), EEVec_im(m,n);
	isOK=Matrix<T>::mes_schur_vecs(TT,QQ,EEVec_re,EEVec_im);
	// k'th eigenvector is k'th column of (X_re + i*X_im)

	// Arrange the results into the compact format.
	EEVal.Resize(m,n); EEVal=T();
	EEVec.Resize(m,n);
	int i,j;
	for(j=1;j<=n;j++){
		if(j+1>n || evals_im[j]==0){ // real eigenvalue
			EEVal[j][j]=evals_re[j];
			for(i=1;i<=m;i++) EEVec[i][j]=EEVec_re[i][j];
		}
		else{ // complex eigenvalue
			isOK=2; // Indicate the results are in the compact format.
			EEVal[j][j]=evals_re[j];
			EEVal[j+1][j]=evals_im[j];
			for(i=1;i<=m;i++) EEVec[i][j]=EEVec_re[i][j];
			j++;
			if(j>n){
				cerr << "ERROR: Eigenvalues are not in conjugate pairs.\n";
				assert(0);
				isOK=0;
				break;
			}
			EEVal[j][j]=evals_re[j];
			EEVal[j-1][j]=evals_im[j];
			if(EEVal[j-1][j]+EEVal[j][j-1] != 0){
				cerr << "ERROR: Eigenvalues are not in conjugate pairs.\n";
				assert(0);
				isOK=0;
				break;
			}
			for(i=1;i<=m;i++) EEVec[i][j]=EEVec_im[i][j-1];
		}
	}
	return isOK;
}


inline int EigReal(const Matrix<float>& A, Matrix<float>& EVec, Matrix<float>& EVal)
{
	Matrix<double> AA=A.ToDouble();
	Matrix<double> EEVec, EEVal;
	int isOK=EigReal(AA,EEVec, EEVal);
	EVec.FromDouble(EEVec);
	EVal.FromDouble(EEVal);
	return isOK;
}


#ifndef DISABLE_COMPLEX
inline int EigReal(const Matrix<complex<float> >& A, Matrix<complex<float> >& EVec, Matrix<complex<float> >& EVal)
{
	cerr << "ERROR Must not arrive here!\n";
	assert(0);
	return 0;
}


inline int EigReal(const Matrix<complex<double> >& A, Matrix<complex<double> >& EVec, Matrix<complex<double> >& EVal)
{
	cerr << "ERROR Must not arrive here!\n";
	assert(0);
	return 0;
}
#endif


///////////////////////////////////////////////////////////////////////////////
// Compute eigenvectors and eigenvalues of a Hermitian matrix.
//   C*U = D*U
// Let C=A+iB, U=u+iv
// (A+iB).(u+iv) = lambda(u+iv)
// The above problem is equivalent to a 2M x 2N real symmetric eigenvector problem. 
// 1. Form 2M x 2N real matrix.
//      [ A  -B ] [ u ] = lambda [ u ]
//      [ -B  A ] [ v ]          [ v ]
// 2. Solve the real symmetric eigenvector problem.
//    The eigenvalues are repeated as d1, d1, d2, d2, ..., dn, dn.
//    The eigenvectors are repeated as u+iv, i(u+iv).
// 3. Select eigenvectors and eigenvalues.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int EigHermitian(const Matrix<T>& AA, Matrix<T>& EEVec, Matrix<T>& EEVal)
{
	int m=NumRows(AA), n=NumCols(AA);

	// Form a 2M x 2N real matrix
	Matrix<double> AA_re(2*m,2*n);
	Matrix<double> EVec_re(2*m,2*n);
	Matrix<double> EVal_re(2*m,2*n);
	int i,j;
	for(i=1;i<=m;i++) for(j=1;j<=n;j++) {
		double t=real(AA[i][j]);
		AA_re[i][j]=t;
		AA_re[i+m][j+n]=t;
	}
	for(i=1;i<=m;i++) for(j=1;j<=n;j++) {
		double t=imag(AA[i][j]);
		AA_re[i+m][j]=t;
		AA_re[i][j+n]=-t;
	}

	// Call eigenvector routine for real symmetric case.
	int isOK=EigS(AA_re, EVec_re, EVal_re);

	EEVec.Resize(m,n); EEVec=T();
	EEVal.Resize(m,n); EEVal=T();
	for(j=1;j<=n;j++){
		EEVal[j][j]=T(EVal_re[2*j][2*j]);
		for(i=1;i<=m;i++) {
			EEVec[i][j]=T(EVec_re[i][2*j],EVec_re[i+m][2*j]);
		}
	}

	return isOK;
}


inline int EigHermitian(const Matrix<double>& AA, Matrix<double>& EEVec, Matrix<double>& EEVal)
{
	cerr << "ERROR: Must not arrive here! EigHermitian() is only defined for complex matrices.\n";
	assert(0);
	return 0;
}


inline int EigHermitian(const Matrix<float>& AA, Matrix<float>& EEVec, Matrix<float>& EEVal)
{
	cerr << "ERROR: Must not arrive here! EigHermitian() is only defined for complex matrices.\n";
	assert(0);
	return 0;
}


template <typename T> int EigComplex(const Matrix<T>& AA, Matrix<T>& EEVec, Matrix<T>& EEVal){
	static Matrix<T> TT,QQ;
	int isOK;
	int m=NumRows(AA);
	int n=NumCols(AA);
	
	QQ.Resize(m,n);
	TT=AA;
	EEVec.Resize(m,n);
	EEVal.Resize(m,n);
	
	// compute Schur form: A = Q*T*Q^T
	isOK=SchurComplex(AA,TT,QQ);
	if(isOK==0) return isOK;
	
	// extract eigenvalues
	Vector<T>   evals;
	isOK=Matrix<T>::mes_zschur_evals(TT,evals);
	if(isOK==0) return isOK;
	
	EEVal.Resize(m,n); EEVal=T();
	EEVal.Resize(m,n); EEVal=T();
	int i;
	for(i=1;i<=m;i++){
		EEVal[i][i]=evals[i];
	}
	
	// Q not needed for eiegenvalues
	EEVec.Resize(m,n);
	isOK=Matrix<T>::mes_zschur_vecs(TT,QQ,EEVec);
	// k'th eigenvector is k'th column of (X_re + i*X_im)
	
	return isOK;
}


#ifndef DISABLE_COMPLEX
inline int EigComplex(const Matrix<complex<float> >& A, Matrix<complex<float> >& EVec, Matrix<complex<float> >& EVal){
	static Matrix<complex<double> > AA,EEVec,EEVal;
	AA=A.ToCDouble();
	int isOK=EigComplex(AA,EEVec,EEVal);
	EVec.FromCDouble(EEVec);
	EVal.FromCDouble(EEVal);
	return isOK;
}
#endif


// ------------------------------------------------------------------------
//Matrix utilities
// ------------------------------------------------------------------------
template <typename T> Matrix<T> Matrix<T>::Randn(int m_, int n_){
	Matrix<T> C;
	if(NumRows(C)!=m_ || NumCols(C)!=n_) C.Resize(m_,n_);
	int i,j;
	for(i=1;i<=m_;i++) for(j=1;j<=n_;j++) {
		C[i][j]=T(GaussianRandom());
	}
	return C;
}


template <typename T> Matrix<T> Matrix<T>::CRandn(int m_, int n_){
	Matrix<T> C;
	if(NumRows(C)!=m_ || NumCols(C)!=n_) C.Resize(m_,n_);
	int i;
	for(i=1;i<=m_;i++) {
		C[i]=Vector<T>::CRandn(n_);
	}
	return C;
}


template <typename T> Matrix<T> Matrix<T>::Eye(int m_,int n_){
	Matrix<T> C(m_,n_);
	C.SetIdentity();
	return C;
}


template <typename T> Matrix<T> Matrix<T>::Zeros(int m_,int n_){
	Matrix<T> C(m_,n_);
	return C;
}


// ------------------------------------------------------------------------
// From modification of the Numerical Recipes in C
// ------------------------------------------------------------------------
#ifdef USE_NRC_CODE
#ifdef NRC_FILE_NAME
#include NRC_FILE_NAME
#else
#include "./nrc_Matrix.h"
#endif
#endif

// ------------------------------------------------------------------------
// From modification of the Meschach library
// ------------------------------------------------------------------------
#include "./mes_Matrix.h"

// ------------------------------------------------------------------------
// Matrix3
// ------------------------------------------------------------------------
#include "./Matrix3.h"
#include "./MatrixUtil.h"

#undef local_max
#undef local_min
#pragma warning(default: 4244)

#endif

