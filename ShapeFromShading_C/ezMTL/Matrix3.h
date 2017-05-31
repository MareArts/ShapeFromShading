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
// Filename: Matrix3.h
///////////////////////////////////////////////////////////////////////////////
#ifndef _MATRIX3_H_
#define _MATRIX3_H_
#include "./Matrix.h"

template <typename T> class Vector;
template <typename T> class Matrix;
template <typename T> class Matrix3;


// ------------------------------------------------------------------------
// Matrix3 declaration
// ------------------------------------------------------------------------
template<typename T>
class Matrix3 {
private:
	Vector<Matrix<T> > mat3;
		
public:
	Matrix3(int l_=0,int m_=0,int n_=0,const T& v_=T());
	Matrix3(const Matrix<T>& v_,int l_,int m_,int n_);
	virtual ~Matrix3();
#ifdef ALLOW_ACCESS_BY_BRACKET
	Matrix<T>&  operator[](int i);
	const Matrix<T>&  operator[](int i) const;
#endif
	Matrix<T>&  operator()(int i);
	const Matrix<T>&  operator()(int i) const;
	T&  operator()(int i,int j,int k);
	const T&  operator()(int i,int j,int k) const;
	int Print(ostream& s, const char *name=0);
	int Print(ostream& s, string& name);
	int Read(istream& s, string& name);
	int Read(istream& s, const char *name=0);
	int Resize(int new_l, int new_m, int new_n);
	Vector<int> Size();
	int Size(int d);
};

typedef Matrix3<float> FMatrix3;
typedef Matrix3<double> DMatrix3;


// ------------------------------------------------------------------------
// Matrix3 implementation
// ------------------------------------------------------------------------
template<typename T> Matrix3<T>::Matrix3(int l_,int m_,int n_,const T& v_){
	int i;
	mat3.Resize(l_);
	for(i=1;i<=l_;i++) mat3[i].Resize(m_,n_);
	for(i=1;i<=l_;i++) mat3[i]=v_;
}


template<typename T> Matrix3<T>::Matrix3(const Matrix<T>& v_,int l_,int m_,int n_){
	int i;
	mat3.Resize(l_);
	for(i=1;i<=l_;i++) mat3[i]=v_;
	for(i=1;i<=l_;i++) mat3[i].Resize(m_,n_);
}


template<typename T> Matrix3<T>::~Matrix3() {
}


#ifdef ALLOW_ACCESS_BY_BRACKET
template<typename T> inline Matrix<T>&  Matrix3<T>::operator[](int i){
	return mat3[i]; 
}


template<typename T> inline const Matrix<T>&  Matrix3<T>::operator[](int i) const {
	return mat3[i]; 
}
#endif


template<typename T> inline Matrix<T>&  Matrix3<T>::operator()(int i){
	return mat3(i);
}


template<typename T> inline const Matrix<T>&  Matrix3<T>::operator()(int i) const {
	return mat3(i);
}


template<typename T> inline T&  Matrix3<T>::operator()(int i,int j,int k){
	return (mat3(i))(j,k);
}


template<typename T> inline const T&  Matrix3<T>::operator()(int i,int j,int k) const {
	return (mat3(i))(j,k);
}


template<typename T> int Matrix3<T>::Print(ostream& s, const char *name) {
	if(name != 0 && name[0] != 0) s << name << " ";
	s << 3 << " " << mat3.Size() << " " << NumRows(mat3(1)) << " " << NumCols(mat3(1)) << endl;
	for(int i=1;i<=mat3.Size();i++){
		s << mat3[i] << endl;
	}
	return 1;
}


template<typename T> int Matrix3<T>::Print(ostream& s, string& name) {
	return Print(s,name.c_str());
}


template<typename T> int Matrix3<T>::Read(istream& s, const char *name){
	string dummy;
	s >> dummy;
	if(name != 0 && name[0] != 0){
		assert(name==dummy);
	}
	int itmp, ll, mm, nn;
	s >> itmp >> ll >> mm >> nn;
	assert(itmp==3);
	Resize(ll,mm,nn);
	for(int i=1;i<=mat3.Size();i++){
		s >> mat3[i];
	}
	return 1;
}


template<typename T> int Matrix3<T>::Read(istream& s, string& name){
	return Read(s,name.c_str());
}


template<typename T> int Matrix3<T>::Resize(int new_l, int new_m, int new_n) {
	// must preserve data
	int i;
	int l=mat3.Size();
	for(i=1;i<=local_min(l,new_l);i++) mat3[i].Resize(new_m,new_n);
	mat3.Resize(new_l);
	for(i=l+1;i<=new_l;i++) mat3[i].Resize(new_m,new_n);
	return 1;
}


template<typename T> inline Vector<int> Matrix3<T>::Size() {
	Vector<int> C(3,ROW_VECTOR);
	C(1)=mat3.Size(); C(2)=Size(mat3(1),1); C(3)=Size(mat3(1),2);
	return C;
}


template<typename T> int Matrix3<T>::Size(int d) {
	assert(d==1 || d==2 || d==3);
	int dim=0;
	if(d==1) dim=mat3.Size();
	else if(d==2) dim=Size(mat3(1),1);
	else if(d==3) dim=Size(mat3(1),2);
	return dim;
}


#endif
