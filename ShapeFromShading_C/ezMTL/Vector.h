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
// Filename: Vector.h
// Implementation:
//   Template class for vector algebra
//   Vector index starts from 1.
//   Access by A(i) rather than A[i].
// Note: 
//   Access by A(i) loses nothing but may gain something at least in some cases.
//   All return values are of the double type to prevent overflow or underflow.
// If you find any bugs, please send me an e-mail.
///////////////////////////////////////////////////////////////////////////////
#ifndef _VECTOR_H_
#define _VECTOR_H_
#include "./MatrixConfig.h"
#include "./MathUtil.h"
#include "./Random.h"

typedef complex<float> FComplex;
typedef complex<double> DComplex;
template <typename T> class Vector;
template <typename T> class Matrix;
enum VectorType { COL_VECTOR, ROW_VECTOR };


// Template functions with default arguments
#ifdef __GNUC__
template <typename T> T Max(const Vector<T>& A,int* idx=0);
template <typename T> T Min(const Vector<T>& A,int* idx=0);
template <typename T> Matrix<T> ToMatrix(const Vector<T>& A, int m=1);
template <typename T> double Norm(const Vector<T>& A, int p=2);
#else
template <typename T> T Max(const Vector<T>& A,int* idx);
template <typename T> T Min(const Vector<T>& A,int* idx);
template <typename T> Matrix<T> ToMatrix(const Vector<T>& A, int m);
template <typename T> double Norm(const Vector<T>& A, int p);
#endif

// friend template functions
template <typename T> Vector<T> mtl_plus(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<T> mtl_minus(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<T> mtl_plus(const Vector<T>& A, const double b);
template <typename T> Vector<T> mtl_minus(const Vector<T>& A, const double b);
template <typename T> Vector<T> mtl_times(const Vector<T>& A, const double b);
template <typename T> Vector<T> mtl_divide(const Vector<T>& A, const double b);
template <typename T> Vector<T> mtl_plus(const double b, const Vector<T>& A);
template <typename T> Vector<T> mtl_minus(const double b, const Vector<T>& A);
template <typename T> Vector<T> mtl_times(const double b, const Vector<T>& A);
template <typename T> Vector<T> mtl_divide(const double b, const Vector<T>& A);
template <typename T> Vector<T> mtl_times(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<T> mtl_divide(const Vector<T>& A, const Vector<T>& B);
template <typename T> ostream& mtl_ostream(ostream& s, const Vector<T>& A);
template <typename T> istream& mtl_istream(istream& s, Vector<T>& A);

#ifndef DISABLE_COMPLEX
template <typename T> Vector<T> mtl_plus(const Vector<T>& A, const complex<T> b);
template <typename T> Vector<T> mtl_minus(const Vector<T>& A, const complex<T> b);
template <typename T> Vector<T> mtl_times(const Vector<T>& A, const complex<T> b);
template <typename T> Vector<T> mtl_divide(const Vector<T>& A, const complex<T> b);
template <typename T> Vector<T> mtl_plus(const complex<T> b, const Vector<T>& A);
template <typename T> Vector<T> mtl_minus(const complex<T> b, const Vector<T>& A);
template <typename T> Vector<T> mtl_times(const complex<T> b, const Vector<T>& A);
template <typename T> Vector<T> mtl_divide(const complex<T> b, const Vector<T>& A);
#endif

template <typename T> int Length(const Vector<T>& A);
template <typename T> T Sum(const Vector<T>& A);
template <typename T> Vector<T> Transpose(const Vector<T>& A);
template <typename T> Vector<T> ArrayTranspose(const Vector<T>& A);
template <typename T> Matrix<T> OuterProd(const Vector<T>& A, const Vector<T>& B);
template <typename T> T InnerProd(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<T> Cross(const Vector<T>& A, const Vector<T>& B);
template <typename T> T Dot(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<T> ArrayMultiply(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<T> ArrayDivide(const Vector<T>& A, const Vector<T>& B);
template <typename T> Vector<int> Int(const Vector<T>& A);
template <typename T> Vector<float> Float(const Vector<T>& A);
template <typename T> Vector<double> Double(const Vector<T>& A);
template <typename T> Vector<T> Abs(const Vector<T>& A);
template <typename T> Vector<T> Sign(const Vector<T>& A);
template <typename T> Vector<T> Pow(const Vector<T>& A, double b);
template <typename T> Vector<T> Pow(const Vector<T>& A, int b);
template <typename T> Vector<T> Exp(const Vector<T>& A);
template <typename T> Vector<T> Log(const Vector<T>& A);
template <typename T> Vector<T> Digamma(const Vector<T>& A);
template <typename T> Vector<T> Gammaln(const Vector<T>& A);
template <typename T> Matrix<T> Diag(const Vector<T> A);
template <typename T> T Mean(const Vector<T>& A);
template <typename T> Vector<T> Cumsum(const Vector<T>& A);
template <typename T> Vector<T> Cumprod(const Vector<T>& A);
template <typename T> string Num2str(const Vector<T>& A);
template <typename T> int Sort(const Vector<T>& A, Vector<T>& v, Vector<int>& orderA, int order);
template <typename T> Vector<T> FFT(const Vector<T>& A,int N);
template <typename T> Vector<T> IFFT(const Vector<T>& A,int N);
template <typename T> Vector<T> Merge(const Vector<T>& A,const Vector<T>& B);



// ------------------------------------------------------------------------
// Vector declaration
// ------------------------------------------------------------------------
template <typename T>
class Vector {
private:
	// I chose "vector" in STL as a storage because the "vector" preserves
	// data after resizing. I waste one element but it made implementation easy.
	int n;
	VectorType type;
	vector<T> vec;
public:
	friend class Matrix<T>;
	Vector(int n_=0, VectorType type_=(VectorType)COL_VECTOR, const T& v_=T());
	Vector(int n_, VectorType type_, const char* v_);
	Vector(const Vector<T>& v_, int n_, VectorType type_);
	//Vector(const Vector<T>& A);
	virtual ~Vector();

	/////////////////////////////
	// member operator overloading
	/////////////////////////////
	T&			operator()(int i);
	const T&	operator()(int i) const;
	Vector<T>	operator+();
	Vector<T>	operator-();
	Vector<T>& operator+=(const Vector<T>& B);
	Vector<T>& operator-=(const Vector<T>& B);
	Vector<T>& operator+=(const double b);
	Vector<T>& operator-=(const double b);
	Vector<T>& operator*=(const double b);
	Vector<T>& operator/=(const double b);
#ifndef DISABLE_COMPLEX
	Vector<T>& operator+=(const complex<T> b);
	Vector<T>& operator-=(const complex<T> b);
	Vector<T>& operator*=(const complex<T> b);
	Vector<T>& operator/=(const complex<T> b);
#endif
	//Vector<T>& operator=(const Vector<T>& A);
	Vector<T>& operator=(const T& a);
#ifdef ALLOW_ACCESS_BY_BRACKET
	// The [] operator does not check if the arguments are in the valid range.
	// I recommend to use it only in implementation of the Vector class.
	T&			operator[](int i);
	const T&	operator[](int i) const;
#endif

	/////////////////////////////
	// friend operator overloading
	/////////////////////////////
	// GCC-2.95-2 did not allow friend operator overloading definition outside class declaration.
	friend Vector<T> operator+(const Vector<T>& A, const Vector<T>& B) {
		return mtl_plus(A,B);
	}
	friend Vector<T> operator-(const Vector<T>& A, const Vector<T>& B) {
		return mtl_minus(A,B);
	}
	friend Vector<T> operator+(const Vector<T>& A, const double b) {
		return mtl_plus(A,b);
	}
	friend Vector<T> operator-(const Vector<T>& A, const double b) {
		return mtl_minus(A,b);
	}
	friend Vector<T> operator*(const Vector<T>& A, const double b) {
		return mtl_times(A,b);
	}
	friend Vector<T> operator/(const Vector<T>& A, const double b) {
		return mtl_divide(A,b);
	}
	friend Vector<T> operator+(const double b, const Vector<T>& A) {
		return mtl_plus(b,A);
	}
	friend Vector<T> operator-(const double b, const Vector<T>& A) {
		return mtl_minus(b,A);
	}
	friend Vector<T> operator*(const double b, const Vector<T>& A) {
		return mtl_times(b,A);
	}
	friend Vector<T> operator/(const double b, const Vector<T>& A) {
		return mtl_divide(b,A);
	}
	friend ostream& operator<<(ostream& s, const Vector<T>& A) {
		return mtl_ostream(s,A);
	}
	friend ostream& operator<<(ostream& s, Vector<T>& A) {
		return mtl_ostream(s,A);
	}
	friend istream& operator>>(istream& s, Vector<T>& A) {
		return mtl_istream(s,A);
	}

#ifndef DISABLE_COMPLEX
	friend Vector<T> operator+(const Vector<T>& A, const complex<T> b) {
		return mtl_plus(A,b);
	}
	friend Vector<T> operator-(const Vector<T>& A, const complex<T> b) {
		return mtl_minus(A,b);
	}
	friend Vector<T> operator*(const Vector<T>& A, const complex<T> b) {
		return mtl_times(A,b);
	}
	friend Vector<T> operator/(const Vector<T>& A, const complex<T> b) {
		return mtl_divide(A,b);
	}
	friend Vector<T> operator+(const complex<T> b, const Vector<T>& A) {
		return mtl_plus(b,A);
	}
	friend Vector<T> operator-(const complex<T> b, const Vector<T>& A) {
		return mtl_minus(b,A);
	}
	friend Vector<T> operator*(const complex<T> b, const Vector<T>& A) {
		return mtl_times(b,A);
	}
	friend Vector<T> operator/(const complex<T> b, const Vector<T>& A) {
		return mtl_divide(b,A);
	}
#endif

	/////////////////////////////
	// friend functions
	/////////////////////////////
#ifdef __GNUC__
	friend Vector<T> mtl_plus<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> mtl_minus<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> mtl_plus<T>(const Vector<T>& A, const double b);
	friend Vector<T> mtl_minus<T>(const Vector<T>& A, const double b);
	friend Vector<T> mtl_times<T>(const Vector<T>& A, const double b);
	friend Vector<T> mtl_divide<T>(const Vector<T>& A, const double b);
	friend Vector<T> mtl_plus<T>(const double b, const Vector<T>& A);
	friend Vector<T> mtl_minus<T>(const double b, const Vector<T>& A);
	friend Vector<T> mtl_times<T>(const double b, const Vector<T>& A);
	friend Vector<T> mtl_divide<T>(const double b, const Vector<T>& A);
	friend Vector<T> mtl_times<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> mtl_divide<T>(const Vector<T>& A, const Vector<T>& B);
	friend ostream& mtl_ostream<T>(ostream& s, const Vector<T>& A);
	friend istream& mtl_istream<T>(istream& s, Vector<T>& A);

#ifndef DISABLE_COMPLEX
	friend Vector<T> mtl_plus<T>(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_minus<T>(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_times<T>(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_divide<T>(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_plus<T>(const complex<T> b, const Vector<T>& A);
	friend Vector<T> mtl_minus<T>(const complex<T> b, const Vector<T>& A);
	friend Vector<T> mtl_times<T>(const complex<T> b, const Vector<T>& A);
	friend Vector<T> mtl_divide<T>(const complex<T> b, const Vector<T>& A);
#endif

	friend int Length<T>(const Vector<T>& A);
	friend T Sum<T>(const Vector<T>& A);
	friend Vector<T> Transpose<T>(const Vector<T>& A);
	friend Vector<T> ArrayTranspose<T>(const Vector<T>& A);
	friend Matrix<T> OuterProd<T>(const Vector<T>& A, const Vector<T>& B);
	friend T InnerProd<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> Cross<T>(const Vector<T>& A, const Vector<T>& B);
	friend T Dot<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> ArrayMultiply<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> ArrayDivide<T>(const Vector<T>& A, const Vector<T>& B);
	friend Vector<int> Int<T>(const Vector<T>& A);
	friend Vector<float> Float<T>(const Vector<T>& A);
	friend Vector<double> Double<T>(const Vector<T>& A);
	friend Vector<T> Abs<T>(const Vector<T>& A);
	friend Vector<T> Sign<T>(const Vector<T>& A);
	friend Vector<T> Pow<T>(const Vector<T>& A, double b);
	friend Vector<T> Pow<T>(const Vector<T>& A, int b);
	friend Vector<T> Exp<T>(const Vector<T>& A);
	friend Vector<T> Log<T>(const Vector<T>& A);
	friend Vector<T> Digamma<T>(const Vector<T>& A);
	friend Vector<T> Gammaln<T>(const Vector<T>& A);
	friend Matrix<T> Diag<T>(const Vector<T> A);
	friend T Mean<T>(const Vector<T>& A);
	friend Vector<T> Cumsum<T>(const Vector<T>& A);
	friend Vector<T> Cumprod<T>(const Vector<T>& A);
	friend string Num2str<T>(const Vector<T>& A);
	friend int Sort<T>(const Vector<T>& A, Vector<T>& v, Vector<int>& orderA, int order);
	friend Vector<T> FFT<T>(const Vector<T>& A,int N);
	friend Vector<T> IFFT<T>(const Vector<T>& A,int N);
	friend Vector<T> Merge<T>(const Vector<T>& A,const Vector<T>& B);
#else
	friend Vector<T> mtl_plus(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> mtl_minus(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> mtl_plus(const Vector<T>& A, const double b);
	friend Vector<T> mtl_minus(const Vector<T>& A, const double b);
	friend Vector<T> mtl_times(const Vector<T>& A, const double b);
	friend Vector<T> mtl_divide(const Vector<T>& A, const double b);
	friend Vector<T> mtl_plus(const double b, const Vector<T>& A);
	friend Vector<T> mtl_minus(const double b, const Vector<T>& A);
	friend Vector<T> mtl_times(const double b, const Vector<T>& A);
	friend Vector<T> mtl_divide(const double b, const Vector<T>& A);
	friend Vector<T> mtl_times(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> mtl_divide(const Vector<T>& A, const Vector<T>& B);
	friend ostream& mtl_ostream(ostream& s, const Vector<T>& A);
	friend istream& mtl_istream(istream& s, Vector<T>& A);

#ifndef DISABLE_COMPLEX
	friend Vector<T> mtl_plus(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_minus(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_times(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_divide(const Vector<T>& A, const complex<T> b);
	friend Vector<T> mtl_plus(const complex<T> b, const Vector<T>& A);
	friend Vector<T> mtl_minus(const complex<T> b, const Vector<T>& A);
	friend Vector<T> mtl_times(const complex<T> b, const Vector<T>& A);
	friend Vector<T> mtl_divide(const complex<T> b, const Vector<T>& A);
#endif

	friend int Length(const Vector<T>& A);
	friend T Sum(const Vector<T>& A);
	friend Vector<T> Transpose(const Vector<T>& A);
	friend Vector<T> ArrayTranspose(const Vector<T>& A);
	friend Matrix<T> OuterProd(const Vector<T>& A, const Vector<T>& B);
	friend T InnerProd(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> Cross(const Vector<T>& A, const Vector<T>& B);
	friend T Dot(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> ArrayMultiply(const Vector<T>& A, const Vector<T>& B);
	friend Vector<T> ArrayDivide(const Vector<T>& A, const Vector<T>& B);
	friend Vector<int> Int(const Vector<T>& A);
	friend Vector<float> Float(const Vector<T>& A);
	friend Vector<double> Double(const Vector<T>& A);
	friend Vector<T> Abs(const Vector<T>& A);
	friend Vector<T> Sign(const Vector<T>& A);
	friend Vector<T> Pow(const Vector<T>& A, double b);
	friend Vector<T> Pow(const Vector<T>& A, int b);
	friend Vector<T> Exp(const Vector<T>& A);
	friend Vector<T> Log(const Vector<T>& A);
	friend Vector<T> Digamma(const Vector<T>& A);
	friend Vector<T> Gammaln(const Vector<T>& A);
	friend Matrix<T> Diag(const Vector<T> A);
	friend T Mean(const Vector<T>& A);
	friend Vector<T> Cumsum(const Vector<T>& A);
	friend Vector<T> Cumprod(const Vector<T>& A);
	friend string Num2str(const Vector<T>& A);
	friend int Sort(const Vector<T>& A, Vector<T>& v, Vector<int>& orderA, int order);
	friend Vector<T> FFT(const Vector<T>& A,int N);
	friend Vector<T> IFFT(const Vector<T>& A,int N);
	friend Vector<T> Merge(const Vector<T>& A,const Vector<T>& B);

#endif
	
	/////////////////////////////
	// member functions
	/////////////////////////////
	int Size() const;
	VectorType Type() const;
	int SetType(VectorType type_);
	Vector<T> t() const;
	Vector<T> ct() const;
	int Resize(int new_n);
	int Add(const T& a);
	int Erase(int c);
	void Clear();
	bool Empty() const;
	T Sum_() const;
	T Max_(int* idx=0) const;
	T Min_(int* idx=0) const;
	int Print(ostream& s, const char *name=0);
	int Print(ostream& s, string& name);
	int Read(istream& s, string& name);
	int Read(istream& s, const char *name=0);
	Vector<T> SubVector(int n1, int n2);
	int Save(const string& fname);
	int SaveMatlab(const string& fname);
	void SetPermIdentity();

	int IsOne(double tolerance=1e-6) const;
	int IsZero(double tolerance=1e-6) const;
	int IsUnit(double tolerance=1e-6) const;

	void FromInt(const Vector<int>& A);
	void FromFloat(const Vector<float>& A);
	void FromDouble(const Vector<double>& A);
	Vector<int> ToInt();
	Vector<float> ToFloat();
	Vector<double> ToDouble();

public:
	/////////////////////////////
	// exported static functions
	/////////////////////////////
	static Vector<T> Rand(int m);
	static Vector<T> Randn(int m);
	static Vector<T> CRandn(int m);

	/////////////////////////////
	// local static functions
	/////////////////////////////
	static double normp(const Vector<T>& A,int p);
	static int FFTBase(const Vector<T>& A,int N,int isign,Vector<T>& C);
	static int FFTCore(Vector<T>& data,int nn,int isign);

};
typedef Vector<float> FVector;
typedef Vector<double> DVector;
typedef Vector<FComplex> CFVector;
typedef Vector<DComplex> CDVector;

// ------------------------------------------------------------------------
// Vector implementation
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
template <typename T> Vector<T>::Vector(int n_,VectorType type_,const T& v_) : n(n_),type(type_) {
	vec.resize(n_+1);
	for(int i=1;i<=n_;i++) vec[i]=v_;
}


template <typename T> Vector<T>::Vector(const Vector<T>& v_,int n_,VectorType type_) : n(n_),type(type_) {
	vec.resize(n_+1);
	int i;
	for(i=1;i<=n_;i++) {
		if(i<=v_.Size()) vec[i]=v_[i];
		else vec[i]=T();
	}
}


template <typename T> Vector<T>::Vector(int n_,VectorType type_,const char* v_) : n(n_),type(type_) {
	vec.resize(n_+1);
	Vector<string> num=mtl_split(v_);
	for(int i=1;i<=n_;i++) {
		if(i<=num.Size()) vec[i]=T(atof(num[i].c_str()));
		else vec[i]=T();
	}
}


//template <typename T> Vector<T>::Vector(const Vector<T>& A) : vec(A.vec),n(A.n),type(A.type) {
//}


template <typename T> Vector<T>::~Vector() {
}



/////////////////////////////
// operator overloading
/////////////////////////////
//template <typename T> Vector<T>& Vector<T>::operator=(const Vector<T>& A) {
//	if(this == &A) return *this;
//	if(n != A.n) { vec.resize(A.n); n=A.n; }
//	vec=A.vec; type=A.type; return *this;
//}
template <typename T> Vector<T>& Vector<T>::operator=(const T& a) {
	for(int i=1; i<=n; ++i) vec[i] = a; return *this;
}


#ifdef ALLOW_ACCESS_BY_BRACKET
#if defined(_DEBUG)
template <typename T> T& Vector<T>::operator[](int i) {
	if(i<=0 || i>n){
		cerr << "ERROR: index invalid " << i << endl;
		assert(0);
	}
	return vec[i];
}


template <typename T> const T& Vector<T>::operator[](int i) const {
	if(i<=0 || i>n){
		cerr << "ERROR: index invalid " << i << endl;
		assert(0);
	}
	return vec[i];
}


#else // !_DEBUG
template <typename T> inline T& Vector<T>::operator[](int i) {
	return vec[i];
}


template <typename T> inline const T& Vector<T>::operator[](int i) const {
	return vec[i];
}
#endif // _DEBUG
#endif // ALLOW_ACCESS_BY_BRACKET


#if defined(_DEBUG)
template <typename T> const T& Vector<T>::operator()(int i) const {
	if(i<=0 || i>n){
		cerr << "ERROR: index invalid " << i << endl;
		assert(0);
	}
	return vec[i];
}


template <typename T> T& Vector<T>::operator()(int i) {
	if(i<=0 || i>n){
		cerr << "ERROR: index invalid " << i << endl;
		assert(0);
	}
	return vec[i];
}


#else // !_DEBUG
template <typename T> inline const T& Vector<T>::operator()(int i) const {
	return vec[i];
}


template <typename T> inline T& Vector<T>::operator()(int i) {
	return vec[i];
}

#endif // _DEBUG


template <typename T> inline Vector<T> Vector<T>::operator+() {
	return (*this);
}


template <typename T> Vector<T> Vector<T>::operator-() {
	Vector<T> C=(*this); C*=(-1); return C;
}


template <typename T> Vector<T>& Vector<T>::operator+=(const Vector<T>& B) {
	assert(n==B.n);
	assert(type==B.type);
	for(int i=1;i<=n;i++) vec[i]+=B.vec[i]; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator-=(const Vector<T>& B) {
	assert(n==B.n);
	assert(type==B.type);
	for(int i=1;i<=n;i++) vec[i]-=B.vec[i]; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator+=(const double b) {
	for(int i=1;i<=n;i++) vec[i]+=b; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator-=(const double b) {
	for(int i=1;i<=n;i++) vec[i]-=b; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator*=(const double b) {
	for(int i=1;i<=n;i++) vec[i]*=b; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator/=(const double b) {
	for(int i=1;i<=n;i++) vec[i]/=b; return (*this);
}


#ifndef DISABLE_COMPLEX
template <typename T> Vector<T>& Vector<T>::operator+=(const complex<T> b) {
	for(int i=1;i<=n;i++) vec[i]+=b; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator-=(const complex<T> b) {
	for(int i=1;i<=n;i++) vec[i]-=b; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator*=(const complex<T> b) {
	for(int i=1;i<=n;i++) vec[i]*=b; return (*this);
}


template <typename T> Vector<T>& Vector<T>::operator/=(const complex<T> b) {
	for(int i=1;i<=n;i++) vec[i]/=b; return (*this);
}
#endif


/////////////////////////////
// member functions
/////////////////////////////
template <typename T> inline int Vector<T>::Size() const { 
	return n;
}


template <typename T> inline VectorType Vector<T>::Type() const { 
	return type;
}


template <typename T> inline int Vector<T>::SetType(VectorType type_) {
	type=type_;
	return 1;
}


template <typename T> Vector<T> Vector<T>::t() const {
	Vector<T> C(Size(),(Type()==ROW_VECTOR)?COL_VECTOR:ROW_VECTOR,T());
	for(int i=1;i<=Size();i++) C[i]=((*this)[i]);
	return C;
}


template <typename T> Vector<T> Vector<T>::ct() const {
	Vector<T> C(Size(),(Type()==ROW_VECTOR)?COL_VECTOR:ROW_VECTOR,T());
	for(int i=1;i<=Size();i++) C[i]=conj((*this)[i]);
	return C;
}


template <typename T> inline int Vector<T>::Resize(int new_n) { 
	vec.resize(new_n+1);
	n=new_n;
	return 1;
}


template <typename T> int Vector<T>::Add(const T& a){ 
	int n0=n; 
	Resize(n0+1); 
	vec[n0+1]=a; 
	n=n0+1; 
	return n; 
}


template <typename T> int Vector<T>::Erase(int c){ 
	vec.erase(vec.begin()+c); 
	n--;
	return 1;
}


template <typename T> inline void Vector<T>::Clear(){ 
	vec.clear(); 
}


template <typename T> inline bool Vector<T>::Empty() const{ 
	return (Size()==0);
}


template <typename T> T Vector<T>::Sum_() const {
	T sum=0; 
	for(int i=1;i<=n;i++) sum+=vec[i]; 
	return sum; 
}


template <typename T> T Min(const Vector<T>& A,int* idx=0){
	return A.Min_(idx); 
}


template <typename T> T Vector<T>::Max_(int* idx) const {
	T maxA=-mtl_numeric_limits<T>::max();
	if(idx) *idx=0;
	for(int i=1;i<=n;i++){
		if(maxA<vec[i]){
			maxA=vec[i];
			if(idx) *idx=i;
		}
	}
	return T(maxA);
}


template <typename T> T Vector<T>::Min_(int* idx) const {
	T minA=mtl_numeric_limits<T>::max();
	if(idx) *idx=0;
	for(int i=1;i<=n;i++){
		if(minA>vec[i]){
			minA=vec[i];
			if(idx) *idx=i;
		}
	}
	return T(minA);
}


template <typename T> int Vector<T>::Print(ostream& s, const char *name) {
	if(name != 0 && name[0] != 0) s << name << " ";
	s << 1 << " " << n << endl;
	//s << (*this) << endl;
	mtl_ostream(s,(*this));
	return 1;
}


template <typename T> int Vector<T>::Print(ostream& s, string& name) {
	return Print(s,name.c_str());
}


template <typename T> int Vector<T>::Read(istream& s, const char *name) {
	string dummy;
	int itmp, nn;
	s >> dummy;
	if(name != 0 && name[0] != 0){
		assert(name==dummy);
	}
	s >> itmp >> nn;
	assert(itmp==1);
	Resize(nn);
	s >> (*this);
	return 1;
}


template <typename T> int Vector<T>::Read(istream& s, string& name) {
	return Read(s,name.c_str());
}


template <typename T> Vector<T> Vector<T>::SubVector(int n1, int n2){
	assert(n2>=n1);
	Vector<T> C(n2-n1+1); for(int i=1;i<=n2-n1+1;i++) C[i]=vec[i-1+n1];
	return C;
}


template <typename T> int Vector<T>::Save(const string& fname){
	FILE* fp=fopen(fname.c_str(),"wb");
	if(!fp) return 0;
	fprintf(fp,"FMAT 1 %d\n",n);
	for(int j=1;j<=n;j++){
		if((j-1)%10!=0) fprintf(fp," ");
		fprintf(fp,"%f",vec[j]);
		if(j%10==0) fprintf(fp,"\n");
	}
	if(n%10!=0) fprintf(fp,"\n");
	fclose(fp);
	return n;
}


template <typename T> int Vector<T>::SaveMatlab(const string& fname){
	FILE* fp=fopen(fname.c_str(),"wb");
	if(!fp) return 0;
	for(int j=1;j<=n;j++){
		fprintf(fp,"%f ",vec[j]);
	}
	fclose(fp);
	return n;
}


template <typename T> int Vector<T>::IsOne(double tolerance) const {
	for (int i=1; i<=n; i++) if(Abs(vec[i]-1)>tolerance) return 0;
	return 1;
}


template <typename T> int Vector<T>::IsZero(double tolerance) const {
	for (int i=1; i<=n; i++) if(Abs(vec[i])>tolerance) return 0;
	return 1;
}


template <typename T> int Vector<T>::IsUnit(double tolerance) const {
	double norm=Norm(*this);
	if(Abs(norm-1)>tolerance) return 0;
	return 1;
}



/////////////////////////////
// basic friend functions
/////////////////////////////
template <typename T> Vector<T> mtl_plus(const Vector<T>& A, const Vector<T>& B) {
	assert(A.n==B.n);
	assert(A.type==B.type);
	Vector<T> C=A; C+=B; return C;
}


template <typename T> Vector<T> mtl_minus(const Vector<T>& A, const Vector<T>& B) {
	assert(A.n==B.n);
	assert(A.type==B.type);
	Vector<T> C=A; C-=B; return C;
}


template <typename T> Vector<T> mtl_plus(const Vector<T>& A, const double b) {
	Vector<T> C=A; C+=b; return C;
}


template <typename T> Vector<T> mtl_minus(const Vector<T>& A, const double b) {
	Vector<T> C=A; C-=b; return C;
}


template <typename T> Vector<T> mtl_times(const Vector<T>& A, const double b) {
	Vector<T> C=A; C*=b; return C;
}


template <typename T> Vector<T> mtl_divide(const Vector<T>& A, const double b) {
	Vector<T> C=A; C/=b; return C;
}


template <typename T> Vector<T> mtl_plus(const double b, const Vector<T>& A) {
	Vector<T> C=A; C+=b; return C;
}


template <typename T> Vector<T> mtl_minus(const double b, const Vector<T>& A) {
	Vector<T> C=-A; C+=b; return C;
}


template <typename T> Vector<T> mtl_times(const double b, const Vector<T>& A) {
	Vector<T> C=A; C*=b; return C;
}


template <typename T> Vector<T> mtl_divide(const double b, const Vector<T>& A) {
	Vector<T> C(A.n); for(int i=1;i<=A.n;i++) C.vec[i]=b/A.vec[i]; return C;
}


#ifndef DISABLE_COMPLEX
template <typename T> Vector<T> mtl_plus(const Vector<T>& A, const complex<T> b) {
	Vector<T> C=A; C+=b; return C;
}


template <typename T> Vector<T> mtl_minus(const Vector<T>& A, const complex<T> b) {
	Vector<T> C=A; C-=b; return C;
}


template <typename T> Vector<T> mtl_times(const Vector<T>& A, const complex<T> b) {
	Vector<T> C=A; C*=b; return C;
}


template <typename T> Vector<T> mtl_divide(const Vector<T>& A, const complex<T> b) {
	Vector<T> C=A; C/=b; return C;
}


template <typename T> Vector<T> mtl_plus(const complex<T> b, const Vector<T>& A) {
	Vector<T> C=A; C+=b; return C;
}


template <typename T> Vector<T> mtl_minus(const complex<T> b, const Vector<T>& A) {
	Vector<T> C=-A; C+=b; return C;
}


template <typename T> Vector<T> mtl_times(const complex<T> b, const Vector<T>& A) {
	Vector<T> C=A; C*=b; return C;
}


template <typename T> Vector<T> mtl_divide(const complex<T> b, const Vector<T>& A) {
	Vector<T> C(A.n); for(int i=1;i<=A.n;i++) C.vec[i]=b/A.vec[i]; return C;
}
#endif


template <typename T> Vector<T> mtl_times(const Vector<T>& A, const Vector<T>& B){
	assert(A.Size()==B.Size());
	assert(A.type==B.type);
	Vector<T> C=A;
	for(int i=1;i<=C.Size();i++) C.vec[i]*=B.vec[i];
	return C;
}


template <typename T> Vector<T> mtl_divide(const Vector<T>& A, const Vector<T>& B){
	assert(A.Size()==B.Size());
	assert(A.type==B.type);
	Vector<T> C=A;
	for(int i=1;i<=C.Size();i++) C.vec[i]/=B.vec[i];
	return C;
}


template <typename T> ostream& mtl_ostream(ostream& s, const Vector<T>& A) {
	s << "  ";
	for(int i=1; i<=A.n; i++) {
		s << Num2str(A[i]) << " ";
	}
	s << endl;
	return s;
}


template <typename T> istream& mtl_istream(istream& s, Vector<T>& A) {
	for (int i=1; i<=A.n; i++)
		s >> A[i];
	return s;
}


template <typename T> inline int Length(const Vector<T>& A) {
	return A.Size(); 
}


template <typename T> inline T Sum(const Vector<T>& A) {
	return A.Sum_(); 
}


template <typename T> inline T Max(const Vector<T>& A,int* idx=0){
	return A.Max_(idx); 
}


template <typename T> Matrix<T> Diag(const Vector<T> A){
	Matrix<T> C(A.Size(),A.Size());
	C=0;
	for(int i=1;i<=NumRows(C);i++) C[i][i]=A[i];
	return C;
}


template <typename T> double Norm(const Vector<T>& A, int p=2){
	if(p==mtl_numeric_limits<int>::max()) {
		double maxA=0;
		for(int i=1;i<=A.Size();i++) if(maxA<Abs(A[i])) maxA=Abs(A[i]);
		return maxA;		
	}
	else if(p==mtl_numeric_limits<int>::min()) {
		double minA=mtl_numeric_limits<double>::max();
		for(int i=1;i<=A.Size();i++) if(minA>Abs(A[i])) minA=Abs(A[i]);
		return minA;
	}
	else return Vector<T>::normp(A,p);
}


template <typename T> double Vector<T>::normp(const Vector<T>& A, int p) {
	double z=0, sum=0;
	int i;
	for(i=1;i<=A.Size();i++) if(z<Abs(A.vec[i])) z=Abs(A.vec[i]);
	if(z==0) return 0;
	for(i=1;i<=A.Size();i++) sum+=Ipow(Abs(A.vec[i])/z,p);
	return (z*pow(sum,1.0/p));
}


template <typename T> inline T Mean(const Vector<T>& A){
	return A.Sum_()/A.n;
}


template <typename T> Vector<T> Cumsum(const Vector<T>& A){
	Vector<T> C(A.n);
	T s=0;
	for(int i=1;i<=A.n;i++) { C(i)=s+A(i); s+=A(i); }
	return C;
}


template <typename T> Vector<T> Cumprod(const Vector<T>& A){
	Vector<T> C(A.n);
	T p=1;
	for(int i=1;i<=A.n;i++) { C(i)=p*A(i); p*=A(i); }
	return C;
}


template <typename T> string Num2str(const Vector<T>& A){
	string s;
	for(int i=1;i<=A.Size();i++){
		if(i>1) s += " ";
		s += Num2str(A[i]);
	}
	return s;
}


// order=1: ascending order
// order=-1: descending order
template <typename T> int VecSort(const Vector<T>& A, Vector<T>& v, Vector<int>& orderA, int order){
	int i;
	int n_=A.Size();
	v=A;
	orderA.Resize(n_);
	for(i=1;i<=n_;i++) orderA[i]=i;
	for(i=1;i<=n_;i++){
		for(int j=i+1;j<=n_;j++){
			if(order*v[i] > order*v[j]){
				T a=v[i];
				v[i]=v[j];
				v[j]=a;
				int b=orderA[i];
				orderA[i]=orderA[j];
				orderA[j]=b;
			}
		}
	}
	return 1;
}


template <typename T> inline Vector<T> Transpose(const Vector<T>& A){
	return A.ct();
}


template <typename T> Matrix<T> OuterProd(const Vector<T>& A, const Vector<T>& B){
	Matrix<T> C(A.Size(),B.Size());
	for(int i=1;i<=A.Size();i++) for(int j=1;j<=B.Size();j++) C[i][j]=A.vec[i]*B.vec[j];
	return C;
}


template <typename T> T InnerProd(const Vector<T>& A, const Vector<T>& B){
	assert(A.Size()==B.Size());
	T sum=T();
	for(int i=1;i<=A.Size();i++) sum+=A.vec[i]*B.vec[i];
	return sum;
}


template <typename T> Vector<T> Cross(const Vector<T>& A, const Vector<T>& B){
	assert(A.Size()==3 && B.Size()==3);
	Vector<T> C(3,A.Type());
	C[1]=A[2]*B[3]-A[3]*B[2];
	C[2]=A[3]*B[1]-A[1]*B[3];
	C[3]=A[1]*B[2]-A[2]*B[1];
	return C;
}


template <typename T> inline T Dot(const Vector<T>& A, const Vector<T>& B){
	return InnerProd(A,B);
}


template <typename T> inline Vector<T> ArrayMultiply(const Vector<T>& A, const Vector<T>& B){
	return mtl_times(A,B);
}


template <typename T> inline Vector<T> ArrayDivide(const Vector<T>& A, const Vector<T>& B){
	return mtl_divide(A,B);
}


// ToMatrix returns mxn matrix whose rows are repetitions of vector A
template <typename T> Matrix<T> ToMatrix(const Vector<T>& A, int m=1){
	assert(A.Type()==ROW_VECTOR);
	Matrix<T> C(m,A.Size());
	for(int i=1;i<=m;i++) C[i]=A;
	return C;
}


template <typename T> void Vector<T>::FromInt(const Vector<int>& A){
	Resize(A.Size());
	for(int i=1;i<=A.Size();i++) (*this)[i]=T(A[i]);
	return;
}


template <typename T> void Vector<T>::FromFloat(const Vector<float>& A){
	Resize(A.Size());
	for(int i=1;i<=A.Size();i++) (*this)[i]=T(A[i]);
	return;
}


template <typename T> void Vector<T>::FromDouble(const Vector<double>& A){
	Resize(A.Size());
	for(int i=1;i<=A.Size();i++) (*this)[i]=T(A[i]);
	return;
}


template <typename T> Vector<int> Vector<T>::ToInt(){
	Vector<int> C(Size());
	for(int i=1;i<=A.Size();i++) C[i]=(double)(*this)[i];
	return C;
}


template <typename T> Vector<float> Vector<T>::ToFloat(){
	Vector<float> C(Size());
	for(int i=1;i<=A.Size();i++) C[i]=(float)(*this)[i];
	return C;
}


template <typename T> Vector<double> Vector<T>::ToDouble(){
	Vector<double> C(Size());
	for(int i=1;i<=A.Size();i++) C[i]=(double)(*this)[i];
	return C;
}


template <typename T> Vector<int> Int(const Vector<T>& A){
	Vector<int> C(A.Size());
	for(int i=1;i<=A.Size();i++) C[i]=(int)A[i];
	return C;
}


template <typename T> Vector<float> Float(const Vector<T>& A){
	Vector<float> C(A.Size());
	for(int i=1;i<=A.Size();i++) C[i]=(float)A[i];
	return C;
}


template <typename T> Vector<double> Double(const Vector<T>& A){
	Vector<double> C(A.Size());
	for(int i=1;i<=A.Size();i++) C[i]=(double)A[i];
	return C;
}


template <typename T> Vector<T> Abs(const Vector<T>& A){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=T(Abs(A.vec[i]));
	return C;
}


template <typename T> Vector<T> Sign(const Vector<T>& A){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=T( ( (A.vec[i]>0) ? 1 : ((A.vec[i]<0) ? (-1) : 0) ) );
	return C;
}


template <typename T> Vector<T> Pow(const Vector<T>& A, double b){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=pow((double)A.vec[i],b);
	return C;
}


template <typename T> Vector<T> Pow(const Vector<T>& A, int b){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=Ipow(A.vec[i],b);
	return C;
}


template <typename T> Vector<T> Exp(const Vector<T>& A){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=exp(A.vec[i]);
	return C;
}


template <typename T> Vector<T> Log(const Vector<T>& A){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=Log(A.vec[i]);
	return C;
}


template <typename T> Vector<T> Digamma(const Vector<T>& A){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=Digamma(A.vec[i]);
	return C;
}


template <typename T> Vector<T> Gammaln(const Vector<T>& A){
	Vector<T> C(A.Size());
	for(int i=1;i<=A.Size();i++) C.vec[i]=Gammaln(A.vec[i]);
	return C;
}


template <typename T> Vector<T> FFT(const Vector<T>& A,int N){
	Vector<T> C;
	Vector<T>::FFTBase(A,N,1,C);
	return C;
}


template <typename T> Vector<T> IFFT(const Vector<T>& A,int N){
	Vector<T> C;
	Vector<T>::FFTBase(A,N,-1,C);
	return C;
}


template <typename T> int Vector<T>::FFTBase(const Vector<T>& A,int N,int isign,Vector<T>& C){
	int log2n=0; int pow2n=1;
	for(log2n=0;pow2n<N;log2n++) pow2n*=2;
	C.Resize(2*pow2n);
	for(int i=1;i<=A.Size();i++) C[i]=A[i];
	Vector<T>::FFTCore(C,pow2n,isign);
	if(isign<0) C /= T(pow2n);
	return pow2n;
}


template <typename T> inline void FFT_SWAP(T& a,T& b){
	T temp=a;
	a=b;
	b=temp;
}


template <typename T> int Vector<T>::FFTCore(Vector<T>& data,int nn,int isign)
{
	int n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			FFT_SWAP(data[j],data[i]);
			FFT_SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=2*mmax;
		theta=6.28318530717959/(isign*mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
	return 1;
}


template <typename T> Vector<T> Merge(const Vector<T>& A, const Vector<T>& B)
{
	Vector<T> C(A.Size()+B.Size());
	int i=1,n=1;
	for(i=1;i<=A.Size();i++,n++) C[n]=A[i];
	for(i=1;i<=B.Size();i++,n++) C[n]=B[i];
	return C;
}


// ------------------------------------------------------------------------
// Vector utilities
// ------------------------------------------------------------------------
template <typename T> Vector<T> Vector<T>::Rand(int m)
{
	Vector<T> C(m);
	for(int i=1;i<=m;i++) C[i]=(float)RANDF();
	return C;
}


template <typename T> Vector<T> Vector<T>::Randn(int m)
{
	Vector<T> C(m);
	for(int i=1;i<=m;i++) C[i]=T(GaussianRandom());
	return C;
}


template <typename T> Vector<T> Vector<T>::CRandn(int m)
{
	Vector<T> C(m);
	for(int i=1;i<=m;i++) {
		C[i]=T(GaussianRandom(),GaussianRandom());
	}
	return C;
}


template <typename T> void Vector<T>::SetPermIdentity()
{
	for(int i=1;i<=n;i++) (*this)[i]=T(i);
}


#ifndef DISABLE_COMPLEX
#ifdef __GNUC__
template <typename T> ostream& operator<<(ostream& s, const complex<T>& c) {
   s << "(" << c.real() << "," << c.imag() << ")";
   return s;
}
#endif
#endif

#undef local_max
#undef local_min


#endif

