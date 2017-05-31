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
// Filename: mes_math.h
// Revision:
//    1. Added Sqrtm(), Powm() and Logm().
// Note:
//    1. These routines do NOT support complex matrices yet.
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_MATH_H_
#define	_MES_MATH_H_

/* ------------- Mathematical functions -------------------- */

///////////////////////////////////////////////////////////////////////////////
// Expm(A)
// Computes matrix exponential using Pade approximation
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Expm(const Matrix<T>& A){
	if(A.IsComplex()){ // complex
		return Matrix<T>::ExpmComplex(A);
	}
	else if(A.IsSymmetric() && A.IsPD()){ // all eigenvalues are real positive
		Matrix<T> Atmp=A;
		Matrix<T> C;
		int q_out,j_out;
		int isOK=Matrix<T>::mes_exp(Atmp,Abs(mtl_numeric_limits<T>::min()),C,&q_out,&j_out);
		return C;
	}
	else{ // real general
		return Matrix<T>::ExpmReal(A);
	}
}


///////////////////////////////////////////////////////////////////////////////
// ExpmComplex(A)
// Computes matrix exponential for real general matrices using eigenvectors and eigenvalues
//    A = E * D * E^(-1)
//    expm(A) = E * expm(D) * E^(-1)
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Matrix<T>::ExpmReal(const Matrix<T>& A){
	Matrix<T> E,D;
	int isOK=Eig(A,E,D);
	if(isOK==0) {
		cerr << "ERROR: Failed in ExpmReal().\n";
		return Matrix<T>();
	}
	else if(isOK==1){ // real eigenvalues
		int i;
		for(i=1;i<=NumCols(D);i++){
			D[i][i]=exp(D[i][i]);
		}
		return mtl_mrdivide(E*D,E);
	}
	else{ // complex eigenvalues
		cerr<<"ERROR (ExpmReal): Input matrix has complex eigenvalues\n";
		cerr<<"      Use a complex matrix.\n";
		assert(0);
		return Matrix<T>();
	}
}


///////////////////////////////////////////////////////////////////////////////
// ExpmComplex(A)
// Computes matrix exponential for complex matrices using eigenvectors and eigenvalues
//    A = E * D * E^(-1)
//    expm(A) = E * expm(D) * E^(-1)
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Matrix<T>::ExpmComplex(const Matrix<T>& A){
	Matrix<T> E,D;
	int isOK=Eig(A,E,D);
	int i;
	for(i=1;i<=NumCols(D);i++){
		D[i][i]=exp(D[i][i]);
	}
	return mtl_mrdivide(E*D,E);
}


///////////////////////////////////////////////////////////////////////////////
// Logm(A)
// Computes matrix logarithm
// Note: Sqrtm is very fragile!
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Logm(const Matrix<T>& A){
	Matrix<T> ans;
	if(A.IsComplex()){ // complex
		ans=Matrix<T>::LogmComplex(A);
	}
	else{ // real general
		ans=Matrix<T>::LogmReal(A);
	}
	double tol=1e-6*Norm(A,mtl_numeric_limits<int>::max());
	if((A-Expm(ans)).IsZero(tol) == false){
		cerr<<"Warning (Logm): The answer may be inaccurate.\n";
		//assert(0);
	}
	return ans;
}


///////////////////////////////////////////////////////////////////////////////
// LogmComplex(A)
// Computes matrix exponential for complex matrices using eigenvectors and eigenvalues
//    A = E * D * E^(-1)
//    logm(A) = E * logm(D) * E^(-1)
// Note: Sqrtm is very fragile!
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Matrix<T>::LogmComplex(const Matrix<T>& A){
	Matrix<T> E,D;
	int isOK=Eig(A,E,D);
	int i;
	for(i=1;i<=NumCols(D);i++){
		D[i][i]=log(D[i][i]);
	}
	return mtl_mrdivide(E*D,E);
}


///////////////////////////////////////////////////////////////////////////////
// LogmReal(A)
// Computes matrix logarithm using the inverse scaling and squaring method
// introduced by Kenney and Laub.
//   Use the identity equation: log(A)=2^k log A^{1/2^k}
//   When k is large, A^{1/2^k} approaches to the identity matrix.
//   Then use the Taylor expansion of log(I-W) where W=I-A^{1/2^k}
//   log(I-W)= W - W^2/2 - W^3/3 - W^4/4 - ...
//   Final answer: log(A) = 2^k * log(I-W)
// Note: This implementation is very sensitive to the condition number of 
//   the input matrix.
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Matrix<T>::LogmReal(const Matrix<T>& A){
	int i,j,k;
	double eps=real(mtl_numeric_limits<T>::epsilon());
	// X <-- X^{1/2^k}
	Matrix<T> X=A;
	double powC=1.0;
	for(k=1;k<=16;k++){
		X=Sqrtm(X);
		powC *= 2;
		if(X.IsDiagonal(eps)){
			break;
		}
	}
	// X <-- I-X
	X *= -1;
	for(i=1;i<=X.m;i++) X[i][i] += 1;

	// compute log(I-X) = - X - X^2/2 - X^3/3 - ...
	Matrix<T> powX=X;
	Matrix<T> logX(X.m,X.n);
	for(k=1;k<=100;k++){
		double maxabs=0;
		for(i=1;i<=powX.m;i++){
			for(j=1;j<=powX.n;j++){
				if(maxabs<Abs(powX[i][j])) maxabs=Abs(powX[i][j]);
			}
		}
		if(maxabs <= eps)
			break;
		logX -= powX/(double)k;
		powX = powX * X;
	}
	logX *= powC;
	return logX;
}


///////////////////////////////////////////////////////////////////////////////
// Sqrtm(A)
// Computes matrix square root of A, X which is defined as
//   A = X * X
// where no transpose is involved.
// For symmetric postive definite matrices
//   [E D] = eig(A); sqrtm(A) = E * sqrt(D) * E';
// For general matrices
//   [E D] = eig(A); sqrtm(A) = E * sqrt(D) / E = mtl_mrdivide(E*sqrt(D), E);
// Note: Sqrtm is very fragile!
///////////////////////////////////////////////////////////////////////////////
template <typename T> Matrix<T> Sqrtm(const Matrix<T>& A) {
	int dim=NumRows(A);
	assert(NumRows(A)==NumCols(A));
	if(A.IsZero()) return Matrix<T>(dim,dim);
	Matrix<T> ans;
	if(A.IsSymmetric() && A.IsPD()){
		Matrix<T> E,D;
		int retCode=EigS(A,E,D);	
		for(int i=1;i<=dim;i++) {
			if(real(D[i][i])<0) D[i][i]=0;
			else D[i][i]=T( sqrt(real(D[i][i])) );
		}
		ans=E*D*Transpose(E);
	}
	else if(A.IsReal()){
		ans=SqrtmGenReal(A);
	}
	else{
		ans=SqrtmComplex(A);
	}
	double tol=1e-6*Norm(A,mtl_numeric_limits<int>::max());
	if((A-ans*ans).IsZero(tol) == false){
		cerr<<"Warning (Sqrtm): The answer may be inaccurate.\n";
		//assert(0);
	}
	return ans;
}


template <typename T> Matrix<T> SqrtmComplex(const Matrix<T>& A) {
	int i;
	int dim=NumRows(A);
	Matrix<T> V,D;
	int isOK=Eig(A,V,D);
	for(i=1;i<=dim;i++){
		D[i][i]=sqrt(D[i][i]);
	}
	return mtl_mrdivide(V*D,V);
}


inline Matrix<float> SqrtmComplex(const Matrix<float>& A) {
	cerr << "ERROR: Must not arrive here.\n";
	assert(0);
	return Matrix<float>();
}


inline Matrix<double> SqrtmComplex(const Matrix<double>& A) {
	cerr << "ERROR: Must not arrive here.\n";
	assert(0);
	return Matrix<double>();
}


#ifndef DISABLE_COMPLEX
inline Matrix<complex<float> > SqrtmComplex(const Matrix<complex<float> >& A) {
	Matrix<complex<double> > AA=A.ToCDouble();
	AA = SqrtmComplex(AA);
	Matrix<complex<float> > C;
	C.FromCDouble(AA);
	return C;
}
#endif


template <typename T> Matrix<T> SqrtmGenReal(const Matrix<T>& A) {
	int dim=NumRows(A);
	int i;
	Matrix<T> V,D;
	int isOK=EigReal(A,V,D);
	if(D.IsDiagonal()==false){
		cerr<<"ERROR (SqrtmGenReal): Input matrix has complex eigenvalues\n";
		cerr<<"      Use a complex matrix.\n";
		assert(0);
		return Matrix<T>();
	}
	for(i=1;i<=dim;i++){
		if(D[i][i] < 0) {
			cerr<<"ERROR (SqrtmGenReal): Input matrix has negative eigenvalues\n";
			cerr<<"      Use a complex matrix.\n";
			assert(0);
			return Matrix<T>();
		}
	}
	// Now A has real positive eigenvalues and real eigenvectors.
	for(i=1;i<=dim;i++){
		D[i][i]=sqrt(D[i][i]);
	}
	return mtl_mrdivide(V*D,V);
}


#ifndef DISABLE_COMPLEX
inline Matrix<complex<float> > SqrtmGenReal(const Matrix<complex<float> >& A) {
	cerr << "ERROR: Must not arrive here.\n";
	assert(0);
	return 0;
}


inline Matrix<complex<double> > SqrtmGenReal(const Matrix<complex<double> >& A) {
	cerr << "ERROR: Must not arrive here.\n";
	assert(0);
	return 0;
}
#endif


inline Matrix<float> SqrtmGenReal(const Matrix<float>& A) {
	Matrix<double> AA=A.ToDouble();
	AA=SqrtmGenReal(AA);
	Matrix<float> C;
	C.FromDouble(AA);
	return C;
}


/**************************************************************************
**
** Copyright (C) 1993 David E. Stewart & Zbigniew Leyk, all rights reserved.
**
**			     Meschach Library
** 
** This Meschach Library is provided "as is" without any express 
** or implied warranty of any kind with respect to this software. 
** In particular the authors shall not be liable for any direct, 
** indirect, special, incidental or consequential damages arising 
** in any way from use of the software.
** 
** Everyone is granted permission to copy, modify and redistribute this
** Meschach Library, provided:
**  1.  All copies contain this copyright notice.
**  2.  All modified copies shall carry a notice stating who
**      made the last modification and the date of such modification.
**  3.  No charge is made for this software or works derived from it.  
**      This clause shall not be construed as constraining other software
**      distributed on the same medium as this software, nor is a
**      distribution fee considered a charge.
**
***************************************************************************/


/* _m_exp -- compute matrix exponential of A and save it in out
   -- uses Pade approximation followed by repeated squaring
   -- eps is the tolerance used for the Pade approximation 
   -- A is not changed
   -- q_out - degree of the Pade approximation (q_out,q_out)
   -- j_out - the power of 2 for scaling the matrix A
              such that ||A/2^j_out|| <= 0.5
*/
template <typename T> int Matrix<T>::mes_exp(Matrix<T>& A,double eps,Matrix<T>& out,int* q_out,int* j_out)
{
	static Matrix<T> D, Apow, N, Y;
	static Vector<T> tmp;
	static Vector<double> c1;
	static Vector<int> pivot;
	int j, k, l, q, r, s, j2max, t, n;
	double inf_norm, eqq, power2, c, sign;
	
	out.Resize(A.m,A.n);
	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SIZES,"_m_exp");
	if ( &A[1] == &out[1] )
		mes_error(E_INSITU,"_m_exp");
	if ( eps < 0.0 )
		mes_error(E_RANGE,"_m_exp");
	else if (eps == 0.0)
		eps = MACH_EPS;
	
	N.Resize(NumRows(A),NumCols(A));
	D.Resize(NumRows(A),NumCols(A));
	Apow.Resize(NumRows(A),NumCols(A));
	out.Resize(NumRows(A),NumCols(A));
	
	/* normalise A to have ||A||_inf <= 1 */
	inf_norm = Norm(A,mtl_numeric_limits<int>::max());
	if (inf_norm <= 0.0) {
		out.SetIdentity();
		*q_out = -1;
		*j_out = 0;
		return 1;
	}
	else {
		j2max = (int)(floor(1+log(inf_norm)/log(2.0))+MACH_EPS);
		j2max = local_max(0, j2max);
		//j2max = 2*j2max; // increase precision
	}
	
	power2 = 1.0;
	for ( k = 1; k <= j2max; k++ )
		power2 *= 2;
	power2 = 1.0/power2;
	if ( j2max > 0 )
		A*=power2;
	
	/* compute order for polynomial approximation */
	eqq = 1.0/6.0;
	for ( q = 1; eqq > eps; q++ )
		eqq /= 16.0*(2.0*q+1.0)*(2.0*q+3.0);
	
	/* construct vector of coefficients */
	c1.Resize(q+1);
	c1[1] = 1.0;
	for ( k = 1; k <= q; k++ ) 
		c1[k+1] = c1[k]*(q-k+1)/((2*q-k+1)*(double)k);
	
	tmp.Resize(NumCols(A));
	
	s = (int)floor(sqrt((double)q/2.0));
	if ( s <= 0 )  s = 1;
	Apow=Pow(A,s);
	r = q/s;
	
	Y.Resize(s,NumCols(A));
	/* y0 and y1 are pointers to rows of Y, N and D */
	Y=T();
	N=T();
	D=T();

	for( j = 1; j <= NumCols(A); j++ ){
		if (j > 1)
			Y[1][j-1] = 0.0;
		Y[1][j]=1.0;
		for ( k = 1; k <= s-1; k++ ){
			tmp=Y.Row(k);
			tmp.SetType(COL_VECTOR);
			tmp=(A)*(tmp);
			tmp.SetType(ROW_VECTOR);
			Y.SetRow(k+1,tmp);
		}

		t = s*r;
		for ( l = 1; l <= q-t+1; l++ ){
			assert(t+l <= q+1);
			c = c1[t+l];
			sign = ((t+l-1) & 1) ? -1.0 : 1.0;
			for(n=1;n<=NumCols(Y);n++)
				N[j][n] += T(c)*Y[l][n];
			for(n=1;n<=NumCols(Y);n++)
				D[j][n] += T(c*sign)*Y[l][n];
		}

		for (k=1; k <= r; k++){
			tmp=N.Row(j);
			tmp.SetType(COL_VECTOR);
			tmp=Apow*tmp;
			tmp.SetType(ROW_VECTOR);
			N.SetRow(j,tmp);
			tmp=D.Row(j);
			tmp.SetType(COL_VECTOR);
			tmp=Apow*tmp;
			tmp.SetType(ROW_VECTOR);
			D.SetRow(j,tmp);
			t = s*(r-k);
			for (l=1; l <= s; l++){
				c = c1[t+l];
				sign = ((t+l-1) & 1) ? -1.0 : 1.0;
				for(n=1;n<=NumCols(Y);n++)
					N[j][n] += T(c)*Y[l][n];
				for(n=1;n<=NumCols(Y);n++)
					D[j][n] += T(c*sign)*Y[l][n];
			}
		}
	}
	
	pivot.Resize(NumRows(A));
	
	/* note that N and D are transposed, therefore we use LUTsolve;
	out is saved row-wise, and must be transposed after this */
	
	mes_LUfactor(D,pivot);
	for (k=1; k <= NumCols(A); k++){
		mes_LUTsolve(D,pivot,N[k],out[k]);
	}
	out=Transpose(out);	

	/* Use recursive squaring to turn the normalised exponential to the
	true exponential */
	
	for( k = 1; k <= j2max; k++){
		if(k&1) Apow=out*out;
		else    out=Apow*Apow;
	}
	
	if( ( ((k)&1)? (&Apow[1]) : (&out[1]) ) == &out[1]){
		out=Apow;
	}

	/* output parameters */
	*j_out = j2max;
	*q_out = q;
	
	/* restore the matrix A */
	A *= (1.0/power2);
	return 1;
}


#endif
