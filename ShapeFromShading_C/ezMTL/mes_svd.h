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
// Filename: mes_svd.h
// Revision:
//    1. Debugged an error in fixsvd().
//    2. Added maxIter in bisvd().
//    3. Added mes_rotsvd() and revised to support complex matrices.
//       Modified the contribution from Gary. D. Brushe.
// Note:
//    1. I strongly recommend the routines in this file be called 
//       with double precision arguments.
///////////////////////////////////////////////////////////////////////////////

#ifndef _MES_SVD_H_
#define _MES_SVD_H_

/* ------------- Singular Value Decomposition --------------------------- */
/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
**
**                          Meschach Library
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


///////////////////////////////////////////////////////////////////////////////
// External routines
//   Just for your information.
///////////////////////////////////////////////////////////////////////////////

// from Givens
//template <typename T> void	mes_givens(T x,T y,double* c,T* s);
//template <typename T> void	mes_rot_rows(Matrix<T>& mat,int i,int k,double c,T s,Matrix<T>& out);

// from Householder
//template <typename T> void	mes_hhvec(Vector<T>& vec,int i0,double* beta,Vector<T>& out,T* newval);
//template <typename T> void	mes_hhtrvec(Vector<T>& hh,double beta,int i0,Vector<T>& in,Vector<T>& out);
//template <typename T> void	mes_hhtrrows(Matrix<T>& M,int i0,int j0,Vector<T>& hh,double beta);
//template <typename T> void	mes_hhtrcols(Matrix<T>& M,int i0,int j0,Vector<T>& hh,double beta);

///////////////////////////////////////////////////////////////////////////////
// fixsvd -- fix minor details about SVD
//	-- make singular values non-negative
//	-- sort singular values in decreasing order
//	-- variables as for bisvd()
//	-- no argument checking
// Note:
//   Removed an error in quicksort.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_fixsvd(Vector<T>& d, Matrix<T>& U, Matrix<T>& V)
{
	const int MAX_STACK=100;
    int		i, j, k, l, r, stack[MAX_STACK], sp;
    T	tmp, dtmp;
	double v;
	int dim = d.Size();
	
    /* make singular values real non-negative */
    for ( i = 1; i <= dim; i++ ) {
		if ( real(d[i]) < 0.0 ){
			d[i] = - d[i];
			if ( NumCols(U) ){
				for ( j = 1; j <= NumRows(U); j++ )
					U[i][j] = - U[i][j];
			}
		}
	}
		
	/* sort singular values */
	/* nonrecursive implementation of quicksort due to R.Sedgewick,
	"Algorithms in C", p. 122 (1990) */
	sp = -1;
	l = 1;	r = dim;
	while(true){
		while ( r > l ){
			/* i = partition(d,l,r) */
			v = real(d[r]);
			
			i = l - 1;	    j = r;
			for ( ; ; ){	/* inequalities are "backwards" for **decreasing** order */
#if 0
				/* the old code with errors */
				while ( d[++i] > v )
					;
				while ( d[--j] < v )
					;
#endif
				/* the new code */
				while(1){
					++i;
					if(i>dim || real(d[i])<=v) break;
				}
				while(1){
					--j;
					if(j<1 || real(d[j])>=v) break;
				}

				if ( i >= j )
					break;
				/* swap entries in d */
				dtmp = d[i];
				d[i] = d[j];
				d[j] = dtmp;
				/* swap rows of U & V as well */
				if ( NumCols(U) ){
					for ( k = 1; k <= NumCols(U); k++ ) {
						tmp = U[i][k];
						U[i][k] = U[j][k];
						U[j][k] = tmp;
					}
				}
				if ( NumCols(V) ){
					for ( k = 1; k <= NumCols(V); k++ ){
						tmp = V[i][k];
						V[i][k] = V[j][k];
						V[j][k] = tmp;
					}
				}
			}
			dtmp = d[i];    d[i] = d[r];    d[r] = dtmp;
			if ( NumCols(U) ){
				for ( k = 1; k <= NumCols(U); k++ ){
					tmp = U[i][k];
					U[i][k] = U[r][k];
					U[r][k] = tmp;
				}
			}
			if ( NumCols(V) ){
				for ( k = 1; k <= NumCols(V); k++ ){
					tmp = V[i][k];
					V[i][k] = V[r][k];
					V[r][k] = tmp;
				}
			}
			/* end i = partition(...) */
			if(sp>=MAX_STACK-3) {
				cerr << "ERROR: Stack overflow in fixsvd()\n";
				assert(0);
				return 0; 
			} 
			if ( i - l > r - i ){
				stack[++sp] = l;
				stack[++sp] = i-1;
				l = i+1;
			}
			else{
				stack[++sp] = i+1;  
				stack[++sp] = r;	
				r = i-1;    
			}
		}
		if ( sp < 0 )
			break;
		r = stack[sp--];
		l = stack[sp--];
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// bisvd -- svd of a bidiagonal m x n matrix represented by d (diagonal) and
//			f (super-diagonals)
//	-- returns with d set to the singular values, f zeroed
//	-- if the orthogonal operations are accumulated
//		in U, V; if U, V == I on entry, then SVD == U^T.A.V
//		where A is initial matrix
//	-- returns d on exit */
// Note:
//    The accuracy of sqrt(x^2+y^2) is critical in stability of the algorithm.
//    Must use 'pithagoras()' to be safe.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_bisvd(Vector<T>& dd, Vector<T>& ff, Matrix<T>& U, Matrix<T>& V, int maxIter)
{
	int	i, j, n;
	int	i_min, i_max, split, iter, isOK=1;
	T	s, z, dtmp, t11, t12, t21, t22;
	double c, size;
	bool isReal=U.IsReal();
	
	if (dd.Size()==0 || ff.Size()==0 )
		mes_error(E_NULL,"bisvd");
	if ( dd.Size() != ff.Size() + 1 )
		mes_error(E_SIZES,"bisvd");
	
    n = dd.Size();
	if ( (NumCols(U) < n ) || ( NumRows(V) < n ) )
		mes_error(E_SIZES,"bisvd");
	if ( (NumRows(U) != NumCols(U)) || 
		(NumRows(V) != NumCols(V)) )
		mes_error(E_SQUARE,"bisvd");
	   
	if ( n == 1 )
		return isOK;
	   
	size = mes_norm_inf(dd) + mes_norm_inf(ff);
	   
	i_min = 1;
	while ( i_min <= n ){	/* outer while loop */
	   /* find i_max to suit;
		submatrix i_min..i_max should be irreducible */
		i_max = n;
		for ( i = i_min; i <= n-1; i++ ){
			if ( mtl_iszero(dd[i]) || mtl_iszero(ff[i]) ) {
				i_max = i;
				if ( mtl_iszero(ff[i]) == false ) {
					/* have to ``chase'' ff[i] element out of matrix */
					z = ff[i];
					ff[i] = 0;
					for ( j = i; j <= n-1 && mtl_iszero(z)==false; j++ ) {
						mes_givens(dd[j+1],z, &c, &s);
						dd[j+1] =  c*dd[j+1] - s*z;
						if ( j+1 <= n-1 ) {
							z       = conj(s)*ff[j+1];
							ff[j+1] = c*ff[j+1];
						}
						if(NumCols(U))
							mes_rot_rows(U,i,j+1,c,s,U);
					}
				}
				break;
			}
		}
		
		
		if ( i_max <= i_min ) {
			i_min = i_max + 1;
			continue;
		}
		
		split = 0;
		iter=0;
		while ( ! split ){
			double shift;
			iter++;

			/* compute shift */
			t11 = dd[i_max-1]*conj(dd[i_max-1]) + (i_max > i_min+1 ? ff[i_max-2]*conj(ff[i_max-2]) : 0.0);
			t12 = conj(dd[i_max-1])*ff[i_max-1];
			t21 = conj(t12);
			t22 = dd[i_max]*conj(dd[i_max]) + ff[i_max-1]*conj(ff[i_max-1]);

			if(isReal){
				/* use e-val of [[t11,t12],[t12,t22]] matrix closest to t22 */
				double diff = real(t11-t22)/2;
				shift = real(t22 - t12*t12/(diff + ((diff>=0) ? 1 : -1)*pythagoras(diff,t12)));
			}
			else{	
				/* use the smallest eigenvalue of the matrix */
				double ztmp=0.5*(real(t11)-real(t22));
				double diff=sqrt(ztmp*ztmp+real(t12*t21));
				double sum=0.5*(real(t11)+real(t22));
				double det = real(t11)*real(t22)-real(t12*t21);
				if(Abs(sum+diff) > Abs(sum-diff))
					shift=det/(sum+diff);
				else
					shift=det/(sum-diff);
			}

			/* perturb shift if convergence is slow */
			if(iter%10==0){
				shift += iter*0.02;
			}
			
			/* initial Givens' rotation */
			mes_givens(conj(dd[i_min])*dd[i_min]-shift, conj(dd[i_min])*ff[i_min], &c, &s);
			
			/* do initial Givens' rotations */
			dtmp        = c*dd[i_min] - s*ff[i_min];
			ff[i_min]   = c*ff[i_min] + conj(s)*dd[i_min];
			dd[i_min]   = dtmp;
			z           = -s*dd[i_min+1];
			dd[i_min+1] =  c*dd[i_min+1];
			if ( NumCols(V) )
				mes_rot_rows(V,i_min,i_min+1,c,conj(s),V);
			
			/* 2nd Givens' rotation */
			mes_givens(dd[i_min],z, &c, &s);
			dd[i_min]   = c*dd[i_min] - s*z;
			dtmp        = c*dd[i_min+1] + conj(s)*ff[i_min];
			ff[i_min]   = c*ff[i_min] - s*dd[i_min+1];
			dd[i_min+1] = dtmp;
			if ( i_min+1 < i_max ){
				z           = -s*ff[i_min+1];
				ff[i_min+1] =  c*ff[i_min+1];
			}
			if ( NumCols(U) )
				mes_rot_rows(U,i_min,i_min+1,c,s,U);
			
			for ( i = i_min+1; i < i_max; i++ ){
				/* get Givens' rotation for zeroing z */
				mes_givens(ff[i-1],z, &c, &s);
				ff[i-1] = c*ff[i-1] - s*z;
				dtmp    = c*dd[i] - s*ff[i];
				ff[i]   = c*ff[i] + conj(s)*dd[i];
				dd[i]   = dtmp;
				z       = -s*dd[i+1];
				dd[i+1] =  c*dd[i+1];
				if ( NumCols(V) )
					mes_rot_rows(V,i,i+1,c,conj(s),V);
				
				/* get 2nd Givens' rotation */
				mes_givens(dd[i],z, &c, &s);
				dd[i]   = c*dd[i] - s*z;
				dtmp    = c*dd[i+1] + conj(s)*ff[i];
				ff[i]   = c*ff[i] - s*dd[i+1];
				dd[i+1] = dtmp;
				if ( i+1 < i_max ){
					z       = -s*ff[i+1];
					ff[i+1] =  c*ff[i+1];
				}
				if ( NumCols(U) )
					mes_rot_rows(U,i,i+1,c,s,U);
			}
			
			/* should matrix be split? */
			for ( i = i_min; i < i_max; i++ ){
				if ( Abs(ff[i]) < (MACH_EPS*(Abs(dd[i])+Abs(dd[i+1]))) ){
					split = 1;
					ff[i] = 0.0;
				}
				else if ( Abs(dd[i]) < (MACH_EPS*size) ){
					split = 1;
					dd[i] = 0.0;
				}
			}
			
			if(iter==maxIter){
				cerr << "ERROR: No convergence in " << maxIter << " mes_bisvd iterations" << endl;
				isOK=0;
				assert(0);
				return isOK;
			}
		}
	}

	/* rotate dd[i] so it is real and rotate U & V appropriately */
	isOK=mes_rotsvd(dd,U,V);

	/* sort singular values and singular vectors in the descending order of singular values */
	isOK=mes_fixsvd(dd,U,V);
	
	return isOK;
}


template <typename T> int Matrix<T>::mes_rotsvd(Vector<T>& dd, Matrix<T>& U, Matrix<T>& V)
{
	int i,j;
	double size = mes_norm_inf(dd);
	for(i=1; i<=dd.Size();i++){
		if((Abs(dd[i])>MACH_EPS*size) && (Abs(imag(dd[i]))>MACH_EPS*size)){
			double theta=atan2( Abs(imag(dd[i])), Abs(real(dd[i])) );
			T dtmp1=T(mtl_sign(real(dd[i]))*cos(theta/2), -mtl_sign(imag(dd[i]))*sin(theta/2));
			T dtmp2=-dtmp1;
			dd[i]=dd[i]*T(dtmp1*dtmp2);
			for(j=1;j<=NumCols(U);j++){
				U[i][j]=U[i][j]*dtmp1;
			}
			for(j=1;j<=NumCols(V);j++){
				V[i][j]=V[i][j]*conj(dtmp2);
			}
		}
		if(Abs(imag(dd[i]))<MACH_EPS*size){
			dd[i]=T(real(dd[i]),0);
		}
		else{
			cerr << "ERROR: singular values must be real.\n";
			assert(0);
		}
	}
	return 1;
}


inline int Matrix<float>::mes_rotsvd(Vector<float>& dd, Matrix<float>& U, Matrix<float>& V)
{
	return 1;
}


inline int Matrix<double>::mes_rotsvd(Vector<double>& dd, Matrix<double>& U, Matrix<double>& V)
{
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// bifactor -- perform preliminary factorisation for bisvd
//	-- updates U and/or V, which ever is defined
// Note:
//   Do not change beta to single precision (float type).
//   The stability of the Householder transformation is very sensitive
//   to the precision of beta in hhvec() and hhtrrows().
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_bifactor(Matrix<T>& A, Matrix<T>& U, Matrix<T>& V)
{
	int compU=1, compV=1;
	int	i, k, m, n;
	static Vector<T> tmp1, tmp2;
	double	beta;

	m=NumRows(A);
	n=NumCols(A);
	if ( ( NumRows(U) != NumCols(U) ) || ( NumRows(V) != NumCols(V) ) )
		mes_error(E_SQUARE,"bifactor");
	if ( ( NumRows(U) != m ) || ( NumRows(V) != n ) )
		mes_error(E_SIZES,"bifactor");
	tmp1.Resize(m);
	tmp2.Resize(n);
	
	if(m>=n){
		for ( k = 1; k <= n; k++ ){
			mes_get_col(A,k,tmp1);
			mes_hhvec(tmp1,k,&beta,tmp1,&A[k][k]);
			mes_hhtrcols(A,k,k+1,tmp1,beta);
			if ( NumCols(U) )
				mes_hhtrcols(U,k,1,tmp1,beta);
			if ( k >= n-1 )
				continue;
			mes_get_row(A,k,tmp2);
			mes_hhvec(tmp2,k+1,&beta,tmp2,&A[k][k+1]);
			for(i=1;i<=tmp2.Size();i++) // conujugate householder vector for V
				tmp2[i]=conj(tmp2[i]);
			mes_hhtrrows(A,k+1,k+1,tmp2,beta);
			if ( NumCols(V) )
				mes_hhtrcols(V,k+1,1,tmp2,beta);
		}
	}
	else{
		for ( k = 1; k <= m; k++ ){
			mes_get_row(A,k,tmp1);
			for(i=1;i<=tmp1.Size();i++) // conujugate column vector of A
				tmp1[i]=conj(tmp1[i]);
			mes_hhvec(tmp1,k,&beta,tmp1,&A[k][k]);
			mes_hhtrrows(A,k+1,k,tmp1,beta);
			if ( NumCols(V) )
				mes_hhtrcols(V,k,1,tmp1,beta);
			if ( k >= m-1 )
				continue;
			mes_get_col(A,k,tmp2);
			mes_hhvec(tmp2,k+1,&beta,tmp2,&A[k+1][k]);
			mes_hhtrcols(A,k+1,k+1,tmp2,beta);
			if ( NumCols(U) )
				mes_hhtrcols(U,k+1,1,tmp2,beta);
		}
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// svd -- returns vector of singular values in d
//	-- also updates U and/or V, if one or the other is defined
//	-- destroys A
// Note:
//    Singular values for MxN (M<N) can be computed from SVD of A'
//    by conjugating A and interchanging the roles of U and V.
//      A = U' * S * V
//      A' = V' * S * U
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_svd(const Matrix<T>& A, Matrix<T>& U, Matrix<T>& V, Vector<T>& dd)
{
	static Vector<T>	ff;
	int	i, limit, tryX, isOK;
	
#if 0
	if(NumRows(A)<NumCols(A)){ // This implementation wastes memory.
		return mes_svd(Transpose(A),V,U,dd);
	}
#endif
	if(A.IsDouble()==false)
		mes_warn(WARN_SINGLE_PRECISION,"mes_svd");
	int m = NumRows(A);
	int n = NumCols(A);
	U.Resize(m,m);
	V.Resize(n,n);
	if ( (NumRows(U) != NumCols(U)) || 
		   (NumRows(V) != NumCols(V)) )
	   mes_error(E_SQUARE,"svd");
	if ( (NumRows(U) != NumRows(A) ) || ( NumRows(V) != NumCols(A) ) )
		mes_error(E_SIZES,"svd");

	int maxIter=100;
	Matrix<T> A_tmp;
	limit = local_min(m,n);
	dd.Resize(limit);
	ff.Resize(limit-1);
	for(tryX=1;tryX<=2;tryX++){
		A_tmp=A;
		U.SetIdentity();
		V.SetIdentity();		
		isOK=mes_bifactor(A_tmp,U,V);
		if(isOK==0)
			break;
		if(m>=n){
			for ( i = 1; i <= limit; i++ ){
				dd[i] = (A_tmp[i][i]);
				if ( i+1 <= limit )
					ff[i] = (A_tmp[i][i+1]);
			}
		}
		else{
			for ( i = 1; i <= limit; i++ ){
				dd[i] = (A_tmp[i][i]);
				if ( i+1 <= limit )
					ff[i] = conj(A_tmp[i+1][i]);
			}
		}
		
		if(tryX==2) maxIter=5*maxIter;
		if(m>=n)
			isOK=mes_bisvd(dd,ff,U,V,maxIter);
		else
			isOK=mes_bisvd(dd,ff,V,U,maxIter);

		if(isOK != 0) // Succeeded!
			break;
	}
	if(isOK == 0){
		cerr << "ERROR: No convergence in " << maxIter << " bisvd iterations" << endl;
		assert(0);
	}
	return isOK;
}


#endif


