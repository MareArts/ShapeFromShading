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
// Filename: mes_eig.h
// Revision:
//    1. Revised to support complex matrices.
// Note:
//    1. I strongly recommend the routines in this file be called 
//       with double precision arguments.
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_EIG_H_
#define	_MES_EIG_H_	

/* ----- Eigen Analysis for Symmetric Matrices and General Matrices ----- */
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


template <typename T> inline T Matrix<T>::mes_in_prod(const Vector<T>& a,const Vector<T>& b) {
	return mes_in_prod_offset(a,b,1);
}


/* in_prod_ -- inner product of two vectors from i0 downwards */
template <typename T> T Matrix<T>::mes_in_prod_offset(const Vector<T>& a,const Vector<T>& b,int i0){
	int	i,limit;
	limit = local_min(a.Size(),b.Size());
	assert ( i0 <= limit );
	T sum=T();
	for(i=i0;i<=limit;i++) sum += conj(a[i])*b[i];
	return sum;
}


/* mes_norm2_offset -- norm of vectors from i0 downwards */
template <typename T> double Matrix<T>::mes_norm2_offset(const Vector<T>& a,int i0){
	int	i,limit;
	limit = a.Size();
	assert ( i0 >= 1 && i0 <= limit );
	double sum=0, z=0;
	for(i=i0;i<=limit;i++) if(z<Abs(a[i])) z=Abs(a[i]);
	if(z==0) return 0;
	for(i=i0;i<=limit;i++) {
		double x=Abs(a[i]/T(z));
		sum += x*x;
	}
	return (z*sqrt(sum));
}


/* _v_norm_inf -- computes (scaled) infinity-norm (supremum norm) of vectors */
template <typename T> double	Matrix<T>::mes_norm_inf(const Vector<T>& x){
	return Norm(x,mtl_numeric_limits<int>::max());
}


/* v_norm2 -- computes (scaled) 2-norm (Euclidean norm) of vectors */
template <typename T> double	Matrix<T>::mes_norm2(const Vector<T>& x){
	return Norm(x,2);
}


/* get_row -- gets a specified row of a matrix and retruns it as a vector */
template <typename T> void	Matrix<T>::mes_get_row(const Matrix<T>& mat,int row,Vector<T>& vec){
   int	i;
   assert ( row <= NumRows(mat) );
   if ( vec.Size() != NumCols(mat) )
     vec.Resize(NumCols(mat));
   
   for ( i=1; i<=NumCols(mat); i++ )
     vec[i] = mat[row][i];
   
   return;
}


/* get_col -- gets a specified column of a matrix and retruns it as a vector */
template <typename T> void	Matrix<T>::mes_get_col(const Matrix<T>& mat,int col,Vector<T>& vec){
   int	i;
   assert ( col <= NumCols(mat) );
   if ( vec.Size() != NumRows(mat) )
     vec.Resize(NumRows(mat));
   
   for ( i=1; i<=NumRows(mat); i++ )
     vec[i] = mat[i][col];
   
   return;
}


/* set_row -- sets row of matrix to values given in vec (in situ) */
template <typename T> void	Matrix<T>::mes_set_row(Matrix<T>& mat,int row,Vector<T>& vec){
   int	i,lim;
   
   assert( row <= NumRows(mat) );
   lim = local_min(NumCols(mat),vec.Size());
   for ( i=1; i<=lim; i++ )
     mat[row][i] = vec[i];

}


/* set_col -- sets column of matrix to values given in vec (in situ) */
template <typename T> void	Matrix<T>::mes_set_col(Matrix<T>& mat,int col,Vector<T>& vec){
   int	i,lim;
   
   assert( col <= NumCols(mat) );
   lim = local_min(NumRows(mat),vec.Size());
   for ( i=1; i<=lim; i++ )
     mat[i][col] = vec[i];

}


/* v_copy_ -- copies vector into new area */
template <typename T> void	Matrix<T>::mes_copy_offset(const Vector<T>& in,Vector<T>& out,int i0){
	if ( &in[1]==&out[1] )
		return;
	if ( out.Size() != in.Size() )
		out.Resize(in.Size());
	for (int i=i0; i <= in.Size(); i++ )
		out[i] = in[i];
}


///////////////////////////////////////////////////////////////////////////////
// trieig -- finds eigenvalues of symmetric tridiagonal matrices
//	-- matrix represented by a pair of vectors a (diag entries)
//		and b (sub- & super-diag entries)
//	-- eigenvalues in a on return	
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_trieig(Vector<T>& a,Vector<T>& b,Matrix<T>& Q,int maxIter){
	int	i, i_min, i_max, n, split, iter, isOK=1;
	double	bk, ak1, bk1, ak2, bk2, z;
	double	c, c2, cs, s, s2, d, mu;
	static double  local_SQRT2=1.4142135623730949;

	if ( a.Size() != b.Size() + 1 || ( NumRows(Q) != a.Size() ) )
		mes_error(E_SIZES,"trieig");
	if ( NumRows(Q) != NumCols(Q) )
		mes_error(E_SQUARE,"trieig");

	n = a.Size();

	i_min = 1;
	while ( i_min <= n ){		/* outer while loop */
		/* find i_max to suit;
			submatrix i_min..i_max should be irreducible */
		i_max = n;
		for ( i = i_min; i < n; i++ )
			if ( b[i] == 0.0 ){
				i_max = i;	break;
			}
		if ( i_max <= i_min ){
		    i_min = i_max + 1;
		    continue;	/* outer while loop */
		}

		/* repeatedly perform QR method until matrix splits */
		split = false;
		iter=0;
		while ( ! split ){		/* inner while loop */
			int sign_d;
			iter++;

			/* find Wilkinson shift */
			d = (a[i_max-1] - a[i_max])/2;
			sign_d = (d >= 0 ? 1 : -1 );
			mu = a[i_max] - (b[i_max-1]*b[i_max-1])/(d + sign_d*pythagoras(d,b[i_max-1]));
			
			/* initial Givens' rotation */
			mes_givens(a[i_min]-mu,b[i_min],&c,&s);

			if ( Abs(c) < local_SQRT2 ){
				c2 = c*c;	s2 = 1-c2;
			}
			else{	
				s2 = s*s;	c2 = 1-s2;
			}
			cs = c*s;
			ak1 = c2*a[i_min]+s2*a[i_min+1]-2*cs*b[i_min];
			bk1 = cs*(a[i_min]-a[i_min+1]) + (c2-s2)*b[i_min];
			ak2 = s2*a[i_min]+c2*a[i_min+1]+2*cs*b[i_min];
			bk2 = ( i_min < i_max-1 ) ? c*b[i_min+1] : 0.0;
			z  = ( i_min < i_max-1 ) ? -s*b[i_min+1] : 0.0;
			a[i_min] = ak1;
			a[i_min+1] = ak2;
			b[i_min] = bk1;
			if ( i_min < i_max-1 )
				b[i_min+1] = bk2;
			if ( NumCols(Q) )
				mes_rot_cols(Q,i_min,i_min+1,c,s,Q);
			
			for ( i = i_min+1; i < i_max; i++ ){
				/* get Givens' rotation for sub-block -- k == i-1 */
				mes_givens(b[i-1],z,&c,&s);
				
				/* perform Givens' rotation on sub-block */
				if ( Abs(c) < local_SQRT2 ){
					c2 = c*c;	s2 = 1-c2;
				}
				else{
					s2 = s*s;	c2 = 1-s2;
				}
				cs = c*s;
				bk  = c*b[i-1] - s*z;
				ak1 = c2*a[i]+s2*a[i+1]-2*cs*b[i];
				bk1 = cs*(a[i]-a[i+1]) +
					(c2-s2)*b[i];
				ak2 = s2*a[i]+c2*a[i+1]+2*cs*b[i];
				bk2 = ( i+1 < i_max ) ? c*b[i+1] : 0.0;
				z  = ( i+1 < i_max ) ? -s*b[i+1] : 0.0;
				a[i] = ak1;	a[i+1] = ak2;
				b[i] = bk1;
				if ( i < i_max-1 )
					b[i+1] = bk2;
				if ( i > i_min )
					b[i-1] = bk;
				if ( NumCols(Q) )
					mes_rot_cols(Q,i,i+1,c,s,Q);
			}
			
			/* test to see if matrix should be split */
			for ( i = i_min; i < i_max; i++ ){
				if ( (Abs(b[i])) < (MACH_EPS*(Abs(a[i])+Abs(a[i+1]))) ){
					b[i] = 0.0;
					split = true;
				}
			}
	
			if(iter==maxIter){
				cerr << "ERROR: Too many iterations in mes_trieig().\n";
				isOK=0;
				assert(0);
				return isOK;
			}
		}
	}
	
	return isOK;
}


///////////////////////////////////////////////////////////////////////////////
// symmeig -- computes eigenvalues of a dense symmetric matrix
//	-- A **must** be symmetric on entry
//	-- eigenvalues stored in out
//	-- Q contains orthogonal matrix of eigenvectors
//	-- returns vector of eigenvalues
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_symmeig(const Matrix<T>& A,Matrix<T>& QQ,Vector<T>& oo){
	int	i, isOK;
	static Vector<double> bb, diag, beta;

	if(A.IsComplex())
		mes_error(E_NEED_REAL,"mes_symmeig");
	if(A.IsDouble()==false)
		mes_warn(WARN_SINGLE_PRECISION,"mes_symmeig");

	Matrix<double> AA=A;
	if ( NumRows(AA) != NumCols(AA) )
		mes_error(E_SQUARE,"symmeig");
	if ( oo.Size() != NumRows(AA) )
		oo.Resize(NumRows(AA));

	bb  .Resize(NumRows(AA) - 1);
	diag.Resize((int)NumRows(AA));
	beta.Resize((int)NumRows(AA));

	isOK=Matrix<double>::mes_Hfactor(AA,diag,beta);
	if(isOK==0) return isOK;
	if ( NumCols(QQ) ) {
		Matrix<double>::mes_makeHQ(AA,diag,beta,QQ);
		if(isOK==0) return isOK;
	}

	for ( i = 1; i <= NumRows(AA) - 1; i++ ){
		oo[i] = AA[i][i];
		bb[i] = AA[i][i+1];
	}
	oo[i] = AA[i][i];
	isOK=Matrix<double>::mes_trieig(oo,bb,QQ);

	return isOK;
}


///////////////////////////////////////////////////////////////////////////////
// Schur decomposition of a real non-symmetric matrix
///////////////////////////////////////////////////////////////////////////////

template <typename T> inline void	Matrix<T>::mes_hhldr3(double x, double y, double z, double *nu1, double *beta, double *newval){
	double	alpha;

	if ( x >= 0.0 )
		alpha = pythagoras3(x,y,z);
	else
		alpha = -pythagoras3(x,y,z);
	*nu1 = x + alpha;
	*beta = 1.0/(alpha*(*nu1));
	*newval = alpha;
}


template <typename T> 	int	Matrix<T>::mes_hhldr3cols(Matrix<T>& A, int k, int j0, double beta, double nu1, double nu2, double nu3){
	double	ip, prod;
	int	j, n;

	if ( k <= 0 || k+3 > NumRows(A)+1 || j0 <= 0 )
		mes_error(E_BOUNDS,"hhldr3cols");
	n = NumCols(A);

	for ( j = j0; j <= n; j++ ){
	    ip = nu1*A(k,j)+nu2*A(k+1,j)+nu3*A(k+2,j);
	    prod = ip*beta;	    
		A[k][j] += (- prod*nu1);
	    A[k+1][j] += (- prod*nu2);
	    A[k+2][j] += (- prod*nu3);
	}
	return 1;
}


template <typename T> 	int	Matrix<T>::mes_hhldr3rows(Matrix<T>& A, int k, int i0, double beta, double nu1, double nu2, double nu3){
	double	ip, prod;
	int	i, m;

	if ( k <= 0 || k+3 > NumCols(A)+1 )
		mes_error(E_BOUNDS,"hhldr3rows");
	m = NumRows(A);
	i0 = local_min(i0,m);

	for ( i = 1; i <= i0; i++ ){
	    ip = nu1*A(i,k)+nu2*A(i,k+1)+nu3*A(i,k+2);
	    prod = ip*beta;
		A[i][k] += (- prod*nu1);
	    A[i][k+1] += (- prod*nu2);
	    A[i][k+2] += (- prod*nu3);
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// schur -- computes the Schur decomposition of the matrix A in situ
//	-- optionally, gives Q matrix such that Q^T.A.Q is upper triangular
//	-- returns upper triangular Schur matrix
// Note: The precision of the eigenvector and temporary vectors is
//       critical in stability of the algorithm.
//       Must be 'double'.
//       Must be more careful in handling mathematical expressions
//       like x/sqrt(x^2+y^2).
// Note: This implementation is very sensitive to the condition number of 
//       the input matrix.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_schur(Matrix<T>& AA,Matrix<T>& QQ,int maxIter){
	int compQ=1;
    int		i, j, k, k_min, k_max, k_tmp, n, split, iter, isOK=1;
    double	beta2, c, discrim, dummy, nu1, s, tmp, x, y, z;
    double	sqrt_macheps;
    static	Vector<double>	diag, beta;
 
	if(AA.IsComplex())
		mes_error(E_NEED_REAL,"mes_schur");
	if(AA.IsDouble()==false)
		mes_warn(WARN_SINGLE_PRECISION,"mes_schur");

    if ( NumRows(AA) != NumCols(AA) || ( NumRows(QQ) != NumCols(QQ) ) )
		mes_error(E_SQUARE,"schur");
    if ( NumRows(QQ) != NumRows(AA) )
		mes_error(E_SIZES,"schur");
    n = NumCols(AA);
    diag.Resize(NumCols(AA));
    beta.Resize(NumCols(AA));
    /* compute Hessenberg form */
    isOK=Matrix<double>::mes_Hfactor(AA,diag,beta);
	if(isOK==0) return isOK;

    /* save QQ if necessary */
    if ( compQ ){
		isOK=Matrix<double>::mes_makeHQ(AA,diag,beta,QQ);
		if(isOK==0) return isOK;
	}
    isOK=Matrix<double>::mes_makeH(AA,AA);
	if(isOK==0) return isOK;
	
    sqrt_macheps = sqrt(MACH_EPS);
	
    k_min = 1;

    while ( k_min <= n ){
		double	a00, a01, a10, a11;
		double	scale, t, numer, denom;
		
		/* find k_max to suit: submatrix k_min..k_max should be irreducible */
		k_max = n;
		for ( k = k_min; k < k_max; k++ ){
			if ( AA(k+1,k) == 0.0 ){
				k_max = k;	break;	
			}
		}
		
		if ( k_max <= k_min ){
			k_min = k_max + 1;
			continue;		/* outer loop */
		}
		
		/* check to see if we have a 2 x 2 block with complex eigenvalues */
		if ( k_max == k_min + 1 ){
			a00 = AA(k_min,k_min);
			a01 = AA(k_min,k_max);
			a10 = AA(k_max,k_min);
			a11 = AA(k_max,k_max);
			tmp = a00 - a11;
			discrim = tmp*tmp + 4*a01*a10;
			if ( discrim < 0.0 ){	// yes -- e-vals are complex
				// -- put 2 x 2 block in form [a b; c a];
				// then eigenvalues have real part a & imag part sqrt(|bc|)
				numer = - tmp;
				double atmp = a01+a10;
				denom = ( atmp >= 0.0 ) ? (atmp + pythagoras(atmp,tmp)) : (atmp - pythagoras(atmp,tmp));
				if ( denom != 0.0 ){   // t = s/c = numer/denom
					t = numer/denom;
					scale = c = 1.0/pythagoras(1,t);
					s = -c*t;
				}
				else{
					c = 1.0;
					s = 0.0;
				}
				mes_rot_cols(AA,k_min,k_max,c,s,AA);
				mes_rot_rows(AA,k_min,k_max,c,s,AA);
				if ( compQ )
					mes_rot_cols(QQ,k_min,k_max,c,s,QQ);
				k_min = k_max + 1;
				continue;
			}
			else {// discrim >= 0; i.e. block has two real eigenvalues
				  // no -- e-vals are not complex;
				// split 2 x 2 block and continue
				// s/c = numer/denom
				numer = ( tmp >= 0.0 ) ? - tmp - sqrt(discrim) : - tmp + sqrt(discrim);
				denom = 2*a01;
				if ( Abs(numer) < Abs(denom) ){   // t = s/c = numer/denom
					t = numer/denom;
					scale = c = 1.0/pythagoras(1,t);
					s = -c*t;
				}
				else if ( numer != 0.0 ){   // t = c/s = denom/numer
					t = denom/numer;
					scale = 1.0/pythagoras(1,t);
					c = Abs(t)*scale;
					s = ( t >= 0.0 ) ? -scale : scale;
				}
				else{ // numer == denom == 0
					c = 0.0;
					s = -1.0;
				}
				mes_rot_cols(AA,k_min,k_max,c,s,AA);
				mes_rot_rows(AA,k_min,k_max,c,s,AA);
				if ( compQ )
					mes_rot_cols(QQ,k_min,k_max,c,s,QQ);
				k_min = k_max + 1;	// go to next block
				continue;
			}
		}
		
		/* now have r x r block with r >= 2: apply Francis QR step until block splits */
		split = 0;
		iter = 0;
		while ( ! split ){
			iter++;
			
			/* set up Wilkinson/Francis complex shift */
			k_tmp = k_max - 1;
			
			a00 = AA(k_tmp,k_tmp);
			a01 = AA(k_tmp,k_max);
			a10 = AA(k_max,k_tmp);
			a11 = AA(k_max,k_max);
			
			/* treat degenerate cases differently
			-- if there are still no splits after five iterations
			and the bottom 2 x 2 looks degenerate, force it to split */
			if ( iter >= 5 &&
				Abs(a00-a11) < sqrt_macheps*(Abs(a00)+Abs(a11)) &&
				(Abs(a01) < sqrt_macheps*(Abs(a00)+Abs(a11)) ||
				Abs(a10) < sqrt_macheps*(Abs(a00)+Abs(a11))) ) {
				if ( Abs(a01) < sqrt_macheps*(Abs(a00)+Abs(a11)) )
					AA[k_tmp][k_max]=0.0;
				if ( Abs(a10) < sqrt_macheps*(Abs(a00)+Abs(a11)) ){
					AA[k_max][k_tmp]=0.0;
					split = 1;
					continue;
				}
			}
			
			s = a00 + a11;
			t = a00*a11 - a01*a10;
			
			/* break loop if a 2 x 2 complex block */
			if ( k_max == k_min + 1 && s*s < 4.0*t ){
				split = 1;
				continue;
			}
			
			/* perturb shift if convergence is slow */
			if ( (iter % 10) == 0 ) {
				s += iter*0.02;
				t += iter*0.02;
			}
			
			/* set up Householder transformations */
			k_tmp = k_min + 1;
			a00 = AA(k_min,k_min);
			a01 = AA(k_min,k_tmp);
			a10 = AA(k_tmp,k_min);
			a11 = AA(k_tmp,k_tmp);
			x = a00*a00 + a01*a10 - s*a00 + t;
			y = a10*(a00+a11-s);
			if ( k_min + 2 <= k_max )
				z = a10*AA[k_min+2][k_tmp];
			else
				z = 0.0;
			
			for ( k = k_min; k <= k_max-1; k++ ){
				if ( k < k_max - 1 ){
					mes_hhldr3(x,y,z,&nu1,&beta2,&dummy);
					mes_hhldr3cols(AA,k,local_max(k-1,1),beta2,nu1,y,z);
					mes_hhldr3rows(AA,k,local_min(n,k+3),beta2,nu1,y,z);
					if ( compQ )
						mes_hhldr3rows(QQ,k,n,beta2,nu1,y,z);
				}
				else{
					mes_givens(x,y,&c,&s);
					mes_rot_cols(AA,k,k+1,c,s,AA);
					mes_rot_rows(AA,k,k+1,c,s,AA);
					if ( compQ )
						mes_rot_cols(QQ,k,k+1,c,s,QQ);
				}
				x = AA(k+1,k);
				if ( k <= k_max - 2 )
					y = AA(k+2,k);
				else
					y = 0.0;
				if ( k <= k_max - 3 )
					z = AA(k+3,k);
				else
					z = 0.0;
			}

			for ( k = k_min; k <= k_max-2; k++ ){
				/* zero appropriate sub-diagonals */	
				AA[k+2][k]=0.0;
				if ( k < k_max-2 )
					AA[k+3][k]=0.0;
			}
			
			/* test to see if matrix should split */
			for ( k = k_min; k < k_max; k++ ){
				if ( (Abs(AA[k+1][k])) < (MACH_EPS*(Abs(AA[k][k])+Abs(AA[k+1][k+1]))) ){
					AA[k+1][k] = 0.0;
					split = 1;
				}
			}
			if(iter==maxIter){
				cerr << "ERROR: Too many iterations in mes_schur().\n";
				isOK=0;
				assert(0);
				return isOK;
			}
		}	
    }
	
    /* polish up AA by zeroing strictly lower triangular elements
	   and small sub-diagonal elements */
	for ( i = 1; i <= NumRows(AA); i++ ){
		for ( j = 1; j < i-1; j++ ){
			AA[i][j] = 0.0;
		}
	}
	for ( i = 1; i <= NumRows(AA) - 1; i++ ){
		if ( (Abs(AA[i+1][i])) < (MACH_EPS*(Abs(AA[i][i])+Abs(AA[i+1][i+1]))) ){
			AA[i+1][i] = 0.0;
		}
	}
	
	return isOK;
}


///////////////////////////////////////////////////////////////////////////////
// schur_vals -- compute real & imaginary parts of eigenvalues
//	-- assumes T contains a block upper triangular matrix
//		as produced by schur()
//	-- real parts stored in real_pt, imaginary parts in imag_pt
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_schur_evals(const Matrix<T>& TT,Vector<T>& real_pt,Vector<T>& imag_pt){
	int	i, n;
	double	discrim;
	double	diff, sum, tmp;

	if ( NumRows(TT) != NumCols(TT) )
		mes_error(E_SQUARE,"schur_evals");
	n = NumCols(TT);
	real_pt.Resize((int)n);
	imag_pt.Resize((int)n);

	i = 1;
	while ( i <= n ){
		if ( i <= n-1 && TT[i+1][i] != 0.0 ){   /* should be a complex eigenvalue */
			sum  = 0.5*(TT[i][i]+TT[i+1][i+1]);
			diff = 0.5*(TT[i][i]-TT[i+1][i+1]);
			discrim = diff*diff + TT[i][i+1]*TT[i+1][i];
			if ( discrim < 0.0 ){	/* yes -- complex e-vals */
				real_pt[i] = real_pt[i+1] = sum;
				imag_pt[i] = sqrt(-discrim);
				imag_pt[i+1] = - imag_pt[i];
			}
			else{	/* no -- actually both real */
				tmp = sqrt(discrim);
				real_pt[i]   = sum + tmp;
				real_pt[i+1] = sum - tmp;
				imag_pt[i]   = imag_pt[i+1] = 0.0;
			}
			i += 2;
		}
		else{   /* real eigenvalue */
			real_pt[i] = TT[i][i];
			imag_pt[i] = 0.0;
			i++;
		}
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// schur_vecs -- returns eigenvectors computed from the real Schur
//		decomposition of a matrix
//	-- T is the block upper triangular Schur matrix
//	-- Q is the orthognal matrix where A = Q.T.Q^T
//	-- if Q is null, the eigenvectors of T are returned
//	-- X_re is the real part of the matrix of eigenvectors,
//		and X_im is the imaginary part of the matrix.
//	-- X_re is returned
// Note: The precision of the eigenvector and temporary vectors is
//       critical in stability of the algorithm.
//       Must be 'double'.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_schur_vecs(const Matrix<T>& TT,const Matrix<T>& QQ,Matrix<T>& XX_re,Matrix<T>& XX_im){
	int	i, j, limit;
	double	t11_re, t11_im, t12, t21, t22_re, t22_im;
	double	l_re, l_im, det_re, det_im, invdet_re, invdet_im,
		val1_re, val1_im, val2_re, val2_im,
		tmp_val1_re, tmp_val1_im, tmp_val2_re, tmp_val2_im;
	double	sum, diff, discrim, magdet, norm, scale;
	static Vector<double>	tmp1_re, tmp1_im, tmp2_re, tmp2_im;

	if ( NumRows(TT) != NumCols(TT) || NumRows(XX_re) != NumCols(XX_re) || ( NumRows(QQ) != NumCols(QQ) ) || ( NumRows(XX_im) != NumCols(XX_im) ) )
	    mes_error(E_SQUARE,"schur_vecs");
	if ( NumRows(TT) != NumRows(XX_re) || ( NumRows(TT) != NumRows(QQ) ) || ( NumRows(TT) != NumRows(XX_im) ) )
	    mes_error(E_SIZES,"schur_vecs");

	tmp1_re.Resize(NumRows(TT));
	tmp1_im.Resize(NumRows(TT));
	tmp2_re.Resize(NumRows(TT));
	tmp2_im.Resize(NumRows(TT));

	i = 1;
	while ( i <= NumRows(TT) ){
		if ( i <= NumRows(TT)-1 && TT[i+1][i] != 0.0 ){	/* complex eigenvalue */
			sum  = 0.5*(TT[i][i]+TT[i+1][i+1]);
			diff = 0.5*(TT[i][i]-TT[i+1][i+1]);
			discrim = diff*diff + TT[i][i+1]*TT[i+1][i];
			l_re = l_im = 0.0;
			if ( discrim < 0.0 ){	/* yes -- complex e-vals */
				l_re = sum;
				l_im = sqrt(-discrim);
			}
			else /* not correct Real Schur form */
				mes_error(E_RANGE,"schur_vecs");
		}
		else{
			l_re = TT[i][i];
			l_im = 0.0;
		}

	    tmp1_im=0;
		Matrix<double>::mes_v_rand(tmp1_re);
	    tmp1_re=MACH_EPS*tmp1_re;
		
		/* solve (T-l.I)x = tmp1 */
		limit = ( l_im != 0.0 ) ? i+1 : i;
		for ( j = limit+1; j <= NumRows(TT); j++ )
			tmp1_re[j] = 0.0;
		j = limit;
		while ( j >= 1 ){
			double sum;
			int k;
			if ( j > 1 && TT[j][j-1] != 0.0 ){   /* 2 x 2 diagonal block */
				for(sum=0,k=j+1;k<=limit;k++) sum += tmp1_re[k]*TT[j-1][k]; 
				val1_re = tmp1_re[j-1] - sum;
				for(sum=0,k=j+1;k<=limit;k++) sum += tmp1_im[k]*TT[j-1][k]; 
				val1_im = tmp1_im[j-1] - sum;
				for(sum=0,k=j+1;k<=limit;k++) sum += tmp1_re[k]*TT[j][k]; 
				val2_re = tmp1_re[j] - sum;
				for(sum=0,k=j+1;k<=limit;k++) sum += tmp1_im[k]*TT[j][k]; 
				val2_im = tmp1_im[j] - sum;
				
				t11_re = TT[j-1][j-1] - l_re;
				t11_im = - l_im;
				t22_re = TT[j][j] - l_re;
				t22_im = - l_im;
				t12 = TT[j-1][j];
				t21 = TT[j][j-1];
				
				scale =  Abs(TT[j-1][j-1]) + Abs(TT[j][j]) + Abs(t12) + Abs(t21) + Abs(l_re) + Abs(l_im);
				
				det_re = t11_re*t22_re - t11_im*t22_im - t12*t21;
				det_im = t11_re*t22_im + t11_im*t22_re;
				magdet = det_re*det_re+det_im*det_im;
				if ( sqrt(magdet) < MACH_EPS*scale ){
					det_re = MACH_EPS*scale;
					magdet = det_re*det_re+det_im*det_im;
				}
				invdet_re =   det_re/magdet;
				invdet_im = - det_im/magdet;
				tmp_val1_re = t22_re*val1_re-t22_im*val1_im-t12*val2_re;
				tmp_val1_im = t22_im*val1_re+t22_re*val1_im-t12*val2_im;
				tmp_val2_re = t11_re*val2_re-t11_im*val2_im-t21*val1_re;
				tmp_val2_im = t11_im*val2_re+t11_re*val2_im-t21*val1_im;
				tmp1_re[j-1] = invdet_re*tmp_val1_re - invdet_im*tmp_val1_im;
				tmp1_im[j-1] = invdet_im*tmp_val1_re + invdet_re*tmp_val1_im;
				tmp1_re[j]   = invdet_re*tmp_val2_re - invdet_im*tmp_val2_im;
				tmp1_im[j]   = invdet_im*tmp_val2_re + invdet_re*tmp_val2_im;
				j -= 2;
			}
			else{
				t11_re = TT[j][j] - l_re;
				t11_im = - l_im;
				magdet = t11_re*t11_re + t11_im*t11_im;
				scale = Abs(TT[j][j]) + Abs(l_re);
				if ( sqrt(magdet) < MACH_EPS*scale ){
					t11_re = MACH_EPS*scale;
					magdet = t11_re*t11_re + t11_im*t11_im;
				}
				invdet_re =   t11_re/magdet;
				invdet_im = - t11_im/magdet;
				for(sum=0,k=j+1;k<=limit;k++) sum += tmp1_re[k]*TT[j][k]; 
				val1_re = tmp1_re[j] - sum;
				for(sum=0,k=j+1;k<=limit;k++) sum += tmp1_im[k]*TT[j][k]; 
				val1_im = tmp1_im[j] - sum;
				tmp1_re[j] = invdet_re*val1_re - invdet_im*val1_im;
				tmp1_im[j] = invdet_im*val1_re + invdet_re*val1_im;
				j -= 1;
			}
		}
		
		norm = mes_norm_inf(tmp1_re) + mes_norm_inf(tmp1_im);
		tmp1_re *= (1/norm);
		if ( l_im != 0.0 )
			tmp1_im *= (1/norm);
		tmp2_re=QQ*tmp1_re;
		if ( l_im != 0.0 )
			tmp2_im=QQ*tmp1_im;
		if ( l_im != 0.0 )
			norm = pythagoras(mes_norm2(tmp2_re),mes_norm2(tmp2_im));
		else
			norm = mes_norm2(tmp2_re);
		tmp2_re *= (1/norm);
		if ( l_im != 0.0 )
			tmp2_im *= (1/norm);
		
		if ( l_im != 0.0 ){
			mes_set_col(XX_re,i,tmp2_re);
			mes_set_col(XX_im,i,tmp2_im);
			tmp2_im *= -1.0;
			mes_set_col(XX_re,i+1,tmp2_re);
			mes_set_col(XX_im,i+1,tmp2_im);
			i += 2;
		}
		else{
			mes_set_col(XX_re,i,tmp2_re);
			mes_set_col(XX_im,i,tmp1_im);	/* zero vector */
			i += 1;
		}
	}
		
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// Hessenberg factorisations
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Hfactor -- compute Hessenberg factorisation in compact form.
//	-- factorisation performed in situ
//	-- for details of the compact form see QRfactor.c and matrix2.doc
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_Hfactor(Matrix<T>& A, Vector<T>& diag, Vector<T>& beta){
	static	Vector<double>	tmp1;
	int	k, limit;
	double beta1;

	if(A.IsDouble()==false)
		mes_warn(WARN_SINGLE_PRECISION,"mes_Hfactor");

	if ( diag.Size() < NumRows(A) - 1 || beta.Size() < NumRows(A) - 1 )
		mes_error(E_SIZES,"Hfactor");
	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SQUARE,"Hfactor");
	limit = NumRows(A);
	tmp1.Resize(NumRows(A));
	for ( k = 1; k < limit; k++ ){
		mes_get_col(A,(int)k,tmp1);
		beta1 = beta[k];
		mes_hhvec(tmp1,k+1,&beta1,tmp1,&A[k+1][k]);
		beta[k] = beta1;
		diag[k] = tmp1[k+1];
		mes_hhtrcols(A,k+1,k+1,tmp1,beta1);
		mes_hhtrrows(A,1  ,k+1,tmp1,beta1);
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// makeHQ -- construct the Hessenberg orthogonalising matrix Q;
//	-- i.e. Hess M = Q * M *Q^T
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_makeHQ(const Matrix<T>& H, Vector<T>& diag, Vector<T>& beta, Matrix<T>& Qout){
	int	i, j, limit;
	static	Vector<T>	tmp1, tmp2;

	limit = NumRows(H);
	if ( diag.Size() < limit || beta.Size() < limit )
		mes_error(E_SIZES,"makeHQ");
	if ( NumRows(H) != NumCols(H) )
		mes_error(E_SQUARE,"makeHQ");
	Qout.Resize(NumRows(H),NumRows(H));

	tmp1.Resize(NumRows(H));
	tmp2.Resize(NumRows(H));

	for ( i = 1; i <= NumRows(H); i++ ){
		/* tmp1 = i'th basis vector */
		for ( j = 1; j <= NumRows(H); j++ )
			tmp1[j] = 0.0;
		tmp1[i] = 1.0;

		/* apply H/h transforms in reverse order */
		for ( j = limit-1; j > 0; j-- ){
			mes_get_col(H,(int)j,tmp2);
			tmp2[j+1] = diag[j];
			mes_hhtrvec(tmp2,beta[j],j+1,tmp1,tmp1);
		}

		/* insert into Qout */
		mes_set_col(Qout,(int)i,tmp1);
	}

	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// makeH -- construct actual Hessenberg matrix
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_makeH(const Matrix<T>& H,Matrix<T>& Hout){
	int	i, j, limit;

	if ( NumRows(H) != NumCols(H) )
		mes_error(E_SQUARE,"makeH");
	Hout.Resize(NumRows(H),NumRows(H));
	Hout=H;

	limit = NumRows(H);
	for ( i = 2; i <= limit; i++ )
		for ( j = 1; j < i-1; j++ )
			Hout[i][j] = 0.0;

	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// Householder transformation
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// hhvec -- calulates Householder vector to eliminate all entries after the
//	i0 entry of the vector vec. It is returned as out. May be in-situ */
// Note:
//   Do not change beta to float type.
//   The stability of the Householder transformation is very sensitive
//   to the precision of beta in hhvec() and hhtrrows().
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_hhvec(Vector<T>& vec,int i0,double* beta,Vector<T>& out,T* newval){
	double	norm;

	mes_copy_offset(vec,out,i0);
	norm = mes_norm2_offset(out,i0);
	if ( norm == 0.0 ){
		*beta = 0.0;
		*newval = out[i0];
		return 1;
	}
	double abs_val=Abs(out[i0]);
	*beta = 1.0/(norm * (norm+abs_val));
	if(abs_val == 0.0){
		*newval = norm;
	}
	else{
		*newval = -mtl_sign(out[i0])*T(norm);
	}
	out[i0] -= *newval;

	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// hhtrvec -- apply Householder transformation to vector -- may be in-situ
// hh = Householder vector
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_hhtrvec(Vector<T>& hh,double beta,int i0,Vector<T>& in,Vector<T>& out){
	T	scale;
	int	i;
	if ( in.Size() != hh.Size() )
		mes_error(E_SIZES,"hhtrvec");
	if ( i0 > in.Size() )
		mes_error(E_BOUNDS,"hhtrvec");
	scale = -T(beta)*mes_in_prod_offset(hh,in,i0);
	out = in;
	for ( i=i0; i<=out.Size(); i++ )
		out[i] += scale*hh[i];
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// hhtrrows -- transform a matrix by a Householder vector by rows
//	starting at row i0 from column j0 -- in-situ
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_hhtrrows(Matrix<T>& M,int i0,int j0,Vector<T>& hh,double beta){
	T	ip, scale;
	int	i;

	if ( NumCols(M) != hh.Size() )
		mes_error(E_RANGE,"hhtrrows");
	if ( i0 > NumRows(M)+1 || j0 > NumCols(M)+1 )
		mes_error(E_BOUNDS,"hhtrrows");

	if ( beta == 0.0 )	return 1;

	/* for each row ... */
	for ( i = i0; i <= NumRows(M); i++ ){	/* compute inner product */
		int j;
		ip = T();
		for ( j = j0; j <= NumCols(M); j++ )
			ip += M[i][j]*hh[j];
		scale = -beta*ip;
		if ( mtl_iszero(scale) )
		    continue;

		/* do operation */
		for ( j = j0; j <= NumCols(M); j++ )
			M[i][j] += scale*conj(hh[j]);
	}

	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// hhtrcols -- transform a matrix by a Householder vector by columns
//	starting at row i0 from column j0 -- in-situ
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_hhtrcols(Matrix<T>& M,int i0,int j0,Vector<T>& hh,double beta){
	int	i,k;
	T scale;
	static	Vector<T>	w;

	if ( NumRows(M) != hh.Size() )
		mes_error(E_SIZES,"hhtrcols");
	if ( i0 > NumRows(M)+1 || j0 > NumCols(M)+1 )
		mes_error(E_BOUNDS,"hhtrcols");

	if ( beta == 0.0 )	return 1;

	w.Resize(NumCols(M));
	w=T();

	for ( i = i0; i <= NumRows(M); i++ ){
		if ( !mtl_iszero(hh[i]) ){
			for(k=j0;k<=NumCols(M);k++){
				w[k]+=conj(M[i][k])*hh[i];
			}
		}		
	}
	for ( i = i0; i <= NumRows(M); i++ ){
		if ( !mtl_iszero(hh[i]) ){
			for(k=j0;k<=NumCols(M);k++){
				scale = -T(beta)*hh[i];
				M[i][k]+=conj(w[k])*scale;
			}
		}
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// Givens rotations
//
// Complex Givens rotation matrix:
//   [ c   -s ]
//   [ s*   c ]
// Note that c is real and s is complex
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// givens -- returns c,s parameters for Givens rotation to
//		eliminate y in the vector [ x y ]'
///////////////////////////////////////////////////////////////////////////////
template <typename T> void	Matrix<T>::mes_givens(T x,T y,double* c,T* s){
	double norm = pythagoras(x,y);
	double inv_norm;
	if ( norm == 0.0 ){ /* identity */
		*c = 1.0;
		*s = 0.0;
	}	
	else{ // To support complex
		double absx=Abs(x);
		T xnorm;
		if(absx==0){
			xnorm=1;
		}
		else{
			xnorm = x/absx; // normalize x = x/|x|
		}
		inv_norm = 1/norm; // inv_norm = 1/||[x,y]||
		*c = inv_norm*absx;
		*s = -inv_norm*((xnorm)*conj(y)); // For real signals: c=x/norm, s=-y/norm
	}
}


///////////////////////////////////////////////////////////////////////////////
// rot_vec -- apply Givens rotation to x's i & k components
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_rot_vec(Vector<T>& x,int i,int k,double c,T s,Vector<T>& out){
	if ( i > x.Size() || k > x.Size() )
		mes_error(E_RANGE,"rot_vec");
	if(&x[1] != &out[1])
		out = x;
	T temp=c*out[i] - s*out[k]; // To support complex.
	out[k] = c*out[k] + conj(s)*out[i];
	out[i]=temp;
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// rot_rows -- premultiply mat by givens rotation described by c,s
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_rot_rows(Matrix<T>& mat,int i,int k,double c,T s,Matrix<T>& out){
	int	j;
	if ( i > NumRows(mat) || k > NumRows(mat) )
		mes_error(E_RANGE,"rot_rows");
	out = mat;
	for ( j=1; j<=NumCols(mat); j++ ){
		T temp=c*out[i][j] - s*out[k][j]; // To support complex.
		out[k][j] = c*out[k][j] + conj(s)*out[i][j];
		out[i][j] = temp;
	}

	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// rot_cols -- postmultiply mat by givens rotation described by c,s
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_rot_cols(Matrix<T>& mat,int i,int k,double c,T s,Matrix<T>& out){
	int	j;
	if ( i > NumCols(mat) || k > NumCols(mat) )
		mes_error(E_RANGE,"rot_cols");
	out = mat;
	for ( j=1; j<=NumRows(mat); j++ ){
		T temp=c*out[j][i] - conj(s)*out[j][k]; // To support complex.
		out[j][k] = c*out[j][k] + s*out[j][i];
		out[j][i] = temp;
	}
	return 1;
}


#endif


