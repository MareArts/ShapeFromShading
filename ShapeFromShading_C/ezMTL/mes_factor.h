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
// Filename: mes_factor.h
// Revision:
//    1. Revised to support complex matrix factorization.
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_FACTOR_H_
#define	_MES_FACTOR_H_	

/**************************************************************************
**
** Copyright (C) 1993 David E. Steward & Zbigniew Leyk, all rights reserved.
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

/* Most matrix factorisation routines are in-situ unless otherwise specified */
/* Matrix factorisation routines to work with the other matrix files.     */
/* Cholesky decomposition: A=R'*R                                         */
/* LU factorization: A=L*U                                                */
/* QR factorization: A=Q*R                                                */

///////////////////////////////////////////////////////////////////////////////
// CHfactor -- Cholesky L.L' factorisation of A in-situ
// Computes the Cholesky decompostion for a positive-definite symmetric matrix.
// A = L * L^T
// The L is returned in the lower triangle of A.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_CHfactor(Matrix<T>& A){
	int	i, j, k, n;
	T  sum;
	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SQUARE,"CHfactor");
	n = NumCols(A);
	for ( k=1; k<=n; k++ ){
		double norm=mes_norm_inf(A.Col(k));
		double eps=norm*Abs(mtl_numeric_limits<T>::epsilon());
		/* do diagonal element */
		for ( sum = A[k][k], j=1; j<k; j++ ) sum -= A[k][j]*conj(A[k][j]);
		if ( Abs(sum) < eps || Abs(imag(sum)) > eps ){
			mes_error(E_POSDEF,"CHfactor");
			return 0;
		}
		else{
			sum=(sum+conj(sum))/T(2); // remove imaginary part
		}
		A[k][k] = sqrt(sum);

		/* set values of column k */
		for ( i=k+1; i<=n; i++ ){
			for ( sum = A[i][k], j=1; j<k; j++ ) sum -= A[i][j]*conj(A[k][j]);
			A[i][j] = sum/A[k][k];
			A[j][i] = conj(A[i][j]);
		}
	}
	return 1;
}


/* CHsolve -- given a CHolesky factorization in A, solve A*x=b */
template <typename T> int	Matrix<T>::mes_CHsolve(const Matrix<T>& A,const Vector<T>& b,Vector<T>& x){
	if ( NumRows(A) != NumCols(A) || NumCols(A) != b.Size() )
		mes_error(E_SIZES,"CHsolve");
	x.Resize(b.Size());
	mes_Lsolve(A,b,x,0.0);
	mes_Usolve(A,x,x,0.0);

	return 1;
}


/* LDLfactor -- L.D.L' factorisation of A in-situ */
template <typename T> int	Matrix<T>::mes_LDLfactor(Matrix<T>& A){
	int	i, k, n, p;
	double d, sum;
	static Vector<T>	r;

	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SQUARE,"LDLfactor");
	n = NumCols(A);
	r.Resize(n);
	for ( k = 1; k <= n; k++ ){
		sum = 0.0;
		for ( p = 1; p < k; p++ ){
		    r[p] = A[p][p]*A[k][p];
		    sum += r[p]*A[k][p];
		}
		d = A[k][k] -= sum;

		if ( d < mtl_numeric_limits<double>::epsilon() )
		    mes_error(E_SING,"LDLfactor");
		for ( i = k+1; i <= n; i++ ){
		    sum = 0.0;
		    for ( p = 1; p < k; p++ )
			sum += A[i][p]*r[p];
		    A[i][k] = (A[i][k] - sum)/d;
		}
	}

	return 1;
}


/* LDLsolve -- given an LDL factorisation in A, solve Ax=b */
template <typename T> int	Matrix<T>::mes_LDLsolve(const Matrix<T>& LDL,const Vector<T>& b,Vector<T>& x){
	if ( NumRows(LDL) != NumCols(LDL) )
		mes_error(E_SQUARE,"LDLsolve");
	if ( NumRows(LDL) != b.Size() )
		mes_error(E_SIZES,"LDLsolve");
	x.Resize(b.Size());
	mes_Lsolve(LDL,b,x,1.0);
	mes_Dsolve(LDL,x,x);
	mes_LTsolve(LDL,x,x,1.0);
	return 1;
}


/* LUfactor -- gaussian elimination with scaled partial pivoting
		-- Note: returns LU matrix which is A */
template <typename T> int	Matrix<T>::mes_LUfactor(Matrix<T>& A,Vector<int>& pivot){
	int	i, j, k, k_max, m, n;
	int	i_max;
	double	max1, temp, tiny;
	static	Vector<double> scale;
	T atemp;

	pivot.Resize(NumRows(A));
	if ( pivot.Size() != NumRows(A) )
		mes_error(E_SIZES,"LUfactor");
	m = NumRows(A);	n = NumCols(A);
	scale.Resize(NumRows(A));

	tiny = 10.0/HUGE_VAL;

	/* initialise pivot with identity permutation */
	for ( i=1; i<=m; i++ )
		pivot[i] = i;

	/* set scale parameters */
	for ( i=1; i<=m; i++ ){
		max1 = 0.0;
		for ( j=1; j<=n; j++ ){
			temp = Abs(A[i][j]);
			max1 = local_max(max1,temp);
		}
		scale[i] = max1;
	}

	/* main loop */
	k_max = local_min(m,n);
	for ( k=1; k<=k_max; k++ ){
	    /* find best pivot row */
	    max1 = 0.0;	i_max = 0;
	    for ( i=k; i<=m; i++ )
		if ( Abs(scale[i]) >= tiny*Abs(A[i][k]) ){
		    temp = Abs(A[i][k])/scale[i];
		    if ( temp > max1 ){
				max1 = temp;	i_max = i;
			}
		}
	    
	    /* if no pivot then ignore column k... */
	    if ( i_max == 0 ){
		/* set pivot entry A[k][k] exactly to zero,
			rather than just "small" */
			A[k][k] = 0.0;
			continue;
		}
		
		/* do we pivot ? */
		if ( i_max != k ){	/* yes we do... */
			mes_px_transp(pivot,i_max,k);
			for ( j=1; j<=n; j++ ){
				atemp = A[i_max][j];
				A[i_max][j] = A[k][j];
				A[k][j] = atemp;
			}
		}
		
		/* row operations */
		for ( i=k+1; i<=m; i++ ){	/* for each row do... */
			/* Note: divide by zero should never happen */
			atemp = A[i][k] = A[i][k]/A[k][k];
			if ( k+1 <= n ){
				for ( j=k+1; j<=n; j++ )
					A[i][j] -= atemp*A[k][j];
			}
		}
	    
	}
	return 1;
}


/* LUsolve -- given an LU factorisation in A, solve A*x=b */
template <typename T> int	Matrix<T>::mes_LUsolve(Matrix<T>& A,Vector<int>& pivot,const Vector<T>& b,Vector<T>& x){
	if ( NumRows(A) != NumCols(A) || NumCols(A) != b.Size() )
		mes_error(E_SIZES,"LUsolve");
	Vector<T> btmp=b;
	x.Resize(b.Size());
	mes_px_vec(pivot,btmp,x);	/* x := P.b */
	mes_Lsolve(A,x,x,1.0);	/* implicit diagonal = 1 */
	mes_Usolve(A,x,x,0.0);	/* explicit diagonal */
	return 1;
}


/* LUTsolve -- given an LU factorisation in A, solve A^T*x=b */
template <typename T> int	Matrix<T>::mes_LUTsolve(Matrix<T>& LU,Vector<int>& pivot,const Vector<T>& b,Vector<T>& x)
{
	if ( NumRows(LU) != NumCols(LU) || NumCols(LU) != b.Size() )
		mes_error(E_SIZES,"LUTsolve");
	x = b;
	mes_UTsolve(LU,x,x,0.0);	/* explicit diagonal */
	mes_LTsolve(LU,x,x,1.0);	/* implicit diagonal = 1 */
	mes_pxinv_vec(pivot,x,x);	/* x := P^T.tmp */
	return 1;
}


/* m_inverse -- returns inverse of A, provided A is not too rank deficient
	-- uses LU factorisation */
template <typename T> int	Matrix<T>::mes_inverse(const Matrix<T>& A,Matrix<T>& out)
{
	int	i;
	static Vector<T> tmp, tmp2;
	static Matrix<T> A_cp;
	static Vector<int> pivot;

	if ( NumRows(A) != NumCols(A) )
	    mes_error(E_SQUARE,"m_inverse");
	if ( NumRows(out) != NumRows(A) || NumCols(out) != NumCols(A) )
	    out.Resize(NumRows(A),NumCols(A));

	A_cp = A;
	tmp.Resize(NumRows(A));
	tmp2.Resize(NumRows(A));
	pivot.Resize(NumRows(A));
	mes_LUfactor(A_cp,pivot);
	for ( i = 1; i <= NumCols(A); i++ ){
	    tmp = T();
	    tmp[i] = 1.0;
	    mes_LUsolve(A_cp,pivot,tmp,tmp2);
	    mes_set_col(out,i,tmp2);
	}
	return 1;
}


/* LUcondest -- returns an estimate of the condition number of LU given the
	LU factorisation in compact form */
template <typename T> double	Matrix<T>::mes_LUcondest(Matrix<T>& LU,Vector<int>& pivot)
{
    static	Vector<T> y, z;
    double	cond_est, L_norm, U_norm, norm;
	T sum;
    int		i, j, n;

    if ( NumRows(LU) != NumCols(LU) )
	mes_error(E_SQUARE,"LUcondest");
    if ( NumCols(LU) != pivot.Size() )
	mes_error(E_SIZES,"LUcondest");

    n = NumCols(LU);
    y.Resize(n);
    z.Resize(n);
    for ( i = 1; i <= n; i++ ){
		sum = T();
		for ( j = 1; j < i; j++ )
			sum -= LU[j][i]*y[j];
		sum += mtl_sign(sum); // To support complex
		if ( mtl_iszero(LU[i][i]) )
			return HUGE_VAL;
		y[i] = sum / LU[i][i];
    }
	
    if(mes_LTsolve(LU,y,y,1.0) == 0)
		return HUGE_VAL;
	if(mes_LUsolve(LU,pivot,y,z) == 0)
		return HUGE_VAL;

    /* now estimate norm of A (even though it is not directly available) */
    /* actually computes ||L||_inf.||U||_inf */
    U_norm = 0.0;
    for ( i = 1; i <= n; i++ ){
		norm = 0.0;
		for ( j = i; j <= n; j++ )
			norm += Abs(LU[i][j]);
		if ( norm > U_norm )
			U_norm = norm;
    }
    L_norm = 0.0;
    for ( i = 1; i <= n; i++ ){
		norm = 1.0;
		for ( j = 1; j < i; j++ )
			norm += Abs(LU[i][j]);
		if ( norm > L_norm )
			L_norm = norm;
    }
	
    cond_est = U_norm*L_norm*mes_norm_inf(z)/mes_norm_inf(y);	
    return cond_est;
}


/* Note: The usual representation of a Householder transformation is taken
   to be:
   P = I - beta.u.uT
   where beta = 2/(uT.u) and u is called the Householder vector
*/

/* QRfactor -- forms the QR factorisation of A -- factorisation stored in
   compact form as described above ( not quite standard format ) */
template <typename T> int	Matrix<T>::mes_QRfactor(Matrix<T>& A,Vector<T>& diag)
{
    int	k,limit;
    double	beta;
    static	Vector<T> tmp1;
    
    limit = local_min(NumRows(A),NumCols(A));
	diag.Resize(limit);
    if ( diag.Size() < limit )
	mes_error(E_SIZES,"QRfactor");
    
    tmp1.Resize(NumRows(A));
    
    for ( k=1; k<=limit; k++ ){
		/* get H/holder vector for the k-th column */
		mes_get_col(A,k,tmp1);
		mes_hhvec(tmp1,k,&beta,tmp1,&A[k][k]);
		diag[k] = tmp1[k];
		
		/* apply H/holder vector to remaining columns */
		mes_hhtrcols(A,k,k+1,tmp1,beta);
    }

    return 1;
}


/* QRCPfactor -- forms the QR factorisation of A with column pivoting
   -- factorisation stored in compact form as described above
   ( not quite standard format )				*/
template <typename T> int	Matrix<T>::mes_QRCPfactor(Matrix<T>& A,Vector<T>& diag,Vector<int>& px)
{
    int	i, i_max, j, k, limit;
    static	Vector<T> gamma, tmp1, tmp2;
    double	beta, maxgamma, sum, tmp;
    
    limit = local_min(NumRows(A),NumCols(A));
	diag.Resize(limit);
	px.Resize(NumCols(A));
    if ( diag.Size() < limit || px.Size() != NumCols(A) )
	mes_error(E_SIZES,"QRCPfactor");
    
    tmp1.Resize(NumRows(A));
    tmp2.Resize(NumRows(A));
    gamma.Resize(NumCols(A));
    
    /* initialise gamma and px */
    for ( j=1; j<=NumCols(A); j++ ){
		px[j] = j;
		sum = 0.0;
		for ( i=1; i<=NumRows(A); i++ )
			sum += (A[i][j]*A[i][j]);
		gamma[j] = sum;
    }
    
    for ( k=1; k<=limit; k++ ){
		/* find "best" column to use */
		i_max = k;	maxgamma = gamma[k];
		for ( i=k+1; i<=NumCols(A); i++ ){
		/* Loop invariant:maxgamma=gamma[i_max]
	       >=gamma[l];l=k,...,i-1 */
		   if ( gamma[i] > maxgamma ){
			   maxgamma = gamma[i]; i_max = i;
		   }
		}
		   
		/* swap columns if necessary */
		if ( i_max != k ){
			/* swap gamma values */
			tmp = gamma[k];
			gamma[k] = gamma[i_max];
			gamma[i_max] = tmp;
			
			/* update column permutation */
			mes_px_transp(px,k,i_max);
			
			/* swap columns of A */
			for ( i=1; i<=NumRows(A); i++ ){
				tmp = A[i][k];
				A[i][k] = A[i][i_max];
				A[i][i_max] = tmp;
			}
		}
		
		/* get H/holder vector for the k-th column */
		mes_get_col(A,k,tmp1);
		mes_hhvec(tmp1,k,&beta,tmp1,&A[k][k]);
		diag[k] = tmp1[k];
		
		/* apply H/holder vector to remaining columns */
		mes_hhtrcols(A,k,k+1,tmp1,beta);
		
		/* update gamma values */
		for ( j=k+1; j<=NumCols(A); j++ )
			gamma[j] -= (A[k][j]*A[k][j]);
    }	
    return 1;
}


/* Qsolve -- solves Qx = b, Q is an orthogonal matrix stored in compact
   form a la QRfactor() -- may be in-situ */
template <typename T> int	Matrix<T>::mes_Qsolve_(Matrix<T>& QR,Vector<T>& diag,const Vector<T>& b,Vector<T>& x,Vector<T>& tmp)
{
    int		k, limit;
    double	beta, r_ii, tmp_val;
    
    limit = local_min(NumRows(QR),NumCols(QR));
    if ( diag.Size() < limit || b.Size() != NumRows(QR) )
		mes_error(E_SIZES,"Qsolve_");
    x.Resize(NumRows(QR));
    tmp.Resize(NumRows(QR));
    
    /* apply H/holder transforms in normal order */
    x = b;
    for ( k = 1 ; k <= limit ; k++ ){
		mes_get_col(QR,k,tmp);
		r_ii = Abs(tmp[k]);
		tmp[k] = diag[k];
		tmp_val = (r_ii*Abs(diag[k]));
		beta = ( tmp_val == 0.0 ) ? 0.0 : 1.0/tmp_val;
		mes_hhtrvec(tmp,beta,k,x,x);
    }
    
    return 1;
}


/* makeQ -- constructs orthogonal matrix from Householder vectors stored in
   compact QR form */
template <typename T> int	Matrix<T>::mes_makeQ(Matrix<T>& QR,Vector<T>& diag,Matrix<T>& Qout)
{
    static	Vector<T> tmp1,tmp2;
    int	i, limit;
    double	beta, r_ii, tmp_val;
    int	j;
    
    limit = local_min(NumRows(QR),NumCols(QR));
    if ( diag.Size() < limit )
		mes_error(E_SIZES,"makeQ");
    if ( NumRows(Qout) != NumRows(QR) || NumCols(Qout) != NumRows(QR) )
		Qout.Resize(NumRows(QR),NumRows(QR));
    
    tmp1.Resize(NumRows(QR));	/* contains basis vec & columns of Q */
    tmp2.Resize(NumRows(QR));	/* contains H/holder vectors */
    
    for ( i=1; i<=NumRows(QR) ; i++ ){	/* get i-th column of Q */
		/* set up tmp1 as i-th basis vector */
		for ( j=1; j<=NumRows(QR) ; j++ )
			tmp1[j] = 0.0;
		tmp1[i] = 1.0;
		
		/* apply H/h transforms in reverse order */
		for ( j=limit; j>=1; j-- ){
			mes_get_col(QR,j,tmp2);
			r_ii = Abs(tmp2[j]);
			tmp2[j] = diag[j];
			tmp_val = (r_ii*Abs(diag[j]));
			beta = ( tmp_val == 0.0 ) ? 0.0 : 1.0/tmp_val;
			mes_hhtrvec(tmp2,beta,j,tmp1,tmp1);
		}
		
		/* insert into Q */
		mes_set_col(Qout,i,tmp1);
    }
	
    return 1;
}


/* makeR -- constructs upper triangular matrix from QR (compact form)
   -- may be in-situ (all it does is zero the lower 1/2) */
template <typename T> int	Matrix<T>::mes_makeR(Matrix<T>& QR,Matrix<T>& Rout)
{
    int	i,j;
    
    Rout = QR;
    for ( i=2; i<=NumRows(QR); i++ )
	for ( j=1; j<=NumCols(QR) && j<i; j++ )
	    Rout[i][j] = 0.0;
    
    return 1;
}


/* QRsolve -- solves the system Q.R.x=b where Q & R are stored in compact form
   -- returns x, which is created if necessary */
template <typename T> int	Matrix<T>::mes_QRsolve(Matrix<T>& QR,Vector<T>& diag,const Vector<T>& b,Vector<T>& x)
{
    int	limit;
    static	Vector<T> tmp;
    
    limit = local_min(NumRows(QR),NumCols(QR));
    if ( diag.Size() < limit || b.Size() != NumRows(QR) )
		mes_error(E_SIZES,"QRsolve");
    tmp.Resize(limit);	
    x.Resize(NumCols(QR));
    mes_Qsolve_(QR,diag,b,x,tmp);
    mes_Usolve(QR,x,x,0.0);
    x.Resize(NumCols(QR));	
    return 1;
}


/* QRCPsolve -- solves A.x = b where A is factored by QRCPfactor()
   -- assumes that A is in the compact factored form */
template <typename T> int	Matrix<T>::mes_QRCPsolve(Matrix<T>& QR,Vector<T>& diag,Vector<int>& pivot,const Vector<T>& b,Vector<T>& x)
{
    static	Vector<T> tmp;
    
    if ( (NumRows(QR) > diag.Size() && NumCols(QR) > diag.Size()) || NumCols(QR) != pivot.Size() )
		mes_error(E_SIZES,"QRCPsolve");
    
    mes_QRsolve(QR,diag /* , beta */ ,b,tmp);
    mes_pxinv_vec(pivot,tmp,x);
	
    return 1;
}


/* Umlt -- compute out = upper_triang(U).x
	-- may be in situ */
template <typename T> int	Matrix<T>::mes_Umlt(const Matrix<T>& U,const Vector<T>& x,Vector<T>& out)
{
    int		i, limit;

    limit = local_min(NumRows(U),NumCols(U));
    if ( limit != x.Size() )
		mes_error(E_SIZES,"Umlt");
    if ( out.Size() != limit )
		out.Resize(limit);
	
    for ( i = 1; i <= limit; i++ ){
		T sum = T();
		for(int j=i; j<=limit; j++)
			sum += x[j]*U[i][j];
		out[i] = sum;
	}
    return 1;
}


/* UTmlt -- returns out = upper_triang(U)^T.x */
template <typename T> int	Matrix<T>::mes_UTmlt(const Matrix<T>& U,const Vector<T>& x,Vector<T>& out)
{
    T	sum;
    int		i, j, limit;
	
    limit = local_min(NumRows(U),NumCols(U));
    if ( out.Size() != limit )
		out.Resize(limit);
	
    for ( i = limit; i >= 1; i-- ){
		sum = 0.0;
		for ( j = 1; j <= i; j++ )
			sum += conj(U[j][i])*x[j];
		out[i] = sum;
    }
    return 1;
}


/* QRTsolve -- solve A^T.sc = c where the QR factors of A are stored in
	compact form
	-- returns sc
	-- original due to Mike Osborne modified Wed 09th Dec 1992 */
template <typename T> int	Matrix<T>::mes_QRTsolve(const Matrix<T>& A,Vector<T>& diag,Vector<T>& c,Vector<T>& sc)
{
    int		i, j, k, n, p;
    double	beta, r_ii, s, tmp_val;

    if ( diag.Size() < local_min(NumRows(A),NumCols(A)) )
		mes_error(E_SIZES,"QRTsolve");
    sc.Resize(sc,NumRows(A));
    n = sc.Size();
    p = c.Size();
    if ( n == p )
		k = p-1;
    else
		k = p;
    sc=0;
    sc[1] = c[1]/A[1][1];
    if ( n ==  1)
		return sc;
    if ( p > 2){
		for ( i = 2; i <= p; i++ ){
			s = 0.0;
			for ( j = 1; j < i; j++ )
				s += A[j][i]*sc[j];
			if ( A[i][i] == 0.0 )
				mes_error(E_SING,"QRTsolve");
			sc[i]=(c[i]-s)/A[i][i];
		}
    }
    for (i = k; i >= 1; i--){
		s = diag[i]*sc[i];
		for ( j = i+1; j <= n; j++ )
			s += A[j][i]*sc[j];
		r_ii = Abs(A[i][i]);
		tmp_val = (r_ii*Abs(diag[i]));
		beta = ( tmp_val == 0.0 ) ? 0.0 : 1.0/tmp_val;
		tmp_val = beta*s;
		sc[i] -= tmp_val*diag[i];
		for ( j = i+1; j <= n; j++ )
			sc[j] -= tmp_val*A[j][i];
    }
	
    return 1;
}


/* QRcondest -- returns an estimate of the 2-norm condition number of the
		matrix factorised by QRfactor() or QRCPfactor()
	-- note that as Q does not affect the 2-norm condition number,
		it is not necessary to pass the diag, beta (or pivot) vectors
	-- generates a lower bound on the true condition number
	-- if the matrix is exactly singular, HUGE_VAL is returned
	-- note that QRcondest() is likely to be more reliable for
		matrices factored using QRCPfactor() */
template <typename T> double	Matrix<T>::mes_QRcondest(const Matrix<T>& QR, int p)
{
    static	Vector<T>	y;
    double	norm1, norm2, tmp1, tmp2;
    int		i, j, limit;
	T sum;

	if(p!=2){
		cerr << "ERROR: QRcondest() was not implemented for p!=2.\n";
		assert(0);
		exit(1);
		return HUGE_VAL;
	}
    limit = local_min(NumRows(QR),NumCols(QR));
    for ( i = 1; i <= limit; i++ ){
		if ( mtl_iszero(QR[i][i]) )
			return HUGE_VAL;
	}
	
	y.Resize(limit);
	/* use the trick for getting a unit vector y with ||R.y||_inf small
	from the LU condition estimator */
	for ( i = 1; i <= limit; i++ ){
		sum = T();
		for ( j = 1; j < i; j++ )
			sum -= QR[j][i]*y[j];
		sum += mtl_sign(sum); // to support complex
		y[i] = sum / QR[i][i];
	}
	mes_UTmlt(QR,y,y);
	
	/* now apply inverse power method to R^T.R */
	for ( i = 1; i <= 3; i++ ){
		tmp1 = Norm(y,p);
		y *= (1/tmp1);
		mes_UTsolve(QR,y,y,0.0);
		tmp2 = Norm(y,p);
		y *= (1/tmp2);
		mes_Usolve(QR,y,y,0.0);
	}
	/* now compute approximation for ||R^{-1}||_2 */
	norm1 = sqrt(tmp1)*sqrt(tmp2);
	
	/* now use complementary approach to compute approximation to ||R||_2 */
	for ( i = limit; i >= 1; i-- ){
		sum = T();
		for ( j = i+1; j <= limit; j++ )
			sum += QR[i][j]*y[j];
		y[i]=mtl_sign(sum)*T(Abs(QR[i][i])); // To support complex
	}
	
	/* now apply power method to R^T.R */
	for ( i = 1; i <= 3; i++ ){
		tmp1 = Norm(y,p);
		y *= (1/tmp1);
		mes_Umlt(QR,y,y);
		tmp2 = Norm(y,p);
		y *= (1/tmp2);
		mes_UTmlt(QR,y,y);
	}
	norm2 = sqrt(tmp1)*sqrt(tmp2);
	
	return norm1*norm2;
}


/* interchange -- a row/column swap routine */
template <typename T> void Matrix<T>::mes_interchange(Matrix<T>& A,int i,int j)
{
	T	tmp;
	int	k, n;

	n = NumCols(A);
	if ( i == j )
		return;
	if ( i > j ){	
		k = i;	i = j;	j = k;	
	}
	for ( k = 0; k < i; k++ ){
		tmp = A[k][i];
		A[k][i] = A[k][j];
		A[k][j] = tmp;
	}
	for ( k = j+1; k < n; k++ ){
		tmp = A[j][k];
		A[j][k] = A[i][k];
		A[i][k] = tmp;
	}
	for ( k = i+1; k < j; k++ ){
		tmp = A[k][j];
		A[k][j] = A[i][k];
		A[i][k] = tmp;
	}
	tmp = A[i][i];
	A[i][i] = A[j][j];
	A[j][j] = tmp;
}


/* BKPfactor -- Bunch-Kaufman-Parlett factorisation of A in-situ
	-- A is factored into the form P'AP = MDM' where 
	P is a permutation matrix, M lower triangular and D is block
	diagonal with blocks of size 1 or 2
	-- P is stored in pivot; blocks[i]==i iff D[i][i] is a block */
template <typename T> int Matrix<T>::mes_BKPfactor(Matrix<T>& A,Vector<int>& pivot,Vector<int>& blocks)
{
	const double mes_alpha = 0.6403882032022076; /* = (1+sqrt(17))/8 */
	int	i, j, k, n, onebyone, r;
	double	aii, aip1, aip1i, lambda, sigma, tmp;
	double	det, s, t;

	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SQUARE,"BKPfactor");
	pivot.Resize(NumRows(A));
	blocks.Resize(NumRows(A));
	if ( NumRows(A) != pivot.Size() || pivot.Size() != blocks.Size() )
		mes_error(E_SIZES,"BKPfactor");

	n = NumCols(A);
	pivot.SetPermIdentity();
	blocks.SetPermIdentity();
	
	for ( i = 1; i <= n; i = onebyone ? i+1 : i+2 ){
		aii = Abs(A(i,i));
		lambda = 0.0;	r = (i+1 < n) ? i+1 : i;
		for ( k = i+1; k <= n; k++ ){
			tmp = Abs(A(i,k));
			if ( tmp >= lambda ){
				lambda = tmp;
				r = k;
			}
		}
		/* determine if 1x1 or 2x2 block, and do pivoting if needed */
		if ( aii >= mes_alpha*lambda ){
			onebyone = true;
			goto dopivot;
		}
		/* compute sigma */
		sigma = 0.0;
		for ( k = i; k <= n; k++ ){
			if ( k == r )
				continue;
			tmp = ( k > r ) ? Abs(A(r,k)) :
			Abs(A(k,r));
			if ( tmp > sigma )
				sigma = tmp;
		}
		if ( aii*sigma >= mes_alpha*(lambda*lambda) )
			onebyone = true;
		else if ( Abs(A(r,r)) >= mes_alpha*sigma ){
			mes_interchange(A,i,r);
			mes_px_transp(pivot,i,r);
			onebyone = true;
		}
		else{
			mes_interchange(A,i+1,r);
			mes_px_transp(pivot,i+1,r);
			mes_px_transp(blocks,i,i+1);
			onebyone = false;
		}
		
dopivot:
		if ( onebyone ){   /* do one by one block */
			if ( A(i,i) != 0.0 ){
				aii = A(i,i);
				for ( j = i+1; j <= n; j++ ){
					tmp = A(i,j)/aii;
					for ( k = j; k <= n; k++ )
						A[j][k] -= tmp*A(i,k);
					A[i][j]=tmp;
				}
			}
		}
		else{ /* onebyone == false */
		    /* do two by two block */
			det = A(i,i)*A(i+1,i+1)-(A(i,i+1)*A(i,i+1));
			/* Must have det < 0 */
			aip1i = A(i,i+1)/det;
			aii = A(i,i)/det;
			aip1 = A(i+1,i+1)/det;
			for ( j = i+2; j <= n; j++ ){
				s = - aip1i*A(i+1,j) + aip1*A(i,j);
				t = - aip1i*A(i,j) + aii*A(i+1,j);
				for ( k = j; k <= n; k++ )
					A[j][k] -= A(i,k)*s + A(i+1,k)*t;
				A[i][j]=s;
				A[i+1][j]=t;
			}
		}
	}
	
	/* set lower triangular half */
	for ( i = 1; i <= NumRows(A); i++ ){
		for ( j = 1; j < i; j++ )
			A[i][j]=A[j][i];
	}
	
	return 1;
}


/* BKPsolve -- solves A.x = b where A has been factored a la BKPfactor()
	-- returns x, which is created if NULL */
template <typename T> int Matrix<T>::mes_BKPsolve(const Matrix<T>& A,Vector<int>& pivot,Vector<int>& block,const Vector<T>& b,Vector<T>& x)
{
	static Vector<T> tmp;	/* dummy storage needed */
	int	i, j, n, onebyone;
	double	a11, a12, a22, b1, b2, det, sum, tmp_diag;

	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SQUARE,"BKPsolve");
	n = NumCols(A);
	if ( b.Size() != n || pivot.Size() != n || block.Size() != n )
		mes_error(E_SIZES,"BKPsolve");
	x.Resize(n);
	tmp.Resize(n);

	Vector<T> btmp=b;
	mes_px_vec(pivot,btmp,tmp);
	/* solve for lower triangular part */
	for ( i = 1; i <= n; i++ ){
		sum = tmp(i);
		if ( block[i] < i ){
			for ( j = 1; j < i-1; j++ )
				sum -= A(i,j)*tmp(j);
		}
		else{
			for ( j = 1; j < i; j++ )
				sum -= A(i,j)*tmp(j);
		}
		tmp[i]=sum;
	}

	/* solve for diagonal part */
	for ( i = 1; i <= n; i = onebyone ? i+1 : i+2 ){
		onebyone = ( block[i] == i );
		if ( onebyone ){
			tmp_diag = A(i,i);
			if ( tmp_diag == 0.0 )
				mes_error(E_SING,"BKPsolve");
			tmp[i]=tmp[i]/tmp_diag;
		}
		else{
			a11 = A(i,i);
			a22 = A(i+1,i+1);
			a12 = A(i+1,i);
			b1 = tmp(i);	b2 = tmp(i+1);
			det = a11*a22-a12*a12;	/* < 0 : see BKPfactor() */
			if ( det == 0.0 )
				mes_error(E_SING,"BKPsolve");
			det = 1/det;
			tmp[i]=det*(a22*b1-a12*b2);
			tmp[i+1]=det*(a11*b2-a12*b1);
		}
	}

	/* solve for transpose of lower traingular part */
	for ( i = n; i >= 1; i-- ){	/* use symmetry of factored form to get stride 1 */
		sum = tmp(i);
		if ( block[i] > i ){
			for ( j = i+2; j <= n; j++ )
				sum -= A(i,j)*tmp(j);
		}
		else{
			for ( j = i+1; j <= n; j++ )
				sum -= A(i,j)*tmp(j);
		}
		tmp[i]=sum;
	}

	/* and do final permutation */
	mes_pxinv_vec(pivot,tmp,x);

	return 1;
}


#endif
