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
// Filename: mes_zeig.h
// Revision:
//    1. Added mes_zschur_evals() and mes_zschur_vecs() to support
//       complex matrices.
// Note:
//    1. I strongly recommend the routines in this file be called 
//       with double precision arguments.
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_ZEIG_H_
#define	_MES_ZEIG_H_	

/* --------------- Eigen Analysis for Complex Matrices ----------------- */

///////////////////////////////////////////////////////////////////////////////
// mes_zschur_evals -- compute eigenvalues
//   Note: This routine assumes complex<double> type.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_zschur_evals(const Matrix<T>& TT,Vector<T>& EEVal){
	int	i, n;
	if ( NumRows(TT) != NumCols(TT) )
		mes_error(E_SQUARE,"schur_evals");
	n = NumCols(TT);
	EEVal.Resize((int)n);
	for ( i=1; i <= n; i++ ){
		EEVal[i] = TT[i][i];
	}
	return 1;
}

///////////////////////////////////////////////////////////////////////////////
// mes_zschur_vecs -- returns eigenvectors computed from the real Schur
//		decomposition of a matrix
//	-- T is the block upper triangular Schur matrix
//	-- Q is the orthognal matrix where A = Q.T.Q^T
//	-- if Q is null, the eigenvectors of T are returned
//	-- X is the matrix of column eigenvectors
// Note:
//     This routine assumes complex<double> type.
//     The precision of the eigenvector and temporary vectors is
//     critical in stability of the algorithm.
//     Must use 'double' type.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_zschur_vecs(const Matrix<T>& TT,const Matrix<T>& QQ,Matrix<T>& XX){
	int	i, j, k, limit;
	double magdet, norm, scale;
	T	l, invdet, val1, t11;	
	static Vector<T>	tmp1, tmp2;

	if ( NumRows(TT) != NumCols(TT) || NumRows(XX) != NumCols(XX) || ( NumRows(QQ) != NumCols(QQ) ) )
	    mes_error(E_SQUARE,"zschur_vecs");
	if ( NumRows(TT) != NumRows(XX) || ( NumRows(TT) != NumRows(QQ) ) )
	    mes_error(E_SIZES,"zschur_vecs");

	tmp1.Resize(NumRows(TT));
	tmp2.Resize(NumRows(TT));

	for(i=1; i<=NumRows(TT);i++){
		l=TT[i][i]; /* complex eigenvalue */
		for(j=1;j<=tmp1.Size();j++){
			tmp1[j]=T(MACH_EPS*mes_mrand());
		}
		
		/* solve (T-l.I)x = tmp1 */
		limit = i;
		for ( j = limit+1; j <= NumRows(TT); j++ )
			tmp1[j] = 0.0;

		for (j = limit;j >= 1; j-- ){
			T sum;
			t11 = TT[j][j] - l;
			magdet = real(t11)*real(t11)+imag(t11)*imag(t11);
			scale = Abs(TT[j][j]) + Abs(l);
			if ( sqrt(magdet) < MACH_EPS*scale ){
				t11 = MACH_EPS*scale;
				magdet = real(t11)*real(t11)+imag(t11)*imag(t11);
			}
			invdet = conj(t11)/magdet;
			for(sum=T(),k=j+1;k<=limit;k++) sum += tmp1[k]*TT[j][k]; 
			val1 = tmp1[j] - sum;
			tmp1[j] = invdet*val1;
		}
		
		norm = Matrix<T>::mes_norm_inf(tmp1);
		tmp1 *= (1/norm);
		tmp2=QQ*tmp1;
		norm = Matrix<T>::mes_norm2(tmp2);
		tmp2 *= (1/norm);	
		Matrix<T>::mes_set_col(XX,i,tmp2);
	}
	return 1;
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


///////////////////////////////////////////////////////////////////////////////
// Hfactor -- compute Hessenberg factorisation in compact form.
//	-- factorisation performed in situ
//	-- for details of the compact form see QRfactor.c and matrix2.doc
//   Note: This routine assumes complex<double> type.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int	Matrix<T>::mes_zHfactor(Matrix<T>& A, Vector<T>& diag){
	static	Vector<T> tmp1;
	int	k, limit;
	double beta;

	if ( diag.Size() < NumRows(A) - 1 )
		mes_error(E_SIZES,"Hfactor");
	if ( NumRows(A) != NumCols(A) )
		mes_error(E_SQUARE,"Hfactor");
	limit = NumRows(A);
	tmp1.Resize(NumRows(A));
	for ( k = 1; k < limit; k++ ){
		mes_get_col(A,k,tmp1);
		mes_hhvec(tmp1,k+1,&beta,tmp1,&A[k+1][k]);
		diag[k] = tmp1[k+1];
		mes_hhtrcols(A,k+1,k+1,tmp1,beta);
		mes_hhtrrows(A,1  ,k+1,tmp1,beta);
	}
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// zHQunpack -- unpack the compact representation of H and Q of a
//	Hessenberg factorisation
//	-- if either H or Q is NULL, then it is not unpacked
//	-- it can be in situ with HQ == H
//	-- returns HQ
// Note: This routine assumes complex<double> type.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_zHQunpack(Matrix<T>& HQ,Vector<T>& diag,Matrix<T>& Q,Matrix<T>& H)
{
	int compQ=1;
	int	i, j, limit;
	double	beta, r_ii, tmp_val;
	static	Vector<T> tmp1, tmp2;

	H.Resize(NumRows(HQ),NumCols(HQ));
	Q.Resize(NumRows(HQ),NumRows(HQ));
	if ( &HQ[1] == &Q[1] || &H[1] == &Q[1] )
	    mes_error(E_INSITU,"zHQunpack");
	limit = NumRows(HQ);
	if ( diag.Size() < limit )
		mes_error(E_SIZES,"zHQunpack");
	if ( NumRows(HQ) != NumCols(HQ) )
		mes_error(E_SQUARE,"zHQunpack");


	if ( compQ ){
	    tmp1.Resize(NumRows(H));
	    tmp2.Resize(NumRows(H));
	    
		for ( i = 1; i <= NumRows(H); i++ ){
			/* tmp1 = i'th basis vector */
			for ( j = 1; j <= NumRows(H); j++ )
				tmp1[j] = 0.0;
			tmp1[i] = 1.0;
			
			/* apply H/h transforms in reverse order */
			for ( j = limit-1; j >= 1; j-- )
			{
				mes_get_col(HQ,j,tmp2);
				r_ii = Abs(tmp2[j+1]);
				tmp2[j+1] = diag[j];
				tmp_val = (r_ii*Abs(diag[j]));
				beta = ( tmp_val == 0.0 ) ? 0.0 : 1.0/tmp_val;
				mes_hhtrvec(tmp2,beta,j+1,tmp1,tmp1);
			}
			
			/* insert into Q */
			mes_set_col(Q,i,tmp1);
		}
	}
	
	if ( NumCols(H)>0 ) {
		H = HQ;
		limit = NumRows(H);
		for ( i = 2; i <= limit; i++ )
			for ( j = 1; j < i-1; j++ )
				H[i][j] = 0.0;
	}
	
	return 1;
}


///////////////////////////////////////////////////////////////////////////////
// zschur -- computes the Schur decomposition of the complex matrix A in situ
//	-- optionally, gives Q matrix such that Q^T.A.Q is upper triangular
//	-- returns upper triangular Schur matrix
// Note: This routine assumes complex<double> type.
//       This implementation is very sensitive to the condition number of 
//       the input matrix.
///////////////////////////////////////////////////////////////////////////////
template <typename T> int Matrix<T>::mes_zschur(Matrix<T>& AA,Matrix<T>& QQ, int maxIter){
	int compQ=1;
    int		i, j, k, k_min, k_max, k_tmp, n, split, iter, isOK=1;
	double c;
	T det, discrim, lambda, lambda0, lambda1, s, sum, tmp;
	T x, y;
    static	Vector<T> diag;

	if(AA.IsReal())
		mes_error(E_NEED_COMPLEX,"mes_zschur");
	if(AA.IsDouble()==false)
		mes_warn(WARN_SINGLE_PRECISION,"mes_zschur");

    if ( NumRows(AA) != NumCols(AA) || ( NumRows(QQ) != NumCols(QQ) ) )
		mes_error(E_SQUARE,"schur");
    if ( NumRows(QQ) != NumRows(AA) )
		mes_error(E_SIZES,"schur");
    n = NumCols(AA);
    diag.Resize(NumCols(AA));
    /* compute Hessenberg form */
    isOK=Matrix<T>::mes_zHfactor(AA,diag);
	if(isOK==0) return isOK;

    /* save QQ if necessary, and make AA explicitly Hessenberg */
    isOK=Matrix<T>::mes_zHQunpack(AA,diag,QQ,AA);
	if(isOK==0) return isOK;
		
    k_min = 1;
    while ( k_min <= n ){
		
		/* find k_max to suit:
		submatrix k_min..k_max should be irreducible */
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
		
		/* now have r x r block with r >= 2:
		apply Francis QR step until block splits */
		split = 0;
		iter = 0;
		while ( ! split ){
			T	a00, a01, a10, a11;
			iter++;
			
			/* set up Wilkinson/Francis complex shift */
			/* use the smallest eigenvalue of the bottom 2 x 2 submatrix */
			k_tmp = k_max - 1;
			
			a00 = AA(k_tmp,k_tmp);
			a01 = AA(k_tmp,k_max);
			a10 = AA(k_max,k_tmp);
			a11 = AA(k_max,k_max);
			tmp=0.5*(a00-a11);
			discrim=sqrt(tmp*tmp+a01*a10);
			sum=0.5*(a00+a11);
			lambda0=sum+discrim;
			lambda1=sum-discrim;
			det=a00*a11-a01*a10;
			if(Abs(lambda0) > Abs(lambda1))
				lambda=det/lambda0;
			else
				lambda=det/lambda1;

			/* perturb shift if convergence is slow */
			if((iter%10)==0){
				lambda += iter*0.02;
			}

			/* set up Householder transformations */
			k_tmp = k_min + 1;

			x=AA[k_min][k_min]-lambda;
			y=AA[k_min+1][k_min];

			/* use Givens' rotations to "chase" off-Hessenberg entry */	
			for ( k = k_min; k <= k_max-1; k++ ){
				mes_givens(x,y,&c,&s);
				mes_rot_cols(AA,k,k+1,c,s,AA);
				mes_rot_rows(AA,k,k+1,c,s,AA);
				if ( compQ )
					mes_rot_cols(QQ,k,k+1,c,s,QQ);

				/* zero thing that should be zero */
				if(k>k_min)
					AA[k+1][k-1]=0;

				/* get next entry to chase along sub-diagonal */
				x = AA(k+1,k);
				if ( k <= k_max - 2 )
					y = AA(k+2,k);
				else
					y = 0.0;
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
				cerr << "ERROR: Too many iterations in mes_zschur().\n";
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


#endif
