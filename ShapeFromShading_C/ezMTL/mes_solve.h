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
// Filename: mes_solve.h
// Revision:
//    1. Revised to support complex matrices.
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_SOLVE_H_
#define	_MES_SOLVE_H_	

/* Matrix factorisation routines to work with the other matrix files. */
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


/* Usolve -- back substitution with optional over-riding diagonal
		-- can be in-situ but doesn't need to be */
template <typename T> int Matrix<T>::mes_Usolve(const Matrix<T>& mat,const Vector<T>& b,Vector<T>& out,double diag)
{
	int	dim, j;
	int	i, i_lim;
	double tiny;
	T	sum;

	dim = local_min(NumRows(mat),NumCols(mat));
	if ( b.Size() < dim )
		mes_error(E_SIZES,"Usolve");
	if ( out.Size() < dim )
		out.Resize(NumCols(mat));

	tiny = 10.0/HUGE_VAL;

	for ( i=dim; i>=1; i-- ){
		if ( !mtl_iszero(b[i]) )
		    break;
		else
		    out[i] = 0.0;
	}
	i_lim = i;

	for ( ; i>=1; i-- ){
		sum = b[i];
		for ( j=i+1; j<=i_lim; j++ )
			sum -= mat[i][j]*out[j];
		if ( diag==0.0 ){
			if ( Abs(mat[i][i]) <= tiny*Abs(sum) )
				mes_error(E_SING,"Usolve");
			else
				out[i] = sum/mat[i][i];
		}
		else
			out[i] = sum/T(diag);
	}

	return 1;
}


/* Lsolve -- forward elimination with (optional) default diagonal value */
template <typename T> int Matrix<T>::mes_Lsolve(const Matrix<T>& mat,const Vector<T>& b,Vector<T>& out,double diag)
{
	int	dim, i, i_lim, j;
	double	tiny;
	T sum;

	dim = local_min(NumRows(mat),NumCols(mat));
	if ( b.Size() < dim )
		mes_error(E_SIZES,"Lsolve");
	if ( out.Size() < dim )
		out.Resize(NumCols(mat));

	for ( i=1; i<=dim; i++ ){
		if ( !mtl_iszero(b[i]) )
		    break;
		else
		    out[i] = 0.0;
	}
	i_lim = i;

	tiny = 10.0/HUGE_VAL;

	for ( ; i<=dim; i++ ){
		sum = b[i];
		for ( j=i_lim; j<i; j++ )
			sum -= mat[i][j]*out[j];
		if ( diag==0.0 ){
			if ( Abs(mat[i][i]) <= tiny*Abs(sum) )
				mes_error(E_SING,"Lsolve");
			else
				out[i] = sum/mat[i][i];
		}
		else
			out[i] = sum/T(diag);
	}

	return 1;
}


/* UTsolve -- forward elimination with (optional) default diagonal value
		using UPPER triangular part of matrix */
template <typename T> int Matrix<T>::mes_UTsolve(const Matrix<T>& U,const Vector<T>& b,Vector<T>& out,double diag)
{
    int	dim, i, i_lim, j;
    double	invdiag, tiny;
	T tmp;
    
    dim = local_min(NumRows(U),NumCols(U));
    if ( b.Size() < dim )
	mes_error(E_SIZES,"UTsolve");
    out.Resize(NumCols(U));

    tiny = 10.0/HUGE_VAL;

    for ( i=1; i<=dim; i++ ){
		if ( !mtl_iszero(b[i]) )
			break;
		else
			out[i] = 0.0;
	}

    i_lim = i;
    if ( &b[1] != &out[1] ){
		out = T();
		for(j=i_lim;j<=dim;j++)
			out[j]=b[j];
    }
	
    if ( diag == 0.0 ){
		for ( ; i<=dim; i++ ){
			tmp = conj(U[i][i]);
			if ( Abs(tmp) <= tiny*Abs(out[i]) )
				mes_error(E_SING,"UTsolve");
			out[i] /= tmp;
			for(j=i+1;j<=dim;j++){
				out[j] += (conj(U[i][j])*(-out[i]));
			}
		}
    }
    else{
		invdiag = 1.0/diag;
		for (    ; i<=dim; i++ ){
			out[i] *= invdiag;
			for(j=i+1;j<=dim;j++){
				out[j] += (conj(U[i][j])*(-out[i]));
			}
		}
    }
    return 1;
}


/* Dsolve -- solves Dx=b where D is the diagonal of A -- may be in-situ */
template <typename T> int Matrix<T>::mes_Dsolve(const Matrix<T>& A,const Vector<T>& b,Vector<T>& x)
{
    int	dim, i;
    double	tiny;
    
    dim = local_min(NumRows(A),NumCols(A));
    if ( b.Size() < dim )
	mes_error(E_SIZES,"Dsolve");
    x.Resize(NumCols(A));

    tiny = 10.0/HUGE_VAL;

    dim = b.Size();
    for ( i=1; i<=dim; i++ )
	if ( Abs(A[i][i]) <= tiny*Abs(b[i]) )
	    mes_error(E_SING,"Dsolve");
	else
	    x[i] = b[i]/A[i][i];
    
    return 1;
}


/* LTsolve -- back substitution with optional over-riding diagonal
		using the LOWER triangular part of matrix
		-- can be in-situ but doesn't need to be */
template <typename T> int Matrix<T>::mes_LTsolve(const Matrix<T>& L,const Vector<T>& b,Vector<T>& out,double diag)
{
    int	dim;
    int		i, i_lim, j;
    double	invdiag, tiny;
	T tmp;
    
    dim = local_min(NumRows(L),NumCols(L));
    if ( b.Size() < dim )
	mes_error(E_SIZES,"LTsolve");
    out.Resize(NumCols(L));

    tiny = 10.0/HUGE_VAL;
    
    for ( i=dim; i>=1; i-- ){
		if ( !mtl_iszero(b[i]) )
			break;
	}
	i_lim = i;
	
	if ( &b[1] != &out[1] ){
		out=T();
		for(j=1;j<=i_lim;j++)
			out[j]=b[j];
	}

	if ( diag == 0.0 ){
		for ( ; i>=1; i-- ){
			tmp = conj(L[i][i]);
			if ( Abs(tmp) <= tiny*Abs(out[i]) )
				mes_error(E_SING,"LTsolve");
			out[i] /= tmp;
			for(j=1;j<i;j++){
				out[j] += conj(L[i][j])*(-out[i]);
			}
		}
	}
	else{
		invdiag = 1.0/diag;
		for ( ; i>=1; i-- ){
			out[i] *= invdiag;
			for(j=1;j<i;j++){
				out[j] += conj(L[i][j])*(-out[i]);
			}
		}
	}
	
	return 1;
}


#endif
