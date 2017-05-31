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
// Filename: mes_perm.h
// Revision:
//    1. 
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_PERM_H_
#define	_MES_PERM_H_	

/* Permutation matrix */
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


/**********************************************************************
Note: A permutation is often interpreted as a matrix
		(i.e. a permutation matrix).
	A permutation px represents a permutation matrix P where
		P[i][j] == 1 if and only if px[i] == j
**********************************************************************/


/* px_inv -- invert permutation -- in situ
	-- taken from ACM Collected Algorithms #250 */
template <typename T> int	Matrix<T>::mes_px_inv(Vector<int>& px,Vector<int>& out)
{
    int	i, j, k, n;
    
    out = px;
    for ( n=out.Size(); n>=1; n-- ){
		i = out[n];
		if ( i < 1 )
			out[n] = -1 - i;
		else if ( i != n ){
			k = n;
			while (true){
				if ( i < 1 || i > out.Size() )
					mes_error(E_BOUNDS,"px_inv");
				j = out[i];	out[i] = -1 - k;
				if ( j == n ){	
					out[n] = i;	break;		
				}
				k = i;		i = j;
			}
		}
    }
    return 1;
}


/* px_mlt -- permutation multiplication (composition) */
template <typename T> int	Matrix<T>::mes_px_mlt(Vector<int>& px1,Vector<int>& px2,Vector<int>& out)
{
    int	i,size;
    
    if ( px1.Size() != px2.Size() )
	mes_error(E_SIZES,"px_mlt");
    if ( &px1[1] == &out[1] || &px2[1] == &out[1] )
	mes_error(E_INSITU,"px_mlt");
    if ( out.Size() != px1.Size() )
	out.Resize(px1.Size());
    
    size = px1.Size();
    for ( i=1; i<=size; i++ )
	if ( px2[i] > size )
	    mes_error(E_BOUNDS,"px_mlt");
	else
	    out[i] = px1[px2[i]];
    
    return 1;
}


/* px_vec -- permute vector */
template <typename T> int	Matrix<T>::mes_px_vec(Vector<int>& px,Vector<T>& vec,Vector<T>& out)
{
    int	old_i, i, size, start;
    T	tmp;
    
    if ( px.Size() > vec.Size() )
	mes_error(E_SIZES,"px_vec");
    if ( out.Size() != vec.Size() )
	out.Resize(vec.Size());
    
    size = px.Size();
    if ( size == 0 ){
		out=vec;
		return 1;
	}
	if(&out[1] != &vec[1]){
		for ( i=1; i<=size; i++ ){
			if ( px[i] > size )
				mes_error(E_BOUNDS,"px_vec");
			else
				out[i] = vec[px[i]];
		}
    }
    else {	/* in situ algorithm */
		start = 1;
		while ( start <= size ){
			old_i = start;
			i = px[old_i];
			if ( i > size ){
				start++;
				continue;
			}
			tmp = vec[start];
			while ( true ){
				vec[old_i] = vec[i];
				px[old_i] = i+size;
				old_i = i;
				i = px[old_i];
				if ( i > size )
					break;
				if ( i == start ){
					vec[old_i] = tmp;
					px[old_i] = i+size;
					break;
				}
			}
			start++;
		}
		
		for ( i = 1; i <= size; i++ )
			if ( px[i] <= size )
				mes_error(E_BOUNDS,"px_vec");
			else
				px[i] = px[i]-size;
    }
    
    return 1;
}


/* pxinv_vec -- apply the inverse of px to x, returning the result in out */
template <typename T> int	Matrix<T>::mes_pxinv_vec(Vector<int>& px,Vector<T>& x,Vector<T>& out)
{
    int	i, size;
    
    if ( px.Size() > x.Size() )
		mes_error(E_SIZES,"pxinv_vec");
    if ( out.Size() != x.Size() )
		out.Resize(x.Size());
    
    size = px.Size();
    if ( size == 0 ){
		out=x;
		return 1;
	}
    if ( &out[1] != &x[1] ){
		for ( i=1; i<=size; i++ )
			if ( px[i] > size )
				mes_error(E_BOUNDS,"pxinv_vec");
			else
				out[px[i]] = x[i];
    }
    else{	/* in situ algorithm --- cheat's way out */
		mes_px_inv(px,px);
		mes_px_vec(px,x,out);
		mes_px_inv(px,px);
    }
	
    return 1;
}


/* px_transp -- transpose elements of permutation
             -- Really multiplying a permutation by a transposition
   ps: permutation to transpose
   i1, i2: elements to transpose */
template <typename T> inline int	Matrix<T>::mes_px_transp(Vector<int>& px,int i1,int i2)
{
	int	temp;

	if ( i1 <= px.Size() && i2 <= px.Size() ) {
		temp = px[i1];
		px[i1] = px[i2];
		px[i2] = temp;
	}

	return 1;
}


/* myqsort -- a cheap implementation of Quicksort on integers
		-- returns number of swaps */
template <typename T> int	Matrix<T>::mes_myqsort(int* a,int num)
{
	int	i, j, tmp, v;
	int	numswaps;

	numswaps = 0;
	if ( num <= 1 )
		return 0;

	i = 0;	j = num;	v = a[0];
	while(true) {
		/////////////////////////////////////////////
#if 0
		while ( a[++i] < v )
			;
		while ( a[--j] > v )
			;
#endif
		// new code
		while(1){
			++i;
			if(i>=num || a[i]>=v) break;
		}
		while(1){
			--j;
			if(j<0 || a[j]<=v) break;
		}
		/////////////////////////////////////////////

		if ( i >= j )	break;

		tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
		numswaps++;
	}

	tmp = a[0];
	a[0] = a[j];
	a[j] = tmp;
	if ( j != 0 )
		numswaps++;

	numswaps += mes_myqsort(&a[0],j);
	numswaps += mes_myqsort(&a[j+1],num-(j+1));

	return numswaps;
}


/* px_sign -- compute the ``sign'' of a permutation = +/-1 where
		px is the product of an even/odd # transpositions */
template <typename T> int	Matrix<T>::mes_px_sign(Vector<int>& px)
{
	int	numtransp;
	Vector<int>	px2;

	px2 = px;
	numtransp = mes_myqsort(&px2[1],px2.Size());
	return ( numtransp % 2 ) ? -1 : 1;
}


/* px_cols -- permute columns of matrix A; out = A.px'
	-- May NOT be in situ */
template <typename T> int	Matrix<T>::mes_px_cols(Vector<int>& px,Matrix<T>& A,Matrix<T>& out)
{
	int	i, j, m, n, px_j;

	if ( px.Size() != NumCols(A) )
		mes_error(E_SIZES,"px_cols");
	if ( &A[1] == &out[1] )
		mes_error(E_INSITU,"px_cols");
	m = NumRows(A);	n = NumCols(A);
	if ( NumRows(out) != m || NumCols(out) != n )
		out.Resize(m,n);

	for ( j = 1; j <= n; j++ ){
		px_j = px[j];
		if ( px_j > n )
		    mes_error(E_BOUNDS,"px_cols");
		for ( i = 1; i <= m; i++ )
		    out[i][px_j] = A[i][j];
	}

	return 1;
}


/* px_rows -- permute columns of matrix A; out = px.A
	-- May NOT be in situ */
template <typename T> int	Matrix<T>::mes_px_rows(Vector<int>& px,Matrix<T>& A,Matrix<T>& out)
{
	int	i, j, m, n, px_i;

	if ( px.Size() != NumRows(A) )
		mes_error(E_SIZES,"px_rows");
	if ( &A[1] == &out[1] )
		mes_error(E_INSITU,"px_rows");
	m = NumRows(A);	n = NumCols(A);
	if ( NumRows(out) != m || NumCols(out) != n )
		out = m_get(m,n);

	for ( i = 1; i <= m; i++ ){
		px_i = px[i];
		if ( px_i > m )
		    mes_error(E_BOUNDS,"px_rows");
		for ( j = 1; j <= n; j++ )
		    out[i][j] = A[px_i][j];
	}

	return 1;
}


#endif
