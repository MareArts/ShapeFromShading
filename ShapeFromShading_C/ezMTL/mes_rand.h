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
// Filename: mes_rand.h
// Revison:
//    1. Revised to support a singleton random number generator.
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_RAND_H_
#define	_MES_RAND_H_

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


///////////////////////////////////////////////////////////////////////////
// Pseudo random number generator data structures
// Knuth's lagged Fibonacci-based generator: See "Seminumerical Algorithms:
// The Art of Computer Programming" sections 3.2-3.3
///////////////////////////////////////////////////////////////////////////

struct Meschach_rand_para {
	 int  meschach_started;
	 int  meschach_inext;
	 int  meschach_inextp;
	 long meschach_mrand_list[56];
};


inline Meschach_rand_para* mes_get_meschach_rand_para()
{
	static Meschach_rand_para meschach_rand_para = { 0, 0, 31, 0};
	return &meschach_rand_para;
}


void mes_smrand(int seed); /* seeds mrand() */
template <typename T> void mes_mrandlist(T *x, int len); /* generates len random numbers */


/* mrand -- pseudo-random number generator */
template <typename T> inline double Matrix<T>::mes_mrand()
{
    long	lval;
    static double  factor = 1.0/((double)LONG_MAX);  
	Meschach_rand_para* p=mes_get_meschach_rand_para();

    if ( ! p->meschach_started )
		mes_smrand(3127);
    
    p->meschach_inext = (p->meschach_inext >= 54) ? 0 : p->meschach_inext+1;
    p->meschach_inextp = (p->meschach_inextp >= 54) ? 0 : p->meschach_inextp+1;
	
    lval = p->meschach_mrand_list[p->meschach_inext]-p->meschach_mrand_list[p->meschach_inextp];
    if ( lval < 0L )
		lval += LONG_MAX;
    p->meschach_mrand_list[p->meschach_inext] = lval;
    return (double)lval*factor;
}


/* mrandlist -- fills the array a[] with len random numbers */
template <typename T> inline int	Matrix<T>::mes_mrandlist(T *a, int len)
{
    int		i;
    long	lval;
    static double  factor = 1.0/((double)LONG_MAX);
	Meschach_rand_para* p=mes_get_meschach_rand_para();

    if ( ! p->meschach_started )
		mes_smrand(3127);
    
    for ( i = 0; i < len; i++ ){
		p->meschach_inext = (p->meschach_inext >= 54) ? 0 : p->meschach_inext+1;
		p->meschach_inextp = (p->meschach_inextp >= 54) ? 0 : p->meschach_inextp+1;
		
		lval = p->meschach_mrand_list[p->meschach_inext]-p->meschach_mrand_list[p->meschach_inextp];
		if ( lval < 0L )
			lval += LONG_MAX;
		p->meschach_mrand_list[p->meschach_inext] = lval;
		
		a[i] = T(lval*factor);
    }
	return 1;
}


/* smrand -- set seed for mrand() */
template <typename T> inline void Matrix<T>::mes_smrand(int seed)
{
    int		i;
	Meschach_rand_para* p=mes_get_meschach_rand_para();
	
    p->meschach_mrand_list[0] = (123413*seed) % LONG_MAX;
    for ( i = 1; i < 55; i++ )
		p->meschach_mrand_list[i] = (123413*p->meschach_mrand_list[i-1]) % LONG_MAX;
	
    p->meschach_started = 1;
	
    /* run mrand() through the list sufficient times to
	thoroughly randomise the array */
    for ( i = 0; i < 55*55; i++ )
		mes_mrand();
}


/* v_rand -- initialises x to be a random vector, components
independently & uniformly ditributed between 0 and 1 */
template <typename T> inline void	Matrix<T>::mes_v_rand(Vector<T>& x)
{
	mes_mrandlist(&x[1],x.Size());
}


#endif


