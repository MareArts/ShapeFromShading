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
// Filename: mes_err.h
// Revision:
//    1. Added error and warning messages to support complex matrices
///////////////////////////////////////////////////////////////////////////////

#ifndef	_MES_ERR_H_
#define	_MES_ERR_H_	

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

int mes_ev_error(char *,int,int,char *);  /* main error handler */
int mes_ev_warn(char *,int,int,char *);  /* main warning handler */

#define	mes_error(err_num,fn_name)	mes_ev_error(__FILE__,err_num,__LINE__,fn_name)
#define mes_warn(warn_num,fn_name) mes_ev_warn(__FILE__,warn_num,__LINE__,fn_name) 

/* error types */
enum mes_error_type {
	E_UNKNOWN,E_SIZES,E_BOUNDS,E_MEM,E_SING,E_POSDEF,E_FORMAT,E_INPUT,E_NULL,E_SQUARE,
		E_RANGE,E_INSITU2,E_INSITU,E_ITER,E_CONV,E_START,E_SIGNAL,E_INTERN,E_EOF,E_SHARED_VECS,
		E_NEG,E_OVERWRITE,E_BREAKDOWN,E_NEED_REAL,E_NEED_COMPLEX
};

/* warning types */
enum mes_warn_type {
	WARN_UNKNOWN,WARN_WRONG_TYPE,WARN_NO_MARK,WARN_RES_LESS_0,WARN_SHARED_VEC,
		WARN_SINGLE_PRECISION
};

static	char	*err_mesg[] = {
	  "unknown error",			    /* 0 */
	  "sizes of objects don't match",	    /* 1 */
	  "index out of bounds",		    /* 2 */
	  "can't allocate memory",		    /* 3 */
	  "singular matrix",			    /* 4 */
	  "matrix not positive definite",	    /* 5 */
	  "incorrect format input",		    /* 6 */
	  "bad input file/device",		    /* 7 */
	  "NULL objects passed",		    /* 8 */
	  "matrix not square",			    /* 9 */
	  "object out of range",		    /* 10 */
	  "can't do operation in situ for non-square matrix",   /* 11 */
	  "can't do operation in situ",		    /* 12 */
	  "excessive number of iterations",	    /* 13 */
	  "convergence criterion failed",	    /* 14 */
	  "bad starting value",			    /* 15 */
	  "floating exception",			    /* 16 */
	  "internal inconsistency (data structure)",/* 17 */
	  "unexpected end-of-file",		    /* 18 */
	  "shared vectors (cannot release them)",   /* 19 */  
	  "negative argument",			    /* 20 */
	  "cannot overwrite object",                /* 21 */
	  "breakdown in iterative method",          /* 22 */
	  "must be real matrix",           /* 23 */
	  "must be complex matrix"           /* 24 */
	 };

static char *warn_mesg[] = {
   "unknown warning",				  /* 0 */
   "wrong type number (use macro TYPE_*)",	  /* 1 */
   "no corresponding mem_stat_mark",		  /* 2 */
   "computed norm of a residual is less than 0",  /* 3 */
   "resizing a shared vector",			  /* 4 */
   "insufficient precision (use double precision)"		  /* 5 */
};

/* ev_error -- reports error (err_num) in file "file" at line "line_num" and exit */
inline int	mes_ev_error(char *file,mes_error_type err_num,int line_num,char* fn_name) {
	fprintf(stderr,"ERROR: %s\n",err_mesg[err_num]);
	assert(0);
	exit(1);
	return 0;
}


inline int	mes_ev_warn(char *file,mes_warn_type warn_num,int line_num,char* fn_name) {
	fprintf(stderr,"WARNING: %s\n",warn_mesg[warn_num]);
	assert(0);
	return 0;
}

#endif


