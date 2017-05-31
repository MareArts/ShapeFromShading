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
// Filename: Matrix.h
//
// Version: 1.02 May 24, 2002
//
// Programmer: Oh-Wook Kwon (ohwook@yahoo.com)
//
// Installed Operating Sytems/Compilers:
//     Windows 32/Visual C++ 6.0
//     Windows 32/Visual C++ 6.0 MFC
//     Red Hat Linux 7.1/GNU C++ 2.96
//     FreeBSD 4.1.1/GNU C++ 2.95.2
//
// Installation to Unix or Win32:
//   1. Unzip the distributed package in a temporary directory.
//      Two directories will be created under ./ezMTL-1.0: 'ezmtl' and 'test'.
//   2. Copy or move the './ezMTL-1.0/ezmtl' directory to the target directory
//      (e.g., ./ezmtl).
//
//   Notes for Visual C++ with MFC:
//       It is highly recommended that the ezMTL header file (Matrix.h) should be 
//       included before including <afxwin.h> in the 'StdAfx.h'.
// 
// Usage:
//   1. Include the library in your programs by using
//      #include "./ezmtl/Matrix.h"
//
//   2. Compile your own programs.
//
// Implementation:
//   1. Template library for matrix algebra.
//        - Most of elementary matrix operations (+, -, *, /, transpose, conjugate, ...)
//        - Some of basic statistical operations (sum, mean, var, std, skewness, kurtosis, ...)
//        - Checking matrix properties (symmetric, Hermitian, positive definite, diagonal, triangular, zero, identity, ...)
//        - Cholesky factorization for square real/complex symmetric positive definite matrices
//        - LU factorization for square real/complex matrices
//        - QR factorization for general real/complex matrices
//        - Singular value decompostion for general real matrices
//        - Hessenberg form for general real/complex matrices
//        - Schur decomposition for general real/complex matrices
//        - BKP decomposition for general real matrices
//        - Eigenvectors and Eigenvalues for symmetric/Hermitian matrices
//        - Eigenvectors and Eigenvalues for general real/complex matrices
//        - Inverse, pseudoinverse, determinant, p-norm
//        - Rank, orthogonal/null space
//        - Condition numbers
//        - Solving linear equations
//        - Matrix exponential, logarithm, square root, power
//        - Random number generation (uniform, normal, exponential, gamma, beta, t, fisher, chi-square, binomial, poisson, ...)
//
//   2. This library was designed to make MATLAB codes run in C++.
//      Most of method names are compatible to MATLAB commands
//      except the first character is capital letter
//
//   3. This library is partly based on the Meschach Library.
//
//   4. All matrix and vector indices start from 1 NOT 0.
//
//   5. This library privides matrix element access by both A(i,j) and A[i][j].
//      I recommend to access matrix elements by A(i,j) 
//      because A(i,j) loses nothing but may gain something at least in some cases.
//
// Programming tips:
//   1. Prefer A+=B to A=A+B because the latter saves the result of addition 
//      to a temporary matrix and then copy it to A.
//      Use
//        A+=B; A*=B; A-=B; A/=b;
//      rather than
//        A=A+B; A=A*B; A=A-B; A=A/b;
//
//   2. Avoid initialization in the loop. Get the initialization out of the loop.
//      The C++ compiler might call the copy constructor for each iteration.
//      Use
//        Matrix<float> A
//        for (i=1;i<=1000;i++){
//          A=Matrix<float>::Randn(10,10);
//        }
//      rather than
//        for (i=1;i<=1000;i++){
//          Matrix<float> A=Matrix<float>::Randn(10,10);
//        }
//
// Notes:
//    1. This library is adequate for quick-and-dirty simulation of matrix operations.
//       You have only to pay a little attention when you want to your MATLAB programs
//       faster and more memory effienct.
//
//    2. The performance of this library lies between the MATLAB and C programs
//       in terms of speed and simplicity (not checked exactly).
//       You had better consider other libraries if your applications
//       require very high memory efficiency and speed like computer games.
//
//    3. This library was designed to be used for dense matrices with medium size (< 1000x1000).
//       The memory efficiency shall be poor for sparse or very large matrices (e.g., 10000x10000).
//
// Things to do:
//   1. To support more complex linear algebra functions:
//        - Singular value decomposition for complex matrices
//        - Expm, logm, sqrtm for complex matrices
//           
//   2. To support sparse matrices.
//
//   3. ...
//
// Bug Report:
//   1. If you find bugs or supplement the above functions,
//      please send me an e-mail. ohwook@yahoo.com
//
// Revision History:
//
//   Version 1.02 May 21, 2002
//   MatrixConfig.h 
//   Commented out undefining max and min macros.
//
//   Version 1.01 March 19, 2002
//   Debugged a compile error on Linux.
//
//   Version 1.0 November 29, 2001
//  Released the public version.
//
///////////////////////////////////////////////////////////////////////////////
