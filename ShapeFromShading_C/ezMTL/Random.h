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
// Filename: Random.h
//
// Note:
//  1. Distributions other than uniform and Gaussian were not tested.
//     You need to test the resulting distribution before use.
//     In particular, the Gamma distribution has a poor pdf when a<=1.
//
///////////////////////////////////////////////////////////////////////////////
#ifndef _RANDOM_H_
#define _RANDOM_H_

// Undefine rand and srand to avoid name collisions
// Use RANDF and SRAND instead.
#ifdef USE_LOCAL_RAND
#ifdef rand
#undef rand
#endif
#ifdef srand
#undef srand
#endif
#else
#include <stdlib.h>
#endif


// RANDF() returns a uniform random number [0,1).
#if defined(USE_LOCAL_RAND)
#define RANDF() UniformRandom()
#define SRAND(x) SRand(x)
#define SRANDN(x) SRandn(x)
double Rand();
void SRand(unsigned int seed);
#elif defined(__unix)
// Prototype for C Library functions deand48 and srand48.
double drand48();
void srand48(long);
#define RANDF() drand48()
#define SRAND(x) srand48(x)
#define SRANDN(x) SRandn(x)
#elif defined(WIN32)
// If not unix use ANSI C defaults.
#define __rand() ((double)rand()/RAND_MAX)
#define RANDF() __rand()
#define SRAND(x) srand(x)
#define SRANDN(x) SRandn(x)
#else
#pragma message ( "Your system does not have random number generators " __FILE__ )
#pragma message ( "Define __unix or WIN32 " )
!!error!! // to stop compiling
#endif


double UniformRandom();
double GaussianRandom();
void SRandn(unsigned int seed);


template <typename T> class UniformDist;
template <typename T> class GaussianDist;
template <typename T> class MiscellaneousDist;

///////////////////////////////////////////////////////////////////////////////
// Exported routines
///////////////////////////////////////////////////////////////////////////////

// continuous distributions
#ifdef USE_LOCAL_RAND
inline double UniformRandom();
inline void SRand(unsigned int seed);
#endif
inline double GaussianRandom();
inline void SRandn(unsigned int seed);
inline double ExponentialRandom(double lambda);
inline double TRandom(double n); // Student's t-distribution with n degree of freedom
inline double GammaRandom(double a, double b);
inline double BetaRandom(double a, double b);
inline double FisherRandom(double m, double n); // m degrees of freedom in the numerator and n degrees of freedom in the denominator
inline double ChiSquareRandom(double n); // n degree of freedom
inline double CauchyRandom();

// discrete distributions
inline int BernoulliRandom(double p);
inline int BinomialRandom(int n,double p);
inline int PossionRandom(double lambda);


///////////////////////////////////////////////////////////////////////////////
// Implementation of random number generation
///////////////////////////////////////////////////////////////////////////////

#ifdef USE_LOCAL_RAND
struct Int32__ {
	unsigned short high;
	unsigned short low;
};


///////////////////////////////////////////////////////////////////////////////
// There must exist only one uniform random number generator in the whole system.
// I had to implement a singleton object using a static member function and 
// a static object defined within the static member function.
///////////////////////////////////////////////////////////////////////////////
template <typename T>
class UniformDist {
public:
	UniformDist() {}
	virtual ~UniformDist() {}
	static UniformDist<double>* GetUniformDist();
	static void SRand_(unsigned int seed);
	static double Rand_();
private:
	unsigned int rand_seed;
	Int32__ rand_value;
};


template <typename T> UniformDist<double>* UniformDist<T>::GetUniformDist(){
	static UniformDist<double> uniformDist;
	return &uniformDist;
}


template <typename T> void UniformDist<T>::SRand_(unsigned int seed){
	UniformDist* p=UniformDist::GetUniformDist();
	p->rand_seed=seed;
	p->rand_value.high=seed/65536;
	p->rand_value.low=seed%65536;
}


///////////////////////////////////////////////////////////////////////////////
// Uniform distribution [0,1)
//   x = (69069*x+1)%pow(2,32);
//   u = x/pow(2,32)
// Reference: C.P. Robert, The Bayesian Choice, 2nd ed., Springer, 2001.
///////////////////////////////////////////////////////////////////////////////
template <typename T> double UniformDist<T>::Rand_()
{
	UniformDist* p=UniformDist::GetUniformDist();
	static Int32__ c={1,3533}; // 69069
	static double rand_max=pow((double)2.0,(double)32);

	unsigned int x1=c.low*p->rand_value.low+1;
	unsigned int x2=c.low*p->rand_value.high;
	unsigned int x3=c.high*p->rand_value.low;
	unsigned int x4=c.high*p->rand_value.high;
	p->rand_value.low=x1 & 0xffff;
	unsigned int x1h=(unsigned int)(0xffff & ((x1&0xffff0000)>>16));
	unsigned int x2l=(unsigned int)(x2&0xffff);
	unsigned int x3l=(unsigned int)(x3&0xffff);
	p->rand_value.high=(unsigned int)(x1h+x2l+x3l);
	double u=(unsigned int)(65536*p->rand_value.high+p->rand_value.low)/rand_max;
	return u;
}

#endif


///////////////////////////////////////////////////////////////////////////////
// Gaussian distribution
///////////////////////////////////////////////////////////////////////////////
template <typename T>
class GaussianDist {
public:
	GaussianDist() { gaussSaved=0; }
	virtual ~GaussianDist() {}
	static GaussianDist<double>* GetGaussianDist();
	static void SRandn_(unsigned int seed);
	static double Randn_();
private:
	int gaussSaved; // GaussDeviate generates numbers in pairs
    double gaussSave; // 2nd of pair is remembered here
#ifdef USE_LOCAL_RAND
	unsigned int rand_seed;
	Int32__ rand_value;
#endif

};


template <typename T> GaussianDist<double>* GaussianDist<T>::GetGaussianDist(){
	static GaussianDist<double> gaussianDist;
	return &gaussianDist;
}


template <typename T> void GaussianDist<T>::SRandn_(unsigned int seed){
	GaussianDist* p=GaussianDist<double>::GetGaussianDist();
#ifdef USE_LOCAL_RAND
	p->rand_seed=seed;
	p->rand_value.high=seed/65536;
	p->rand_value.low=seed%65536;
#endif
	p->gaussSaved=0;
	p->gaussSave=0;
}


// random number with a N(0,1) distribution
template <typename T> double GaussianDist<T>::Randn_()
{
	GaussianDist* p=GaussianDist<double>::GetGaussianDist();
   double fac,r,v1,v2,x;
   if (p->gaussSaved) {
      x = p->gaussSave; p->gaussSaved = 0;
   }
   else {
      do {
#ifdef USE_LOCAL_RAND
         v1 = 2.0*UniformRandom() - 1.0;
         v2 = 2.0*UniformRandom() - 1.0;
#else
		 v1 = 2.0*RANDF() - 1.0;
         v2 = 2.0*RANDF() - 1.0;
#endif
         r = v1*v1 + v2*v2;
      }
      while (r>=1.0);
      fac = sqrt(-2.0*log(r)/r);
      p->gaussSaved = 1;
      p->gaussSave = v1*fac;
      x = v2*fac;
   }
   return x;
}


///////////////////////////////////////////////////////////////////////////////
// Other distributions
///////////////////////////////////////////////////////////////////////////////
template <typename T> 
class MiscellaneousDist {
public:
	// continuous distributions
	static double ExponentialRandom_(double lambda=1);
	static double TRandom_(double n);
	static double GammaRandom_(double a, double b=1);
	static double BetaRandom_(double a, double b);
	static double FisherRandom_(double m, double n);
	static double ChiSquareRandom_(double lambda);
	static double CauchyRandom_();

	// discrete distributions
	static int BernoulliRandom_(double p);
	static int BinomialRandom_(int n,double p);
	static int NegativeBinomialRandom_(int n,double p);
	static int PossionRandom_(double lambda);
	static int HyperGeometricalbRandom_(int N,int n,double p);

};


template <typename T> double MiscellaneousDist<T>::ExponentialRandom_(double lambda){
	double u;
	while(1){
		u=UniformRandom();
		if(u!=0) break;
	}
	return -log(u)/lambda;
}


// n degree of freedom
template <typename T> double MiscellaneousDist<T>::TRandom_(double n){
	while(1){
		double u1=UniformRandom();
		double u2=UniformRandom();
		double x,v;
		if(u1 < 0.5){
			x=1/(4*u1-1);
			v=u2/(x*x);
		}
		else{
			x=4*u1-3;
			v=u2;
		}
		if( v < 1-(fabs(x)/2) || v < pow(1+(x*x/n), -(n+1)/2.0) )
			return x;
	}
}


template <typename T> double MiscellaneousDist<T>::GammaRandom_(double a, double b){
	assert(a>0 && b>0);
	double u1,u2;
	if(a > 50){
		return GaussianRandom()*sqrt(a)/b + (a/b);
	}
	else if(a > 1){
		double c1=a-1, c2=(a-(1/(6*a)))/c1, c3=2/c1, c4=1+c3, c5=1/sqrt(a);
		while(1){
			while(1){
				u1=UniformRandom();
				u2=UniformRandom();
				if(a > 2.5){
					u1 = u2 + c5 * (1-1.86*u1);
				}
				if(u1 > 0 && u1 < 1) break;
			}
			double w=c2*u2/u1;
			if(c3*u1+w+1/w <= c4 || c3*log(u1)-log(w)+w <= 1)
				return (b*c1*w);
		}
	}
	else{
		double e=exp(1);
		while(1){
			u1=UniformRandom();
			u2=UniformRandom();
			double x,y;
			if(u1 > e/(e+a)){
				x=-log((a+e)*(1-u1)/(a*e));
				y=pow(x,a-1);
			}
			else{
				x=pow((a+e)*u1*e,1/a);
				y=exp(-a);
			}
			if(u2<y)
				return (b*x);
		}
	}
}


template <typename T> double MiscellaneousDist<T>::BetaRandom_(double a, double b){
	double y1=GammaRandom(a,1);
	double y2=GammaRandom(b,1);
	return y1/(y1+y2);
}


// m degrees of freedom in the numerator and n degrees of freedom in the denominator
template <typename T> double MiscellaneousDist<T>::FisherRandom_(double m, double n){
	double xnom=ChiSquareRandom(m)/m;
	double xdenom;
	while(1){
		xdenom=ChiSquareRandom(n)/n;
		if(xdenom > xnom/mtl_numeric_limits<double>::max()) break;
	}
	return xnom/xdenom;
}


// n degree of freedom
template <typename T> double MiscellaneousDist<T>::ChiSquareRandom_(double n){
	return GammaRandom(n/2.0,2.0);
}


template <typename T> double MiscellaneousDist<T>::CauchyRandom_(){
	static double pi=4*atan(1);
	//f(x)=1/(pi*(1+x^2))
	return tan(pi*(UniformRandom()-0.5));
}


template <typename T> int MiscellaneousDist<T>::BernoulliRandom_(double p){
	if(UniformRandom()<=p) return 1;
	else return 0;
}


template <typename T> int MiscellaneousDist<T>::BinomialRandom_(int n,double p){
	int K=n;
	if(n <= 30){
		int count = 0;
		for(int i=1;i<=n;i++){
			double u=UniformRandom();
			if(u < p) count++;
		}
		return count;
	}
	else{
		int i, k=n, x=0;
		double l=p;
		while(1){
			i=(int)(1+k*l);
			double v=BetaRandom(i,k+1-i);
			if(l>v){
				l=l/v;
				k=i-1;
			}
			else{
				x=x+i;
				l=(l-v)/(1-v);
				k=k-i;
			}
			if(k<=K) break;
		}
		for(i=1;i<=k;i++){
			double u=UniformRandom();
			if(u<p) x++;
		}
		return x;
	}
}


template <typename T> int MiscellaneousDist<T>::PossionRandom_(double lambda){
	if(lambda<30){
		double p=1;
		int N=0;
		double c=exp(-lambda);
		while(1){
			double u=UniformRandom();
			p=p*u;
			N=N+1;
			if(p<c) break;
		}
		return N-1;
	}
	else{
		static double pi=4*atan(1);
		double c=0.767-3.36/lambda;
		double beta=pi*pow((2*lambda),-0.5);
		double alpha=beta*lambda;
		double k=log(c)-lambda-log(beta);
		double x;
		while(1){
			while(1){
				double u1=UniformRandom();
				x=(alpha-log((1-u1)/u1))/beta;
				if(x>-0.5) break;
			}
			double u2=UniformRandom();
			int N=(int)(x+0.5);
			if( alpha-beta*x+log(u2/(1+exp(2*(alpha-beta*x)))) <= k+N*log(lambda)-log(N) )
				return N;
		}
	}
}


////////////////////////////////////////////////////////////////////////////
// Implementation of the exported routines
////////////////////////////////////////////////////////////////////////////

#ifdef USE_LOCAL_RAND
inline double UniformRandom(){
	return UniformDist<double>::Rand_();
}


inline void SRand(unsigned int seed){
	UniformDist<double>::SRand_(seed);
}
#endif


inline double GaussianRandom(){
	return GaussianDist<double>::Randn_();
}


inline void SRandn(unsigned int seed){
#ifdef USE_LOCAL_RAND
	GaussianDist<double>::SRandn_(seed);
#endif
}


inline double ExponentialRandom(double lambda){
	return MiscellaneousDist<double>::ExponentialRandom_(lambda);
}


// n degree of freedom
inline double TRandom(double n){
	return MiscellaneousDist<double>::TRandom_(n);
}


inline double GammaRandom(double a, double b){
	return MiscellaneousDist<double>::GammaRandom_(a,b);
}


inline double BetaRandom(double a, double b){
	return MiscellaneousDist<double>::BetaRandom_(a,b);
}


// m degrees of freedom in the numerator and n degrees of freedom in the denominator
inline double FisherRandom(double m, double n){
	return MiscellaneousDist<double>::FisherRandom_(m,n);
}


// n degree of freedom
inline double ChiSquareRandom(double n){
	return MiscellaneousDist<double>::ChiSquareRandom_(n);
}


inline double CauchyRandom(){
	return MiscellaneousDist<double>::CauchyRandom_();
}


inline int BernoulliRandom(double p){
	return MiscellaneousDist<double>::BernoulliRandom_(p);
}


inline int BinomialRandom(int n,double p){
	return MiscellaneousDist<double>::BinomialRandom_(n,p);
}


inline int PossionRandom(double lambda){
	return MiscellaneousDist<double>::PossionRandom_(lambda);
}


#endif
