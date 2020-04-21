/*
 * BrentsRootFinder.hpp
 *
 *  Created on: Sep 27, 2013
 *      Author: marcel
 */

#ifndef BRENTSROOTFINDER_HPP_
#define BRENTSROOTFINDER_HPP_

#include <cmath>

namespace epital {

namespace BrentsRootFinder {

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}

inline long double SIGN(const long double &a, const long double &b)
	{return b >= 0.0l ? (a >= 0.0l ? a : -a) : (a >= 0.0l ? -a : a);}

template <class F>
long double zbrent(F func, const long double x1, const long double x2, const long double tol)
{
	const int ITMAX=100;
	const long double EPS=numeric_limits<long double>::epsilon();
	long double a=x1,b=x2,c=x2,d,e,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	if ((fa > 0.0l && fb > 0.0l) || (fa < 0.0l && fb < 0.0l))
		throw("Root must be bracketed in zbrent");
	fc=fb;
	for (int iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0l && fc > 0.0l) || (fb < 0.0l && fc < 0.0l)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (std::abs(fc) < std::abs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0l*EPS*std::abs(b)+0.5l*tol;
		xm=0.5l*(c-b);
		if (std::abs(xm) <= tol1 || fb == 0.0l) return b;
		if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0l*xm*s;
				q=1.0l-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0l*xm*q*(q-r)-(b-a)*(r-1.0l));
				q=(q-1.0l)*(r-1.0l)*(s-1.0l);
			}
			if (p > 0.0l) q = -q;
			p=std::abs(p);
			long double min1=3.0l*xm*q-std::abs(tol1*q);
			long double min2=std::abs(e*q);
			if (2.0l*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (std::abs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	throw("Maximum number of iterations exceeded in zbrent");
}


template <class F, class T>
T zbrent(F func, const T x1, const T x2, const T tol)
{
	const int ITMAX=100;
	const T EPS=numeric_limits<T>::epsilon();
	T a=x1,b=x2,c=x2,d,e,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		throw("Root must be bracketed in zbrent");
	fc=fb;
	for (int iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (std::abs(fc) < std::abs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*std::abs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (std::abs(xm) <= tol1 || fb == 0.0) return b;
		if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=std::abs(p);
			T min1=3.0*xm*q-std::abs(tol1*q);
			T min2=std::abs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (std::abs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	throw("Maximum number of iterations exceeded in zbrent");
}


template <class F, class T>
T zbrent_efloat(F func, const T x1, const T x2, const T tol)
{
	const int ITMAX=100;
	const T EPS=numeric_limits<T>::epsilon();
	T a=x1,b=x2,c=x2,d,e,fa=func(a),fb=func(b),fc,p,q,r,s,tol1,xm;
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		throw("Root must be bracketed in zbrent");
	fc=fb;
	for (int iter=0;iter<ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (std::abs(fc) < std::abs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*std::abs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (std::abs(xm) <= tol1 || fb == 0.0) return b;
		if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=std::abs(p);
			T min1=3.0*xm*q-std::abs(tol1*q);
			T min2=std::abs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (std::abs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
			fb=func(b);
	}
	throw("Maximum number of iterations exceeded in zbrent");
}


}


} /* namespace epital */
#endif /* BRENTSROOTFINDER_HPP_ */
