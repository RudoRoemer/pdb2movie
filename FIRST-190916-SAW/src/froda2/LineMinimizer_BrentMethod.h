#ifndef LINEMINIMIZER_BRENTMETHOD_H_
#define LINEMINIMIZER_BRENTMETHOD_H_

#include <cmath>
using namespace std;

template <typename T>
class LineMinimizer_BrentMethod
{
public:
	LineMinimizer_BrentMethod();
	virtual ~LineMinimizer_BrentMethod() {}
  double operator()( 
      T &f, int maxIterations, double fractionalTol, 
      double a, double x, double b,
      double fa, double fx, double fb );
  bool foundMinimum;
  int iteration;
  int maxIterations;
  double fractionalTol;
  double x;
  double fx;
  
private:
  const double tinyAbsoluteTol;
  const double golden;
};

template <typename T>
LineMinimizer_BrentMethod<T>::LineMinimizer_BrentMethod() :
  foundMinimum( false ),
  iteration( 0 ),
  maxIterations( 0 ),
  fractionalTol( 1.0e-3 ),
  tinyAbsoluteTol(1e-10),
  golden(0.5*(3.0 - sqrt(5.0))) // approx 0.382
  {}

template <typename T>
double LineMinimizer_BrentMethod<T>::operator()( 
      T &f, int maxIterations_, double fractionalTol_, 
      double a, double x_, double b,
      double fa, double fx_, double fb )
{
  // initialize the output data
  foundMinimum = false;
  iteration = 0;
  maxIterations = maxIterations_;
  fractionalTol = fractionalTol_;
  x = x_;
  fx = fx_;
  
  double v;
  double fv;
  double w;
  double fw;
  
  // verify that fx is less than both fa and fb,
  // and that a < x < b.
  // (The triple a,x,b brackets the minimum)
  if ( fx > fa || fx > fb || a >= x || x >= b ) {
    cout << "Error in LineMinimizer_BrentMethod: input values do not form a valid bracket " << endl;
    exit( 0 );
  }
  
  // initialize values of v and w
  // w is the second best minimum so far,
  // v is the third best so far
  // ( x is the best so far )
  if ( fa <= fb ) {
    w = a; fw = fa;
    v = b; fv = fb;
  }
  else {
    w = b; fw = fb;
    v = a; fv = fa;
  }
  
  // initialize d and d_oneIterationAgo
  // to big enough values that they won't stop the
  // first step from being a parabolic interpolation.
  double d = b - a;
  double d_oneIterationAgo = d * 2.0;
  
  for ( iteration = 0; iteration < maxIterations; iteration++ ) {
    double bracketMidpoint = 0.5*( b + a );
    double absoluteTol = fractionalTol*abs(x) + tinyAbsoluteTol;
    double absoluteTol_Times2 = absoluteTol * 2.0;
    double uncertaintyOfEstimate = abs(x-bracketMidpoint) + 0.5*( b - a );
    
    if ( uncertaintyOfEstimate <= absoluteTol_Times2 ) {
      //done.  Minimum is located to within tolerance.
      //Minimum is at x.
      foundMinimum = true;
      return x;
    }
    
    // update the stored previous values of d
    double d_twoIterationsAgo = d_oneIterationAgo;
    d_oneIterationAgo = d;
    
    // generate d
    double xw = x-w;
    double xv = x-v;
    double fxw = fx - fw;
    double fxv = fx - fv;
    
    double numerator;
    double denominator;
    if ( abs(d_twoIterationsAgo) > absoluteTol ) {
      numerator = xw*xw*fxv - xv*xv*fxw;
      denominator = 2.0*( xv*fxw - xw*fxv );
      d = numerator/denominator;
      
      // it is possible that the denominator was zero, and that a Nan
      // or Inf was generated.  If so, d is rejected.
      if ( abs(d) >= abs( 0.5*d_twoIterationsAgo ) || x+d <= a || x+d >= b || isnan(d) || isinf(d) ) {
        //reject parabolic guess
        d_oneIterationAgo = (x >= bracketMidpoint) ? a-x : b-x;
        d = golden * d_oneIterationAgo;
      }
      else {
        //accept parabolic guess
        
        //but first, check to see if d will take us
        //too close to one of the endpoints.
        //If so, it is likely that our trial guess is converging
        //onto the minimum, in the vicinity of the endpoint.
        //And, since f(x) is less than the values at a and b,
        //x is probably already close to that endpoint.
        if ( x+d-a < absoluteTol_Times2 || 
              b-x+d < absoluteTol_Times2 ) {
          d = ( x >= bracketMidpoint ) ? -absoluteTol : absoluteTol; 
        }
      }
    }
    else {
      d_oneIterationAgo = (x >= bracketMidpoint) ? a-x : b-x;
      d = golden * d_oneIterationAgo;
    }
    if ( abs(d) < absoluteTol ) d = ( d < 0.0 ) ? -absoluteTol : absoluteTol;
    
    //generate trial point u
    double u = x + d;
    double fu = f(u);
    
    if ( fu <= fx ) {
      // u is now the best minimum so far,
      // better than (or equal to) x.
      // First, we can tighten the bracket endpoints:
      //   if u>=x, (a,b) <- (x,b).  Otherwise (a,b) <- (a,x)
      if ( u >= x ) a = x;
      else b = x;
      
      //Second, we must update v,w,x
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    }
    else {
      // u is not the best minimum so far,
      // but it might be the second best
      // And certainly we can tighten the bracket endpoints:
      //   if u<x, (a,b) <- (u,b).  Otherwise, (a,b) <- (a,u)
      if ( u < x ) a = u;
      else b = u;
      
      // w keeps track of the second best minimum.
      // If u is the new second best so far, then w <- u.
      // ( Note below that if w==x, then w is
      //   not actually the second best but a copy of the best (x),
      //   which means that u would be the new second best. )
      if ( fu <= fw || w==x ) {
        v = w; fv = fw;
        w = u; fw = fu;
      }
      // v keeps track of the third best minimum ( not as good as x or w ).
      // if u is the new third best so far, then v <- u.
      else if ( fu <= fv || v==x || v==w ) {
        v = u; fv = fu;
      }
    }
    // return to the beginning of the for loop  
  }
  return x;
}

#endif /*LINEMINIMIZER_BRENTMETHOD_H_*/
