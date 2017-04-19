#include "MultiVarFunction_1DProjection.h"
#include "MultiVarFunction.h"
#include <cmath>
#include <limits>

using namespace std;

MultiVarFunction_1DProjection::MultiVarFunction_1DProjection() : 
  f(NULL),
  lambda(0),
  N(0)
{}

MultiVarFunction_1DProjection::~MultiVarFunction_1DProjection() {}

void MultiVarFunction_1DProjection::setProjection( MultiVarFunction *f_, const std::vector<double> &p_ ) {
  f = f_;
  f->getQ( q );
  p = p_;
  N = q.size();
  newq.resize( N );
  lambda = 0.0;
}

void MultiVarFunction_1DProjection::update( double lambda_ ) {

  //cout << "Updating 1D model to lambda = " << lambda_ << endl;
  double deltaLambda = lambda_ - lambda;
  if ( abs(deltaLambda) < numeric_limits<double>::epsilon() ) {
    // The requested 1-D update to lambda will have no effect
    // on the underlying Rigid Model.
    //cout << "1D Update to lambda skipped\n";
    return;
  }
  
  lambda = lambda_;
  //cout << "1D Update to " << lambda << "\n";
  for ( size_t i=0; i<N; i++ ) {
    newq[i] = q[i] + lambda*p[i];
  }
  
  f->setQ( newq );
  
}

double MultiVarFunction_1DProjection::operator()( double lambda_ ) {
  update( lambda_ );
  return f->eval();
}
