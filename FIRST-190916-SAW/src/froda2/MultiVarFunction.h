#ifndef MULTIVARFUNCTION_
#define MULTIVARFUNCTION_

#include <vector>

class MultiVarFunction
{
public:
  MultiVarFunction() {}
  virtual ~MultiVarFunction() {}

  virtual void getQ( std::vector<double> &q ) = 0;
  virtual void setQ( const std::vector<double> &q ) = 0;
  virtual double eval() = 0;
  virtual void getNegGrad( std::vector<double> &r ) = 0;
  virtual void getPreconditionedNegGrad( std::vector<double> &s, bool &isPreconditioned )
  {
    getNegGrad( s );
    isPreconditioned = false; 
  }
};

#endif /*MULTIVARFUNCTION_*/
