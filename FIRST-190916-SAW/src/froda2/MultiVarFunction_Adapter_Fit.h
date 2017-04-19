#ifndef MULTIVARFUNCTION_ADAPTER_FIT_H_
#define MULTIVARFUNCTION_ADAPTER_FIT_H_

#include "MultiVarFunction.h"
class Fit;

class MultiVarFunction_Adapter_Fit : public MultiVarFunction
{
public:
	MultiVarFunction_Adapter_Fit( Fit *fit_ );
	virtual ~MultiVarFunction_Adapter_Fit();

  void getQ( std::vector<double> &q );
  void setQ( const std::vector<double> &q );
  double eval();
  void getNegGrad( std::vector<double> &r );
  void getPreconditionedNegGrad( std::vector<double> &s, bool &isPreconditioned );
  
private:
  Fit *fit;
};

#endif /*MULTIVARFUNCTION_ADAPTER_FIT_H_*/
