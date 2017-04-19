#ifndef MULTIVARFUNCTION_1DPROJECTION_H_
#define MULTIVARFUNCTION_1DPROJECTION_H_

class MultiVarFunction;
#include <vector>
#include <cstddef>

class MultiVarFunction_1DProjection {
public:
  MultiVarFunction_1DProjection();
  ~MultiVarFunction_1DProjection();
  void setProjection( MultiVarFunction *f_, const std::vector<double> &p_ );
  double operator()( double lambda_ );
  void update( double lambda_ );
    
private:
  MultiVarFunction *f;
  double lambda;
  size_t N;
  std::vector<double> q;
  std::vector<double> p;
  std::vector<double> newq;
};

#endif /*MULTIVARFUNCTION_1DPROJECTION_H_*/
