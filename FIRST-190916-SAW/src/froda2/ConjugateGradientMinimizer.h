#ifndef CONJUGATEGRADIENTMINIMIZER_H_
#define CONJUGATEGRADIENTMINIMIZER_H_

class MultiVarFunction;
#include "MultiVarFunction_1DProjection.h"
#include "Mismatch.h"
#include "FindBracket.h"
#include "LineMinimizer_BrentMethod.h"
#include <vector>
#include "Observable.h"

class ConjugateGradientMinimizer : public Observable
{
public:
  ConjugateGradientMinimizer();
  virtual ~ConjugateGradientMinimizer();

  void minimize( MultiVarFunction *f_, int nsteps );

  bool isTolReached() const { return reachedTol; }
  int getFinalStepNum() const { return step; }
  int getNumCompletedSteps() const { return step; }

  //settings
  void setOutputPeriod( int outputPeriod_ ) {
    outputPeriod = outputPeriod_;
  } 
  void setMethod_SteepestDescent() { 
    doSteepestDescent_or_doConjGrad_0or1 = 0;
  }
  void setMethod_ConjGrad() {
    doSteepestDescent_or_doConjGrad_0or1 = 1;
  }
  void set_requireRigorousLineMinimization( bool val ) { 
    requireRigorousLineMinimization = val;
  }
  void enable_TolMaxGradComp( double tol ) {
    isEnabled_TolMaxGradComp = true;
    tolMaxGradComp = tol;
  }
  void enable_TolMaxMismatch( Mismatch *mismatch_, double tol ) {
    isEnabled_TolMaxMismatch = true;
    tolMaxMismatch = tol;
    mismatch = mismatch_;
  }
  void disable_TolMaxGradComp() {
    isEnabled_TolMaxGradComp = false;
  }
  void disable_TolMaxMismatch() {
    isEnabled_TolMaxMismatch = false;
    mismatch = NULL;
  }
  void enable_TolMaxPreconditionedGradComp( double tol ) {
    isEnabled_TolMaxPreconditionedGradComp = true;
    tolMaxPreconditionedGradComp = tol;
  }
  void setTrialStepLambda( double lam ) { trial_lambdamax = lam; }
  void setTrialStepDQmax( double dqmax ) { trial_dqmax = dqmax; } 
  void setLimitingParabolicMagnification( double lim ) { 
    findBracket.limitingParabolicMagnification = lim;
  }
private:
  MultiVarFunction *f;
  FindBracket<MultiVarFunction_1DProjection> findBracket;
  LineMinimizer_BrentMethod<MultiVarFunction_1DProjection> linemin;
  MultiVarFunction_1DProjection multiVarFunction_1DProjection;

  double lambda;
  bool reachedTol;
  int step;
  int step_CGreset;
  
  bool isEnabled_TolMaxGradComp;
  bool isEnabled_TolMaxMismatch;
  bool isEnabled_TolMaxPreconditionedGradComp;
  Mismatch *mismatch;
  double tolMaxGradComp;
  double tolMaxMismatch;
  double tolMaxPreconditionedGradComp;
  int outputPeriod;
  bool doSteepestDescent_or_doConjGrad_0or1;
  bool requireRigorousLineMinimization;
  double trial_dqmax;
  double trial_lambdamax;
  bool isPreconditioned;

  double r_dot_s_new;
  double r_dot_s_mid;
  double r_dot_s_old;
  std::vector<double> r;
  std::vector<double> s;
  std::vector<double> p;

  void stepAlongCurrentDirection();
  void setNextDirection_SteepestDescent();
  void setNextDirection_ConjGrad();
  bool toleranceSatisfied();
  void startMinimization( MultiVarFunction *f_ );
  void doStep();
  
};

#endif /*CONJUGATEGRADIENTMINIMIZER_H_*/
