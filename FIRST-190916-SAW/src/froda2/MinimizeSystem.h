#ifndef MINIMIZESYSTEM_H_
#define MINIMIZESYSTEM_H_

class RigidUnitSystem;
class ConstraintEnforcingPotential;
#include "ConjugateGradientMinimizer.h"
#include "MultiVarFunction_Adapter_RigidUnits.h"
#include "Observable.h"

class MinimizeSystem : public Observer, public Observable
{
public:
  MinimizeSystem( 
      RigidUnitSystem *rigidUnitSystem_,
      ConstraintEnforcingPotential *cep_ );
  virtual ~MinimizeSystem();
  void minimize();
  void setToleranceCondition( string tolType, double tol );
  void setNminimizationSteps( int NminimizationSteps_ ) { 
    NminimizationSteps = NminimizationSteps_;
  }
  int getFinalStepNum() const {
    return nSteps;
  }
  int getNumCompletedIterations() const {
    return nSteps;
  }
  void receiveNotification( Observable *obs );
private:
  MultiVarFunction_Adapter_RigidUnits multiVarFunction;
  ConjugateGradientMinimizer cgmin;
  RigidUnitSystem *rigidUnitSystem;
  ConstraintEnforcingPotential *cep;
  int NminimizationSteps;
  int nSteps;
};

#endif /*MINIMIZESYSTEM_H_*/
