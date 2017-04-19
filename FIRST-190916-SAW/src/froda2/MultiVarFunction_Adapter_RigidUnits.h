#ifndef MULTIVARFUNCTION_ADAPTER_RIGIDUNITS_H_
#define MULTIVARFUNCTION_ADAPTER_RIGIDUNITS_H_

#include <vector>
#include "MultiVarFunction.h"
class RigidUnitSystem;
class ConstraintEnforcingPotential;

class MultiVarFunction_Adapter_RigidUnits : public MultiVarFunction
{
public:
	MultiVarFunction_Adapter_RigidUnits(
    RigidUnitSystem *rigidUnitSystem_,
    ConstraintEnforcingPotential *cep_ );
  MultiVarFunction_Adapter_RigidUnits();
	virtual ~MultiVarFunction_Adapter_RigidUnits();

  void setRigidUnits(
    RigidUnitSystem *rigidUnitSystem_,
    ConstraintEnforcingPotential *cep_ );  
  void getQ( std::vector<double> &q );
  void setQ( const std::vector<double> &q );
  double eval();
  void getNegGrad( std::vector<double> &r );
  void getPreconditionedNegGrad( std::vector<double> &s, bool &isPreconditioned );
private:
  RigidUnitSystem *rigidUnitSystem;
  ConstraintEnforcingPotential *cep;
};

#endif /*MULTIVARFUNCTION_ADAPTER_RIGIDUNITS_H_*/
