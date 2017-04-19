#ifndef SYMMETRIC_CEP_BUILDER_H_
#define SYMMETRIC_CEP_BUILDER_H_

#include "ConstraintEnforcingPotentialBuilder.h"

class SymmetryMatrices;

class SymmetricCEPBuilder : public ConstraintEnforcingPotentialBuilder {
public:
  SymmetricCEPBuilder(RigidUnitSystem *rigidUnitSystem, 
    SymmetryMatrices *symmetryMatrices_);
  virtual ~SymmetricCEPBuilder();
private:
  void setSymmetryEnforcerEnergy();
  SymmetryMatrices *symmetryMatrices;
};

#endif /* SYMMETRIC_CEP_BUILDER_H_ */
