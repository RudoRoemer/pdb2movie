#include "SymmetricCEPBuilder.h"
#include "SymmetryMatrices.h"
#include "SymmetryEnforcer.h"
#include "RigidUnitSystem.h"

SymmetricCEPBuilder::SymmetricCEPBuilder(RigidUnitSystem *rigidUnitSystem, 
  SymmetryMatrices *symmetryMatrices_) : ConstraintEnforcingPotentialBuilder(rigidUnitSystem), 
  symmetryMatrices(symmetryMatrices_) {
  setSymmetryEnforcerEnergy();
}

SymmetricCEPBuilder::~SymmetricCEPBuilder() {
}

void SymmetricCEPBuilder::setSymmetryEnforcerEnergy() {
  SymmetryEnforcer *symmetryEnforcer = new SymmetryEnforcer(rigidUnitSystem,
                                                 symmetryMatrices);
  rigidUnitSystem->registerObserver(symmetryEnforcer);
  energy->addTerm(symmetryEnforcer);
  gradient->addTerm(symmetryEnforcer);
  mismatch->addTerm(symmetryEnforcer);
}


  
