#include "ConstraintEnforcingPotentialBuilder.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "SharedPoints.h"
#include "VL.h"
#include "Repulsion.h"
#include "OverlapList.h"
#include "Repulsion.h"
#include "RopesAssignment.h"
#include "Ropes.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include "InequalityAngleConstraints.h"

ConstraintEnforcingPotentialBuilder::ConstraintEnforcingPotentialBuilder(
  RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ )
{
  
  energy = new Energy();
  gradient = new Gradient( rigidUnitSystem );
  mismatch = new Mismatch();
  
  rigidUnitSystem->registerObserver( energy );
  rigidUnitSystem->registerObserver( gradient );
  rigidUnitSystem->registerObserver( mismatch );
  
}

ConstraintEnforcingPotentialBuilder::~ConstraintEnforcingPotentialBuilder()
{
}

ConstraintEnforcingPotential* ConstraintEnforcingPotentialBuilder::getConstraintEnforcingPotential() {
  ConstraintEnforcingPotential *cep = new ConstraintEnforcingPotential( 
    verletList,
    energy,
    gradient,
    mismatch
  );
  return cep;
}

void ConstraintEnforcingPotentialBuilder::setSharedPointsEnergy() {
  SharedPoints *sharedPoints = new SharedPoints( rigidUnitSystem ); 
  energy->addTerm( sharedPoints );
  gradient->addTerm( sharedPoints );
  mismatch->addTerm( sharedPoints );
}

void ConstraintEnforcingPotentialBuilder::setRopesEnergy( 
    const RopesAssignment* ropesAssignment ) {
  Ropes *ropes = new Ropes( *ropesAssignment, rigidUnitSystem );
  energy->addTerm( ropes );
  gradient->addTerm( ropes );
  mismatch->addTerm( ropes );
}

void ConstraintEnforcingPotentialBuilder::setOverlapEnergy( Repulsion *repulsion ) {
  verletList = new VL( rigidUnitSystem, *repulsion );
  rigidUnitSystem->registerObserver( verletList );
  OverlapList *overlapList = new OverlapList( rigidUnitSystem, *verletList );
  verletList->registerObserver( overlapList );
  energy->addTerm( overlapList );
  gradient->addTerm( overlapList );
  mismatch->addTerm( overlapList );
}

void ConstraintEnforcingPotentialBuilder::setInequalityAngleConstraints( InequalityAngleConstraints *con ) {
  energy->addTerm( con );
  gradient->addTerm( con );
  mismatch->addTerm( con );
}
