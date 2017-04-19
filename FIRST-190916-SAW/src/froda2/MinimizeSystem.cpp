#include "MinimizeSystem.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"

MinimizeSystem::MinimizeSystem( 
    RigidUnitSystem *rigidUnitSystem_,
    ConstraintEnforcingPotential *cep_) :
  multiVarFunction( rigidUnitSystem_, cep_ ),
  rigidUnitSystem( rigidUnitSystem_ ),
  cep( cep_ ),
  NminimizationSteps( 200 ),
  nSteps( 0 )
{
  cgmin.registerObserver( this );
  cgmin.set_requireRigorousLineMinimization( false );
}

MinimizeSystem::~MinimizeSystem()
{
}

void MinimizeSystem::setToleranceCondition( string tolType, double tol ) {
  if ( tolType == "mismatch" ) {
    cgmin.enable_TolMaxMismatch( cep->mismatch(), tol );
  }
  else if ( tolType == "maxPreconditionedGradComponent" ) {
    cgmin.enable_TolMaxPreconditionedGradComp( tol );
  }
  else {
    cout << "Error: tolType \"" << tolType << "\" not recognized." << endl;
    exit(0);
  }
}

void MinimizeSystem::receiveNotification( Observable *obs ) {
  nSteps++;
  notifyObservers();
}

void MinimizeSystem::minimize() {
  nSteps = 0;
  cgmin.setMethod_SteepestDescent();
  cgmin.minimize( &multiVarFunction, 2 );
  cgmin.setMethod_ConjGrad();
  rigidUnitSystem->collapseRotors();
  /*
  cgmin.minimize( &multiVarFunction, min( 400, NminimizationSteps - nSteps ) );
  rigidUnitSystem->collapseRotors();
  */
  int nRemainingSteps;
  int nCGsteps = 400;
  while ( !cgmin.isTolReached() && nSteps < NminimizationSteps ) {
    nRemainingSteps = NminimizationSteps - nSteps;
    cgmin.minimize( &multiVarFunction, min( nCGsteps, nRemainingSteps ) );
    rigidUnitSystem->collapseRotors();
  }
}
