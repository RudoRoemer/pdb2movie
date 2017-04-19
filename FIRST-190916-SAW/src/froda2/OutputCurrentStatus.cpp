#include "OutputCurrentStatus.h"
#include "PerturbRelaxCycle.h"
#include "RMSD.h"
#include "RigidUnitSystem.h"
#include "ConstraintEnforcingPotential.h"
#include "MinimizeSystem.h"
#include "PhaseSpacePathLengthIntegrator.h"
#include <iostream>
#include <iomanip>

OutputCurrentStatus::OutputCurrentStatus( 
    const RigidUnitSystem *rigidUnitSystem_,
    const PerturbRelaxCycle *cycle_,
    const MinimizeSystem *minim_,
    ConstraintEnforcingPotential *cep_ ) :
      rigidUnitSystem(rigidUnitSystem_),
      cycle(cycle_),
      minim(minim_),
      cep(cep_),
      path( NULL ),
      pathLengthEnabled( false ),
      pathMultiplier(0),
      outputPeriod( 1 )
{
  initialPoints = rigidUnitSystem->meanPositions();
  rmsdToInitial = new RMSD( &initialPoints, &rigidUnitSystem->meanPositions() );  
}

OutputCurrentStatus::~OutputCurrentStatus()
{
  delete rmsdToInitial;
}

void OutputCurrentStatus::display() {
  cout << fixed << std::setprecision(10) << cycle->getCycleCount() <<
       " " << minim->getFinalStepNum() << 
       " " << cep->maxMismatch() << 
       " " << (*rmsdToInitial)();
  if ( pathLengthEnabled ) {
    cout << " " << path->getSegmentLength() << 
            " " << path->getIntegratedPathLength() * pathMultiplier;
  }
  cout << "\n";
  if ( cycle->getCycleCount() % outputPeriod == 0 ) cout.flush();
}
