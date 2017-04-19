#include "GlobalMotionRemover.h"
#include "Fit.h"
#include "Rotator.h"
#include "RigidUnitSystem.h"

GlobalMotionRemover::GlobalMotionRemover(
      RigidUnitSystem *rigidUnitSystem_ ) :
  rigidUnitSystem( rigidUnitSystem_ )
{
}

GlobalMotionRemover::~GlobalMotionRemover()
{
}

void GlobalMotionRemover::setCurrentPointsAsTarget() {
  fit.setTargetAbsolutePoints( rigidUnitSystem->meanPositions() );
}

void GlobalMotionRemover::fitCurrentPointsToTarget() {
  //calculate the current center 
  Vec3 center(0,0,0); 
  size_t N = rigidUnitSystem->meanPositions().size();
  for ( size_t i=0; i < N; i++ ) { 
    center += rigidUnitSystem->meanPositions()[i];
  }
  center /= N;
  
  //calculate the global fit
  fit.setSourceAbsolutePoints( rigidUnitSystem->meanPositions() );
  fit.simpleFit();
  
  // Apply global translation to each rigid unit
  Vec3 deltaCenter = fit.getFitCenter() - center;
  size_t nRU = rigidUnitSystem->nRigidUnits();
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    rigidUnitSystem->addToCenter( ru, deltaCenter );
  }

  // Apply global rotation to each rigid unit
  Rotator rotator( fit.getFitRotor() );
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    rigidUnitSystem->rotate( ru, rotator, fit.getFitCenter() );
  }
  
  rigidUnitSystem->collapseRotors();
  rigidUnitSystem->update();
}
