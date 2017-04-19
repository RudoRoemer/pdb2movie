#ifndef OUTPUTCURRENTSTATUS_H_
#define OUTPUTCURRENTSTATUS_H_

#include "PerturbRelaxCycle.h"
class RMSD;
#include "Vec3.h"
#include <vector>
class RigidUnitSystem;
class MinimizeSystem;
class ConstraintEnforcingPotential;
class PhaseSpacePathLengthIntegrator;

class OutputCurrentStatus
{
public:
  OutputCurrentStatus( 
      const RigidUnitSystem *rigidUnitSystem_,
      const PerturbRelaxCycle *cycle_,
      const MinimizeSystem *minim_,
      ConstraintEnforcingPotential *cep_ );
  virtual ~OutputCurrentStatus();
  
  void display();
  
  void setPeriod( int outputPeriod_ ) { outputPeriod = outputPeriod_; }
  void enablePathLength( const PhaseSpacePathLengthIntegrator *path_, double pathMultiplier_ ) {
    pathLengthEnabled = true;
    path = path_;
    pathMultiplier = pathMultiplier_; 
  }
  
private:
  const RigidUnitSystem *rigidUnitSystem;
  const PerturbRelaxCycle *cycle;
  const MinimizeSystem *minim;
  ConstraintEnforcingPotential *cep;
  RMSD *rmsdToInitial;
  const PhaseSpacePathLengthIntegrator *path;
  bool pathLengthEnabled;
  double pathMultiplier;
  int outputPeriod;
  std::vector<Vec3> initialPoints;
};

#endif /*OUTPUTCURRENTSTATUS_H_*/
