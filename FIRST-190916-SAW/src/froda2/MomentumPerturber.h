#ifndef MOMENTUMPERTURBER_H_
#define MOMENTUMPERTURBER_H_

#include <vector>
class RigidUnitSystem;
#include "RigidUnitFitter.h"
#include "Vec3.h"

class MomentumPerturber
{
public:
  MomentumPerturber( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~MomentumPerturber();
  void setQ1();
  void determineDeltaQ();
  void perturb();

private:
  RigidUnitSystem *rigidUnitSystem;
  RigidUnitFitter rigidUnitFitter;
  std::vector<Vec3> Q1_rigidUnitPoints;
  bool isQ1set;
  bool readyToPerturb;  
};

#endif /*MOMENTUMPERTURBER_H_*/
