#ifndef FITPERTURBER_H_
#define FITPERTURBER_H_

class RigidUnitSystem;
#include "RigidUnitFitter.h"
#include "Vec3.h"
#include <vector>

class FitPerturber
{
public:
  FitPerturber( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~FitPerturber();
  void setMeanPointsTarget( const std::vector<Vec3>* meanPointsTarget_ );
  void perturb();
private:
  RigidUnitSystem *rigidUnitSystem;
  RigidUnitFitter fitter;
  const std::vector<Vec3>* meanPointsTarget;
};

#endif /*FITPERTURBER_H_*/
