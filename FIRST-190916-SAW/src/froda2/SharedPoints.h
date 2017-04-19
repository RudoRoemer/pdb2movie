#ifndef SHAREDPOINTS_H_
#define SHAREDPOINTS_H_

class RigidUnitSystem;
#include <vector>
#include "GeneralizedCoords.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "Vec3.h"

class SharedPoints : public EnergyTerm, 
                     public GradientTerm,
                     public MismatchTerm
{
public:
	SharedPoints( const RigidUnitSystem *rigidUnitSystem_ );
	~SharedPoints();
  double energy();
  void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint );
  double mismatch();
private:
  std::vector<int> sharedPoints;
  const RigidUnitSystem *rigidUnitSystem;
};

#endif /*SHAREDPOINTS_H_*/
