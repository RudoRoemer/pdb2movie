#ifndef ROPES_H_
#define ROPES_H_

#include <vector>
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "GeneralizedCoords.h"
#include "RopesAssignment.h"
#include "Vec3.h"
class RigidUnitSystem;

class Ropes : public EnergyTerm, 
              public GradientTerm,
              public MismatchTerm
{
public:
	Ropes( const RopesAssignment& ropes_, 
         const RigidUnitSystem *rigidUnitSystem_ );
  
	virtual ~Ropes() {}
  
  double energy();
  void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint );
  double mismatch();

private:
  RopesAssignment ropes;
  const RigidUnitSystem *rigidUnitSystem;
};

#endif /*ROPES_H_*/
