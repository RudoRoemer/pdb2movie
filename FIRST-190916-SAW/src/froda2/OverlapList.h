#ifndef OVERLAPLIST_H_
#define OVERLAPLIST_H_

#include "Observable.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include <vector>
#include "VL.h"
#include "Vec3.h"
class RigidUnitSystem;

class Overlap {
public:
  int p1;
  int p2;
  double overlapDist;
  Vec3 delta21;
};

class OverlapList : public EnergyTerm, 
                    public GradientTerm,
                    public MismatchTerm,
                    public Observer
{
public:
	OverlapList( const RigidUnitSystem *rigidUnitSystem_,
               const VL &verletList_ );
	virtual ~OverlapList() {}
  double energy();
  void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
                      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint );
  double mismatch();
  void receiveNotification( Observable *observable ) {
    update();
  }
private:
  const RigidUnitSystem *rigidUnitSystem;
  const VL &verletList;
  std::vector<Overlap> overlapList;

  void update();
};

#endif /*OVERLAPLIST_H_*/
