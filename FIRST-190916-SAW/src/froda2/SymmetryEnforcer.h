// SymmetryEnforcer.h
// Energy term to enforce symmetry
// Craig Jolley, January 2008
////////////////////////////////////////////////////////////////////////////////

#ifndef SYMMETRY_ENFORCER_H_
#define SYMMETRY_ENFORCER_H_

#include <vector>
#include "GeneralizedCoords.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "Vec3.h"
#include "SymmetryMatrices.h"
#include "Observable.h"
#include "RigidUnitSystem.h"

class SymmetryEnforcer : public EnergyTerm,
                         public GradientTerm,
                         public MismatchTerm,
                         public Observer {
public:
  SymmetryEnforcer(const RigidUnitSystem *rigidUnitSystem_,
                   const SymmetryMatrices *symmetryMatrices_);
  ~SymmetryEnforcer();
  double energy();
  void addToGradient(std::vector<Vec3> &dV_dr_rigidUnitPoint,
                     std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint);
  double mismatch();
  virtual void receiveNotification( Observable *observable ) {
    if (rigidUnitSystem->AbsolutePositionsChanged()) {
      updateAvgPositions();
    }
  }
private:
  const RigidUnitSystem *rigidUnitSystem;
  const SymmetryMatrices *symmetryMatrices;
  vector <Vec3> averagePositions;
  vector <double> d2V_dx2;
  vector <double> d2V_dy2;
  vector <double> d2V_dz2;
  Vec3 transform(Vec3 v, Matrix m);
  Vec3 getAvgPosition(int p);
  void updateAvgPositions();
  void initializeSecondDerivative();
};


#endif /*  SYMMETRY_ENFORCER_H_*/
