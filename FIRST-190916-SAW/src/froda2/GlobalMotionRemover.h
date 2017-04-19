#ifndef GLOBALMOTIONREMOVER_H_
#define GLOBALMOTIONREMOVER_H_

using namespace std;

class RigidUnitSystem;
#include "Fit.h"

class GlobalMotionRemover
{
public:
  GlobalMotionRemover( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~GlobalMotionRemover();
  void setCurrentPointsAsTarget();
  void fitCurrentPointsToTarget();
private:
  RigidUnitSystem *rigidUnitSystem;
  Fit fit;
};

#endif /*GLOBALMOTIONREMOVER_H_*/
