#ifndef RIGIDUNITFITTER_H_
#define RIGIDUNITFITTER_H_

class RigidUnitSystem;
#include <vector>
#include "Vec3.h"
#include "Fit.h"

class RigidUnitFitter
{
public:
	RigidUnitFitter( const RigidUnitSystem *rigidUnitSystem_ );
	virtual ~RigidUnitFitter();
  void setRigidUnits( const RigidUnitSystem *rigidUnitSystem_ );
  void calcFitToMeanPoints( const std::vector<Vec3> &target );
  void calcFitToRigidUnitPoints( const std::vector<Vec3> &absolutePositions );
  const Vec3& getFitRotation( size_t ru ) const { return fitRotation[ru]; }
  const Vec3& getFitTranslation( size_t ru ) const { return fitTranslation[ru]; }
private:
  const RigidUnitSystem *rigidUnitSystem;
  std::vector<Vec3> fitRotation;
  std::vector<Vec3> fitTranslation;
  /*std::vector<Vec3> ruBase;
  std::vector<Vec3> ruTarget;
  Fit fit;*/
};

#endif /*RIGIDUNITFITTER_H_*/
