#ifndef ABSTRACTRIGIDUNITDATA_H_
#define ABSTRACTRIGIDUNITDATA_H_

#include "Vec3.h"
#include <vector>

class AbstractRigidUnitData
{
public:
  AbstractRigidUnitData() {}
  virtual ~AbstractRigidUnitData()=0;
  std::vector<Vec3> rotors;
  std::vector<Vec3> centers;
  std::vector<Vec3> basePositions;
  std::vector<Vec3> absolutePositions;
  std::vector<double> radius;
  std::vector<char> hasZeroRadius;
};

inline AbstractRigidUnitData::~AbstractRigidUnitData() {}

#endif /*ABSTRACTRIGIDUNITDATA_H_*/
