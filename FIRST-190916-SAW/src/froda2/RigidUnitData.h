#ifndef RIGIDUNITDATA_H_
#define RIGIDUNITDATA_H_

#include "Vec3.h"
#include <vector>
#include "AbstractRigidUnitData.h"

class RigidUnitData : public AbstractRigidUnitData
{
public:
  RigidUnitData() {}
  virtual ~RigidUnitData() {}  
  std::vector<Vec3> meanPositions;
};

#endif /*RIGIDUNITDATA_H_*/
