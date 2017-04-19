#ifndef GENERALIZEDCOORDS_H_
#define GENERALIZEDCOORDS_H_

#include <vector>
#include "Vec3.h"

class GeneralizedCoords
{
public:
	GeneralizedCoords() {}
	virtual ~GeneralizedCoords() {}
  std::vector<Vec3> centersRU;
  std::vector<Vec3> rotorsRU;
  size_t nRU() const { return centersRU.size(); }  
};


#endif /*GENERALIZEDCOORDS_H_*/
