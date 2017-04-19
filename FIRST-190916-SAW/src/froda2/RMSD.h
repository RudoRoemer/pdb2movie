#ifndef RMSD_H_
#define RMSD_H_

#include <vector>
#include "Vec3.h"

class RMSD
{
public:
	RMSD( const std::vector<Vec3> *points1_, const std::vector<Vec3> *points2_ );
	virtual ~RMSD();
  double calc() const;
  double operator()() const { return calc(); }
private:
  const std::vector<Vec3> *points1;
  const std::vector<Vec3> *points2;
};

#endif /*RMSD_H_*/
