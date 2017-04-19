#ifndef TESTFOREXCLUSION_H_
#define TESTFOREXCLUSION_H_

#include "RigidUnitSystem.h"

class TestForExclusion
{
public:
	TestForExclusion( const RigidUnitSystem *rigidUnitSystem_ ) :
    rigidUnitSystem( rigidUnitSystem_ ) {}
	virtual ~TestForExclusion() {}
  bool operator()( int p1, int p2 ) const;
private:
  const RigidUnitSystem *rigidUnitSystem;
};

inline bool TestForExclusion::operator()( int p1, int p2 ) const {
  // exclude pair if p1 >= p2 to avoid self-counting
  // and double-counting.
  // Also, exclude pair if the two points belong
  // to the same rigid unit.
  return ( p1 >= p2 ||
           rigidUnitSystem->doPointsBelongToSameRigidUnit( p1, p2 ) );
}

#endif /*TESTFOREXCLUSION_H_*/
