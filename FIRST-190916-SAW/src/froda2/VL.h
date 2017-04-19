#ifndef VL_H_
#define VL_H_

#include "Observable.h"
#include "VerletList.h"
#include "TestForExclusion.h"
#include "Repulsion.h"
class RigidUnitSystem;
#include <list>
#include <vector>

class VL : public Observer,
           public Observable
{
public:
	VL( const RigidUnitSystem *rigidUnitSystem,
      const Repulsion &repulsion );
	virtual ~VL();
  
  const std::vector<VL_NeighborInfo>& getProximityList( int point1 ) const {
    return verletList.getProximityList( point1 );
  }

  void receiveNotification( Observable *observable ) {
    update();
  }
private:  
  TestForExclusion testForExclusion;
  VerletList<const TestForExclusion, const Repulsion> verletList;
  void update();  
};

#endif /*VL_H_*/
