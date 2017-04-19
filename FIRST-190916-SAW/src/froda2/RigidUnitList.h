#ifndef RIGIDUNITLIST_H_
#define RIGIDUNITLIST_H_

#include "NeighborTable.h"
#include <vector>
#include <set>
#include <map>

class RigidUnitList
{
public:
  RigidUnitList(
    const std::vector<int>& mapPtoRC,
    const NeighborTable *neighborTable );

	virtual ~RigidUnitList();

  const std::vector< std::vector<int> >& getAssignment_RUtoPlist() const;

  friend std::ostream& operator<< ( std::ostream& os, const RigidUnitList& rigidUnitList );

private:
  void extendUnitsByOneNeighbor( const NeighborTable *neighborTable );
  void finalize();
  void assignPtoRC( const std::vector<int>& mapPtoRC );
  std::map< int, std::set< int > > mapRCtoP;
  std::vector< std::vector<int> > mapRUtoP_vector;
};

#endif /*RIGIDUNITLIST_H_*/
