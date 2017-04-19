#ifndef RIGIDUNITLOOKUP_H_
#define RIGIDUNITLOOKUP_H_

#include "AbstractRigidUnitLookup.h"
#include <vector>
#include <cstddef>

class RigidUnitLookup: public AbstractRigidUnitLookup
{
public:
  RigidUnitLookup( const std::vector<IDlist>& RUtoPlist_, size_t nP )
  {
    setLookup_RUtoPlist( RUtoPlist_, nP );
  }
  RigidUnitLookup() {}
  virtual ~RigidUnitLookup() {}

  void setLookup_RUtoPlist( 
    const std::vector<IDlist>& RUtoPlist_, size_t nP );

  const IDlist& getPlistFromRU( int ru ) const {
    return RUtoPlist[ru];
  }
  const IDlist& getRUlistFromP( int p ) const {
    return PtoRUlist[p];
  }
  const IDlist& getRUPlistFromP( int p ) const {
    return PtoRUPlist[p];
  }
  const int getPfromRUP( int rup ) const {
    return RUPtoP[rup];
  }
  const size_t nPoints() const {
    return PtoRUlist.size();
  }
  
  bool doPointsBelongToSameRigidUnit( int p1, int p2 ) const;

protected:
  std::vector<IDlist> RUtoPlist;
  std::vector<IDlist> PtoRUlist;

  std::vector<IDlist> PtoRUPlist;
  std::vector<int> RUPtoP;   

};

inline bool RigidUnitLookup::doPointsBelongToSameRigidUnit( int p1, int p2 ) const {
  IDlist::const_iterator first1 = getRUlistFromP(p1).begin();
  IDlist::const_iterator last1 = getRUlistFromP(p1).end();
  IDlist::const_iterator first2 = getRUlistFromP(p2).begin();
  IDlist::const_iterator last2 = getRUlistFromP(p2).end();
  
  while ( first1 != last1 && first2 != last2 ) {
    if ( *first1 < *first2 ) first1++;
    else if ( *first2 < *first1 ) first2++;
    else {
      // We know *first1==*first2.  
      // This means that the two points belong
      // to the same rigid unit.
      return true;
    }
  }
  
  return false;
}

#endif /*RIGIDUNITLOOKUP_H_*/
