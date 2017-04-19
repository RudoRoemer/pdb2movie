#ifndef ABSTRACTRIGIDUNITLOOKUP_H_
#define ABSTRACTRIGIDUNITLOOKUP_H_

#include <vector>
#include <cstddef>

class AbstractRigidUnitLookup
{
public:
  AbstractRigidUnitLookup() {}
  virtual ~AbstractRigidUnitLookup()=0;
  
  typedef std::vector<int> IDlist;
  
  const IDlist& getRUPlistFromRU( int ru ) const {
    return RUtoRUPlist[ru];
  }
  int getRUfromRUP( int rup ) const {
    return RUPtoRU[rup];
  }
  size_t nRigidUnitPoints() const {
    return RUPtoRU.size();
  }
  size_t nRigidUnits() const {
    return RUtoRUPlist.size();
  }
protected:
  std::vector<IDlist> RUtoRUPlist;
  std::vector<int> RUPtoRU;  
};

inline AbstractRigidUnitLookup::~AbstractRigidUnitLookup() {}

#endif /*ABSTRACTRIGIDUNITLOOKUP_H_*/
