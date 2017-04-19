#ifndef AREAPROXIMITY_H_
#define AREAPROXIMITY_H_

#include "global_defs.h"
#include "SiteID.h"
#include "ProximityMonitor.h"

class MolFramework;
class Froda;

class AreaProximity : public ProximityMonitor {
public:
  AreaProximity( Froda *froda_, MolFramework *structure_ );
  ~AreaProximity() {};
  
private:
  Froda *froda;
  bool isPairDisqualified( const ProximityAtom &atom1 , const ProximityAtom &atom2 ) const;
};

#endif /*AREAPROXIMITY_H_*/
