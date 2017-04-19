#ifndef STERICS_H_
#define STERICS_H_

#include "global_defs.h"
#include "SiteID.h"
#include "ProximityMonitor.h"

class MolFramework;
class Froda;

class Sterics : public ProximityMonitor {
public:
  Sterics( Froda *froda_, MolFramework *structure_ );
  ~Sterics() {};
  double getMaxVdwRadius() const {return maxVdwRadius;}
  
private:
  double maxVdwRadius;
  Froda *froda;
  bool isPairDisqualified( const ProximityAtom &atom1 , const ProximityAtom &atom2 ) const;
};

#endif /*STERICS_H_*/
