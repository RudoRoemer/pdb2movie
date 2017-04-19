#include "AreaProximity.h"
#include "Froda.h"
#include "MolFramework.h"
#include "ProximityMonitor.h"
#include "SiteID.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
AreaProximity::AreaProximity( Froda *froda_, MolFramework *structure ) : ProximityMonitor(),
									 froda(froda_) {

  if( parameters.verbose >= 2 ){
    cout << "Setting up Area Calcs Proximity Monitor " << endl;
  }
  
  double cushion = 2;
  double minSpannerBottomEdgeLength = 8.0; //6.0 for distance + 2 for cushion
  double beta = 3.9;
  double c=5.0;
  double unitDistance = minSpannerBottomEdgeLength/c;
  if ( c*unitDistance < minSpannerBottomEdgeLength ) {
    cout << "Warning, bottom level neighbor distance in spanner is too small for looking up possible contacts" << endl;
  }
  
  //store each atom that will be monitored
  for ( SiteID siteID = 1; siteID<=froda->nTotalSites; siteID++ ) {
    double cutoffDist = 6.0;
    //store the atom
    insert( &froda->currentPos[siteID], cutoffDist, siteID );
  }
  
  commit( beta, c, unitDistance, cushion );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool AreaProximity::isPairDisqualified( const ProximityAtom &atom1 , const ProximityAtom &atom2 ) const {

  SiteID id1 = atom1.getID();
  SiteID id2 = atom2.getID();

  if ( id2 == id1 ) {
    return true; //avoid self-counting
  }

  return false;  
}
