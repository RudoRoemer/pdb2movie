#include "Sterics.h"
#include "Froda.h"
#include "MolFramework.h"
#include "ProximityMonitor.h"
#include "SiteID.h"

extern Parameters parameters;

Sterics::Sterics( Froda *froda_, MolFramework *structure ) : ProximityMonitor(),
                                                             froda(froda_) {

  if( parameters.verbose >= 2 ){
    cout << "Setting up Sterics Proximity Monitor " << endl;
  }
  
  //find the maximum VDW radius
  maxVdwRadius=0;
  for ( SiteID siteID = 1; siteID<=froda->nTotalSites; siteID++ ) {
    double vdwRadius = structure->site_info[siteID].vdw_radius;
    if ( vdwRadius > maxVdwRadius ) maxVdwRadius=vdwRadius;
  }

  double cushion = 2;
  double minSpannerBottomEdgeLength = 2*maxVdwRadius + cushion;
  double beta = 3.9;
  double c=4.0;
  double unitDistance = minSpannerBottomEdgeLength/c;
  if ( c*unitDistance < minSpannerBottomEdgeLength ) {
    cout << "Warning, bottom level neighbor distance in spanner is too small for looking up possible contacts" << endl;
  }
  
  //store each atom that will be monitored
  for ( SiteID siteID = 1; siteID<=froda->nTotalSites; siteID++ ) {
    double vdwRadius = structure->site_info[siteID].vdw_radius;
    double cutoffDist = vdwRadius + maxVdwRadius;
    //store the atom
    insert( &froda->currentPos[siteID], cutoffDist, siteID );
  }
  
  commit( beta, c, unitDistance, cushion );
}

bool Sterics::isPairDisqualified( const ProximityAtom &atom1 , const ProximityAtom &atom2 ) const {
  SiteID id1 = atom1.getID();
  SiteID id2 = atom2.getID();
  //if ( id2 <= id1 ) return true; //avoid double-counting, use only i<j
  //SAW jan 10 07: moving check to Mismatch
  if ( id2 == id1 ) return true; //avoid i == j
  //throw out pairs within the same ghost (which includes
  //1st and 2nd neighbors)
  if ( froda->bonded( froda->frodaAtom.at(id1).ghostlist, froda->frodaAtom.at(id2).ghostlist ) ) return true;
  return false;  
}
