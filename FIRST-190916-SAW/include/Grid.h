#ifndef GRID_H_
#define GRID_H_

#include "MolFramework.h"
#include "SiteID.h"
#include "global_defs.h"

class Grid {

public:

  Grid(MolFramework &structure, double myGridLength );
  ~Grid();
  void getNearbyAtoms( SiteID atom1, vector<SiteID>& nearbySiteIDs );
  void getAtomsWithin2GridPoints( SiteID atom1, vector<SiteID>& nearbySiteIDs );
  void getAtomsWithin2LengthsOfPoint( Vector location, vector<SiteID>& nearbySiteIDs );
  vector< int > gridCellsOfLocation( Vector location );
  void update();
  
  MolFramework *structure;
  MolFramework *restart;
  MolFramework *target;

  double myXLength, myYLength, myZLength;
  int myXCells, myYCells, myZCells;
  vector< vector<SiteID> > myGrid;
  
  // contains atoms at each grid point 
  vector< vector<SiteID> > secondGrid;
  
  // contains atoms within two grid points of each point 
  double myGridLength;
};

#endif /*GRID_H_*/
