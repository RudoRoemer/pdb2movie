#include "MolFramework.h"
#include "SiteID.h"
#include "Parameters.h"
#include "global_defs.h"
#include "Grid.h"

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   set up the coarse gridding for, well, whatever we're using it for
////////////////////////////////////////////////////////////////////////////////
Grid::Grid(MolFramework &structure_, double myGridLength_ ) : structure(&structure_),
							      myGridLength(myGridLength_) {
            
  // TODO: make sure electron density, restart, etc. don't require a bigger box; 
  //       expand structure->max_coords if needed
  
  // pad the coarse grid to allow for larger movements
  structure->max_coords[X] += 30.0; structure->min_coords[X] -= 30.0;
  structure->max_coords[Y] += 30.0; structure->min_coords[Y] -= 30.0;
  structure->max_coords[Z] += 30.0; structure->min_coords[Z] -= 30.0;
  cout << "ED Coarse grid minimum: " << structure->min_coords[X] << ' ';
  cout << structure->min_coords[Y] << ' ' << structure->min_coords[Z] << endl;
  cout << "ED Coarse grid maximum: " << structure->max_coords[X] << ' ';
  cout << structure->max_coords[Y] << ' ' << structure->max_coords[Z] << endl;
  
  myXLength = structure->max_coords[X] - structure->min_coords[X];
  myYLength = structure->max_coords[Y] - structure->min_coords[Y];
  myZLength = structure->max_coords[Z] - structure->min_coords[Z];
  
  myXCells = (int) ceil (myXLength/myGridLength);
  myYCells = (int) ceil (myYLength/myGridLength);
  myZCells = (int) ceil (myZLength/myGridLength);

  myXCells +=2;
  myYCells +=2;
  myZCells +=2;

  structure->x_cells = myXCells;
  structure->y_cells = myYCells;
  structure->z_cells = myZCells;
 
  myGrid.resize(myXCells*myYCells*myZCells);
  secondGrid.resize(myXCells*myYCells*myZCells);

  update();

  return;
}

Grid::~Grid()
{
}

////////////////////////////////////////////////////////////////////////////////
// Description: Update coarse grid used during steric detection
// This function relies on the atomic positions being found in the
// MolFramework structure.
////////////////////////////////////////////////////////////////////////////////
void Grid::update(){ 
  vector<vector<SiteID> >::iterator igrid, endgrid;
  endgrid = myGrid.end();
  igrid = myGrid.begin();
  while (igrid != endgrid) {
    (*igrid++).clear();
  }
  endgrid = secondGrid.end();
  igrid = secondGrid.begin();
  while (igrid != endgrid) {
    (*igrid++).clear();
  }
  
  for(unsigned int siteNumber = 1; siteNumber<= structure->total_sites;siteNumber++){

    structure->site_info[siteNumber].grid_X = (int) ((structure->site_info[siteNumber].coords[X] - structure->min_coords[X])/myGridLength);
    structure->site_info[siteNumber].grid_Y = (int) ((structure->site_info[siteNumber].coords[Y] - structure->min_coords[Y])/myGridLength);
    structure->site_info[siteNumber].grid_Z = (int) ((structure->site_info[siteNumber].coords[Z] - structure->min_coords[Z])/myGridLength);

    structure->site_info[siteNumber].grid_X++;
    structure->site_info[siteNumber].grid_Y++;
    structure->site_info[siteNumber].grid_Z++;

    // Error check.
    //////////////////////////////////////////////////////////////////////
    if( structure->site_info[siteNumber].grid_X > myXCells || structure->site_info[siteNumber].grid_X < 0 || 
        structure->site_info[siteNumber].grid_Y > myYCells || structure->site_info[siteNumber].grid_Y < 0 || 
        structure->site_info[siteNumber].grid_Z > myZCells || structure->site_info[siteNumber].grid_Z < 0 ){
      cout << " ERROR: Coordinates of atom not in 3D grid. (" << structure->site_info[siteNumber].grid_X 
           << " " << structure->site_info[siteNumber].grid_Y << " " << structure->site_info[siteNumber].grid_Z << ") = " 
           << structure->gridPosition(siteNumber) << endl;
      exit(1);
    }

    myGrid[structure->gridPosition(siteNumber)].push_back(siteNumber);
    // now update secondGrid, putting this atom into several different slots
    int thisX = structure->site_info[siteNumber].grid_X;
    int thisY = structure->site_info[siteNumber].grid_Y;
    int thisZ = structure->site_info[siteNumber].grid_Z;
    int otherX, otherY, otherZ, otherN;
    for (otherX = thisX-2; otherX < thisX+3; otherX++) {
      for (otherY = thisY-2; otherY < thisY+3; otherY++) {
        for (otherZ = thisZ-2; otherZ < thisZ+3; otherZ++) {

          if (otherX < 0 || otherY < 0 || otherZ < 0 ) continue;
          if (otherX >= myXCells || otherY >= myYCells || otherZ >= myZCells ) continue;

          otherN = structure->gridPosition( otherX, otherY, otherZ );
          secondGrid[otherN].push_back(siteNumber);
        }
      }
    }         
  }
}

void Grid::getNearbyAtoms( SiteID atom1, vector<SiteID>& nearbySiteIDs ) {
  nearbySiteIDs.clear();
  SiteID atom2;
  int gridsToCheck[27];
  
  structure->getGrids( atom1, gridsToCheck);
  for ( int p = 0; p < 27; p++ ) {
    int thisGrid = gridsToCheck[p];
    int thisGridSize = myGrid[thisGrid].size();
    for ( int q=0; q < thisGridSize; q++) {
      atom2 = myGrid[thisGrid][q];
      if ( atom2 <= atom1 ) continue;
      //if ( froda->bonded( froda->isInGhost[atom1], froda->isInGhost[atom2] ) ) continue;
      nearbySiteIDs.push_back( atom2 );
    }
  } 
}

void Grid::getAtomsWithin2GridPoints( SiteID atom1, vector<SiteID>& nearbySiteIDs ) {
  nearbySiteIDs.clear();
  int thisX = structure->site_info[atom1].grid_X;
  int thisY = structure->site_info[atom1].grid_Y;
  int thisZ = structure->site_info[atom1].grid_Z;
  int otherX, otherY, otherZ, otherN;
  SiteID atom2;
  //now check in the cells around;
  for (otherX = thisX-2; otherX < thisX+3; otherX++) {
    for (otherY = thisY-2; otherY < thisY+3; otherY++) {
      for (otherZ = thisZ-2; otherZ < thisZ+3; otherZ++) {

        if (otherX < 0 || otherY < 0 || otherZ < 0 ) continue;
        if (otherX >= myXCells || otherY >= myYCells || otherZ >= myZCells ) continue;

        otherN = structure->gridPosition( otherX, otherY, otherZ );

        int otherSize = myGrid.at(otherN).size();
        for ( int j=0; j < otherSize; j++) {
          atom2 = myGrid.at(otherN).at(j);
          if ( atom2 <= atom1 ) continue;
          nearbySiteIDs.push_back( atom2 );
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Given a location in Cartesian space, returns a list of the atoms that are
//   within two cells of that point on the coarse grid
////////////////////////////////////////////////////////////////////////////////
void Grid::getAtomsWithin2LengthsOfPoint( Vector location, vector<SiteID>& nearbySiteIDs ) {
  nearbySiteIDs.clear();
  int thisX = (int) ((location.x - structure->min_coords[X])/myGridLength);
  int thisY = (int) ((location.y - structure->min_coords[Y])/myGridLength);
  int thisZ = (int) ((location.z - structure->min_coords[Z])/myGridLength);
  /*int otherX, otherY, otherZ, otherN;
  SiteID atom;
  //now check in the cells around;
  for (otherX = thisX-2; otherX < thisX+3; otherX++) {
    for (otherY = thisY-2; otherY < thisY+3; otherY++) {
      for (otherZ = thisZ-2; otherZ < thisZ+3; otherZ++) {

        if (otherX < 0 || otherY < 0 || otherZ < 0 ) continue;
        if (otherX >= myXCells || otherY >= myYCells || otherZ >= myZCells ) continue;

        otherN = structure->gridPosition( otherX, otherY, otherZ );

        int otherSize = myGrid.at(otherN).size();
        for ( int j=0; j < otherSize; j++) {
          atom = myGrid.at(otherN).at(j);
          nearbySiteIDs.push_back( atom );
        }
      }
    }
  }*/
  nearbySiteIDs = secondGrid.at(structure->gridPosition( thisX, thisY, thisZ ));
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Given a location in Cartesian space, returns the cell in the coarse grid
//   to which the point belongs.
////////////////////////////////////////////////////////////////////////////////
vector< int > Grid::gridCellsOfLocation( Vector location ) {
  int thisX = (int) ((location.x - structure->min_coords[X])/myGridLength);
  int thisY = (int) ((location.y - structure->min_coords[Y])/myGridLength);
  int thisZ = (int) ((location.z - structure->min_coords[Z])/myGridLength);

  vector< int > gridback;
  gridback.resize(3);
  gridback.at(0) = thisX;
  gridback.at(1) = thisY;
  gridback.at(2) = thisZ;

  return gridback;
}


