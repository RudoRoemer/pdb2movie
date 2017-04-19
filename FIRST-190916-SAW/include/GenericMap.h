// Author name:  Craig Jolley
// Created:      19 Jun 2006

#ifndef GENERIC_MAP_H_
#define GENERIC_MAP_H_

#include "Vect.h"
#include "MolFramework.h"
#include "EDStructure.h"


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Indices of a single point in the electron density grid.  Used internally
//   by the GenericMap and EDStructure classes.
////////////////////////////////////////////////////////////////////////////////
struct gridPoint {  // single point on grid
  int i;
  int j;
  int k;
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   An abstract electron density map class.  This contains the features that
//   should be common to all types of electron density maps, e.g. storing the
//   electron density data, converting between its crystallographic grid 
//   coordinates and 3D Cartesian coordinates, and calculating its correlation
//   with an EDStructure object.
////////////////////////////////////////////////////////////////////////////////
class GenericMap {
protected:
  double a; // unit cell parameters
  double b;
  double c;
  double alpha;
  double beta;
  double gamma;
  int dataSize;      // equal to extent.i*extent.j*extent.k
  Vector i, j, k;    // unit cell vectors
  Vector offset;     // real-space offset from origin
  gridPoint origin;  // grid origin
  double *gridData;  // array containing ED data
  gridPoint extent;   // grid extent
public:
  GenericMap();
  ~GenericMap();
  gridPoint getExtent() {return extent;}
  Vector gridPointToVector(gridPoint gp);  // convert i,j,k to x,y,z
  double correlate(EDStructure &pdb);
  void writeEZD(string fileName); // output EZD map
  void trimMap(int scalingFactor); 
protected:
  gridPoint gridIndex(int n);
   // converts a position in the 1D data array to a position on the 3D grid
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Derived class of GenericMap -- this is designed to hold data taken from an
//   electron density file in the EZD (E-Z Density) format.  It inherits the
//   general functionalities from GenericMap and adds a few that are specific
//   to this data format.
////////////////////////////////////////////////////////////////////////////////
class EZDMap : public GenericMap {
private:
  char comment[80];
  double scale;        // scaling factor
  gridPoint num;       // number of grid points
public:
  EZDMap(string fileName);   // constructor for EZD_map
  ~EZDMap() {             // destructor    
    delete [] gridData; }             
  void displayHeader();    // display file header info
private:
  void readEZDData(string fileName, double G[]);
    // reads EZD data into grid
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   "Theoretical" electron density map derived from a PDB structure.  The 
//   only addition to GenericMap in this case is the constructor.
////////////////////////////////////////////////////////////////////////////////
class TheoMap : public GenericMap {

public:
  TheoMap(string fileName, double resolutionFactorInput = 0.0, 
          double gridDensity = 1.0, double noiseLevel = 0.0);  
  // constructor; forms map from PDB file
  ~TheoMap();
};

#endif
