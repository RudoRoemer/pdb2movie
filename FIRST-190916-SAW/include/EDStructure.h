// Author name:  Craig Jolley
// Created:      19 Jun 2006

#ifndef ED_STRUCTURE_H_
#define ED_STRUCTURE_H_

#include "MolFramework.h"
#include "Grid.h"

class Froda;
  // this gets defined more completely in Froda.h
class GenericMap;
  // this gets defined more completely in GenericMap.h

////////////////////////////////////////////////////////////////////////////////
// Description:
//   EDStructure objects are used by the routines for Monte-Carlo fitting of 
//   FRODA conformers into electron density maps.  The heart of EDStructure is
//   its vector of pointers to Atom objects; the information contained in these
//   Atom objects can be derived either from a PDB file or from a running 
//   FRODA simulation, depending on which constructor is used.
////////////////////////////////////////////////////////////////////////////////
class EDStructure {
private:
  //////////////////////////////////////////////////////////////////////////////
  // Description:
  //   Private nested class in EDStructure; contains the electron density 
  //   lookup tables.
  //////////////////////////////////////////////////////////////////////////////
  class LookupTable {
  public:
    double ai[7][4];  // first index is element, second is 1..4;
    double bi[7][4];
    double allZero[4];
    double binSize;    // use same bin size for all lookup tables
    double * cDensity;  // electron density lookup table for carbon
    double * nDensity;
    double * oDensity;
    double * sDensity;
    double * pDensity;
    double * feDensity; // iron 2+ ion
    double * mgDensity; // magnesium 2+ ion
    LookupTable(double resolutionFactor, double &cutoff);
    ~LookupTable();
  };
  ////////////////////////////////////////////////////////////////////////////////
  // Description:
  //   Pirvate nested class in EDStructure; holds the information for an 
  //   individual atom.  An Atom object knows its 
  //   location and identity and is able to calculate its distance from a 
  //   given point and the electron density it contributes at a given distance.
  //   Atom objects are used internally by EDStructure.
  ////////////////////////////////////////////////////////////////////////////////
  class Atom {
  public:
    double *x;  // atom coordinates
    double *y;
    double *z;
    double *lookup; // points to a lookup array
    double *ai;     // points to an array of ai (see LookupTable)
    double *bi;     // points to an array of bi
  public:
    Atom(){};
    virtual ~Atom(){};
  public:
    SiteID atomID;
    double distanceSquared(double x2, double y2, double z2);
     // calculates the square of the distance from the atom to a point
    virtual double density(double r2, double binSize);
  };
  // Begin EDStructure data & methods
public:
  Vector max, min;  // maximum and minimum coordinates
  MolFramework *structure;
    // points to the MolFramework object for the molecule being simulated
  double cachedCorrelation;
    // make sure this is initialized before MC starts
  Grid *coarseGrid; // coarse grid for faster evaluation of densityAt()
private:
  enum elem {carbon, nitrogen, oxygen, sulfur, phosphorus, iron, magnesium};
  std::vector<Atom*> atoms; // array of pointers to atoms
  double resolutionFactor;
  double cutoff;  // cutoff for calculation of electron density
  LookupTable *lookupTable;
public:
  EDStructure(string fileName, double resolutionFactorInput);
    // constructor initializes data in class by reading in a PDB file
  EDStructure(Froda &froda, MolFramework &structureInput, GenericMap &gmap, double resolutionFactorInput);
    // constructor initializes pointers in atoms vector by referencing them to 
    // the appropriate data in FIRST/FRODA
  ~EDStructure();
  double densityAt(Vector v);
    // returns the calculated electron density at x,y,z
};

#endif
