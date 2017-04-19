#ifndef SAXS_H_
#define SAXS_H_

#include "MolFramework.h"

class Froda;
// this gets defined more completely in Froda.h

////////////////////////////////////////////////////////////////////////////////
// Description:
//   SAXS objects are used by the routines for Monte-Carlo fitting of 
//   FRODA conformers to SAXS profiles.  The heart of a SAXS object is
//   its vector of pointers to Atom objects; the information contained in these
//   Atom objects can be derived either from a PDB file or from a running 
//   FRODA simulation, depending on which constructor is used.
////////////////////////////////////////////////////////////////////////////////
class SAXS {
private:
  //////////////////////////////////////////////////////////////////////////////
  // Description:
  //   Private nested class in SAXS; contains various lookup tables.
  //////////////////////////////////////////////////////////////////////////////
  class LookupTable {
  public:
    double ai[7][4];  // first index is element, second is 1..4;
    double bi[7][4];
    double allZero[4];
    double binSize;    // use same bin size for all lookup tables
    double * sinXOverX;
    double * eToTheMinusX;
    LookupTable(double resolutionFactor);
    ~LookupTable();
  };
  ////////////////////////////////////////////////////////////////////////////////
  // Description:
  //   Private nested class in SAXS; holds the information for an 
  //   individual atom.  An Atom object knows its 
  //   location and identity and is able to calculate its distance from a 
  //   given point and the electron density it contributes at a given distance.
  ////////////////////////////////////////////////////////////////////////////////
  class Atom {
  public:
    double *x;  // atom coordinates
    double *y;
    double *z;
    double *ai;     // points to an array of ai (see LookupTable)
    double *bi;     // points to an array of bi
  public:
    double distanceSquared(const Atom &a2);
     // calculates the square of the distance between two atoms
  };
  ////////////////////////////////////////////////////////////////////////////////
  // Description:
  //   Private nested class in SAXS; a vector of these forms a SAXS spectrum
  ////////////////////////////////////////////////////////////////////////////////
  struct SAXSPoint {
    double k;
    double intensity;
  };
  // Begin SAXS data & methods
public:
  MolFramework *structure;
    // points to the MolFramework object for the molecule being simulated
  double cachedCorrelation;
    // make sure this is initialized before MC starts
private:
  enum elem {carbon, nitrogen, oxygen, sulfur, phosphorus, iron, magnesium};
  std::vector<Atom*> atoms; // array of pointers to atoms
  double resolutionFactor;
  LookupTable *lookupTable;
  std::vector<SAXSPoint> saxsTarget;
  std::vector<SAXSPoint> saxsCurrent;
public:
  SAXS(Froda &froda, MolFramework &structureInput, double resolutionFactorInput, 
    string filename = "");
    // constructor initializes pointers in atoms vector by referencing them to 
    // the appropriate data in FIRST/FRODA
    // TODO: This should be able to read in a file containing SAXS data
  ~SAXS();
  double saxs(double k);
  // returns the SAXS scattering intensity at wavenumber k
  double restrictedSAXS(double k, int degree); 
  // returns only the contribution to the scattering of nth nearest neighbors
  void loadTarget(string filename);
  // loads a SAXS spectrum from file into saxsTarget
  void writeSAXS(string filename);
  // output calculated SAXS profile to filename
  void updateProfile();
  // update saxsCurrent based on current position of molecule in FIRST
  double correlate();
  // calculates correlation between saxsTarget and saxsCurrent
};

#endif
