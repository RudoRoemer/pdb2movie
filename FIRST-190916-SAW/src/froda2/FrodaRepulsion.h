#ifndef FRODAREPULSION_H_
#define FRODAREPULSION_H_

#include "Repulsion.h"
class PDB;
#include "NeighborTable.h"
#include <vector>
#include <map>
using namespace std;
class FrodaRepulsion : public Repulsion
{
public:
	FrodaRepulsion( const PDB &pdb,
                  const NeighborTable &neighborTable_,
                  double mismatchTol );
	virtual ~FrodaRepulsion();
  virtual double getInteractionCutoff( int p1, int p2 ) const;
  double getMaxInteractionCutoff( int p ) const;
  void printinfo() const;
private:
  bool polarGeometry;
  double polarHRadius;
  double vdwTol;
  double mismatchTol;
  const NeighborTable &neighborTable;
  vector<string> nameOfType;
  vector<double> radiusOfType;
  vector<signed char> chargeOfType;
  vector<bool> isPotentialDonorType;
  int nTypes;
  int nPairs;
  map<string, int> mapNameToType;
  vector<int> atomTypeLookup;
  
  vector<double*> rInnerLookup;
  vector<double*> rInner_ForThirdNeighbors_Lookup;
  vector<double> rInnerLookupStorage;
  vector<double> rInner_ForThirdNeighbors_LookupStorage;
  vector<bool> isHydrogen;
  void setupTypes();
  void setupInteractionCutoffLookupTable();
  void assignTypes( const PDB &pdb );
  double calcInteractionCutoffForTypePair( int t1, int t2 ) const;
};

#endif /*FRODAREPULSION_H_*/
