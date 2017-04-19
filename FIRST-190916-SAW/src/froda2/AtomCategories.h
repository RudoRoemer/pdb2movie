#ifndef ATOMCATEGORIES_H_
#define ATOMCATEGORIES_H_

#include <vector>
#include <string>
#include <map>
#include "Repulsion.h"
class AmberPrmtop;
class PDB;
class NeighborTable;

class AtomCategories : public Repulsion
{
public:
	AtomCategories( const NeighborTable *neighborTable_, const AmberPrmtop *prmtop );
	AtomCategories( const NeighborTable *neighborTable_, const PDB *pdb );
	virtual ~AtomCategories();

	int getPairType( int atomindex1, int atomindex2 ) const;
	const std::string& getPairName( int pairID ) const { return pairName[pairID]; }
	int getNPairTypes() const { return pairName.size(); }
	double getInteractionCutoff( int p1, int p2 ) const {
	  int pairID = getPairType( p1, p2 );
	  return ( pairID == -1 ) ? 0.0 : getInteractionCutoffForPairType( pairID );
	}
	double getInteractionCutoffForPairType( int pairID ) const {
	  return lookupCutoffFromPairID[pairID];
	}
	double getMaxInteractionCutoff( int p ) const {
	  return lookupMaximumCutoffForBasicType[atomIndexToBasicType[p]];
	}
	
	
private:
  const NeighborTable *neighborTable;
  std::vector<std::string> basicTypes;
  std::map<std::string, int> mapBasicTypeNameToTypeIndex;
  std::map<std::string, int> mapPairNameToPairID;
  std::vector<int> basicPairTypeLookupStorage;
  std::vector<int*> basicPairTypeLookup;
  std::vector<std::string> pairName;
  int pairID_OH_O;
  int pairID_NH_O;
  int basicType_oxygen;
  int basicType_nitrogen;
  std::vector<char> potentialDonors;
  std::vector<int> atomIndexToBasicType;
  std::vector<double> lookupCutoffFromPairID;
  std::vector<double> lookupMaximumCutoffForBasicType;
  
  void defineTypes();
  void assignBasicTypesToAtoms( const AmberPrmtop *prmtop );
  void assignBasicTypesToAtoms( const PDB *pdb );
};

#endif /*ATOMCATEGORIES_H_*/
