#ifndef AMBERPRMTOP_H_
#define AMBERPRMTOP_H_

#include <string>
#include <vector>
#include <set>

class AmberPrmtop {
public:
  int natom;
  int ntypes;
  int nexc;
  int nres;
  int ifbox;
  std::vector<int> atomTypeIndex;
  std::vector<int> numExcludedAtoms;
  std::vector<int> nonbondedParmIndexStorage;
  std::vector<int*> nonbondedParmIndex;
  std::vector<std::string> residueLabels;
  std::vector<int> residueFirstAtomIndex;
  std::vector<int> lookupResidueIndexFromAtomIndex;
  std::vector<double> lennardJonesACoef;
  std::vector<double> lennardJonesBCoef;
  std::vector<double> bondEquilValue;
  std::vector<int> excludedAtomsList;
  std::vector< std::set<int> > excludedAtoms;
  std::vector<double> charge;
  std::vector< std::vector<std::pair<int,int> > > neighborTable;
  std::vector<std::string> typeIDtoTypeString;
  std::vector<std::string> amberAtomType;
  std::vector<std::string> atomName;
  AmberPrmtop( std::string filename );
};

#endif
