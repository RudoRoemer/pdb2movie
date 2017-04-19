#include "AtomCategories.h"
#include "AmberPrmtop.h"
#include "PDB.h"
#include "NeighborTable.h"

using namespace std;

AtomCategories::AtomCategories( const NeighborTable *neighborTable_, const AmberPrmtop *prmtop ) :
  neighborTable( neighborTable_ )
{
  defineTypes();
  assignBasicTypesToAtoms( prmtop ); 
}

AtomCategories::AtomCategories( const NeighborTable *neighborTable_, const PDB *pdb ) :
  neighborTable( neighborTable_ )
{
  defineTypes();
  assignBasicTypesToAtoms( pdb ); 
}

AtomCategories::~AtomCategories()
{
}

void AtomCategories::defineTypes() {
  //establish basic atom types
  basicTypes.push_back("C=O");
  basicTypes.push_back("C");
  basicTypes.push_back("N");
  basicTypes.push_back("O");
  basicTypes.push_back("S");
  basicTypes.push_back("Hpolar");
  basicTypes.push_back("Hnonpolar");
  for ( size_t i = 0; i < basicTypes.size(); i++ ) {
    mapBasicTypeNameToTypeIndex[ basicTypes[i] ] = i;
  }
  basicType_oxygen = mapBasicTypeNameToTypeIndex["O"];
  basicType_nitrogen = mapBasicTypeNameToTypeIndex["N"];
  
  
  //establish pair types
  basicPairTypeLookupStorage.resize( basicTypes.size()*basicTypes.size() );
  basicPairTypeLookup.resize( basicTypes.size() );
  for ( size_t i = 0; i < basicTypes.size(); i++ ) {
    basicPairTypeLookup[i] = &basicPairTypeLookupStorage[i*basicTypes.size()];
  }

  int pairID = 0;
  for ( size_t i = 0; i < basicTypes.size(); i++ ) {
    for ( size_t j = i; j < basicTypes.size(); j++ ) {
      pairName.push_back( basicTypes[i] + "_" + basicTypes[j] );
      basicPairTypeLookup[i][j] = pairID;
      basicPairTypeLookup[j][i] = pairID;
      pairID++;
    }
  }

  pairName.push_back( "OH_O" );
  pairID_OH_O = pairID++;
  pairName.push_back( "NH_O" );
  pairID_NH_O = pairID++;

  for ( size_t i = 0; i < pairName.size(); i++ ) {
    mapPairNameToPairID[ pairName[i] ] = i;
  }
  
  lookupCutoffFromPairID.resize( pairName.size() );
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_C=O"]] = 3.11;
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_C"]] = 3.09;
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_N"]] = 3.02;
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_O"]] = 2.86;
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_S"]] = 3.30;
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_Hpolar"]] = 2.26;
  lookupCutoffFromPairID[mapPairNameToPairID["C=O_Hnonpolar"]] = 2.54;
  lookupCutoffFromPairID[mapPairNameToPairID["C_C"]] = 3.38;
  lookupCutoffFromPairID[mapPairNameToPairID["C_N"]] = 3.07;
  lookupCutoffFromPairID[mapPairNameToPairID["C_O"]] = 3.04;
  lookupCutoffFromPairID[mapPairNameToPairID["C_S"]] = 3.38;
  lookupCutoffFromPairID[mapPairNameToPairID["C_Hpolar"]] = 2.43;
  lookupCutoffFromPairID[mapPairNameToPairID["C_Hnonpolar"]] = 2.50;
  lookupCutoffFromPairID[mapPairNameToPairID["N_N"]] = 3.20;
  lookupCutoffFromPairID[mapPairNameToPairID["N_O"]] = 3.25;
  lookupCutoffFromPairID[mapPairNameToPairID["N_S"]] = 3.20;
  lookupCutoffFromPairID[mapPairNameToPairID["N_Hpolar"]] = 2.23;
  lookupCutoffFromPairID[mapPairNameToPairID["N_Hnonpolar"]] = 2.49;
  lookupCutoffFromPairID[mapPairNameToPairID["O_O"]] = 2.96;
  lookupCutoffFromPairID[mapPairNameToPairID["O_S"]] = 3.24;
  lookupCutoffFromPairID[mapPairNameToPairID["O_Hpolar"]] = 1.70;
  lookupCutoffFromPairID[mapPairNameToPairID["O_Hnonpolar"]] = 2.26;
  lookupCutoffFromPairID[mapPairNameToPairID["S_S"]] = 4.00;
  lookupCutoffFromPairID[mapPairNameToPairID["S_Hpolar"]] = 2.50;
  lookupCutoffFromPairID[mapPairNameToPairID["S_Hnonpolar"]] = 2.73;
  lookupCutoffFromPairID[mapPairNameToPairID["Hpolar_Hpolar"]] = 1.87;
  lookupCutoffFromPairID[mapPairNameToPairID["Hpolar_Hnonpolar"]] = 1.88;
  lookupCutoffFromPairID[mapPairNameToPairID["Hnonpolar_Hnonpolar"]] = 2.11;
  lookupCutoffFromPairID[mapPairNameToPairID["OH_O"]] = 2.63;
  lookupCutoffFromPairID[mapPairNameToPairID["NH_O"]] = 2.71;
    
  lookupMaximumCutoffForBasicType.resize( basicTypes.size(), 0 );
  for ( size_t t1 = 0; t1 < basicTypes.size(); t1++ ) {
    double cutoff;
    for ( size_t t2 = 0; t2 < basicTypes.size(); t2++ ) {
      cutoff = lookupCutoffFromPairID[basicPairTypeLookup[t1][t2]];
      if ( cutoff > lookupMaximumCutoffForBasicType[t1] )
        lookupMaximumCutoffForBasicType[t1] = cutoff;
    }
  }
  
}

void AtomCategories::assignBasicTypesToAtoms( const AmberPrmtop *prmtop ) {
  //potential donors
  size_t natoms = prmtop->natom;
  potentialDonors.resize( natoms, 0 );
  for ( size_t i = 0; i < natoms; i++ ) {
    if ( prmtop->amberAtomType[i][0] != 'O' && prmtop->amberAtomType[i][0] != 'N' ) continue;
    
    //cycle through the neighbors, searching for a hydrogen
    const vector<int> *neighborlist = &(*neighborTable)[i];
    for ( size_t j = 0; j < neighborlist->size(); j++ ) {
      int atomindex = (*neighborlist)[j];
      if ( prmtop->amberAtomType[atomindex][0] == 'H') {
        potentialDonors[i] = 1; //1 means true
        break;
      }
    }
    
  }
  
  //assign basic atom types
  atomIndexToBasicType.resize( natoms );
  for ( size_t atom = 0; atom < natoms; atom++ ) {
    string ambertype = prmtop->amberAtomType[atom];
    if ( ambertype[0] == 'C' ) {
      if ( ambertype == "C" ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["C=O"];
      else atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["C"];
    }
    else if ( ambertype[0] == 'H' ) {
      if ( ambertype == "H" || ambertype == "HW" || ambertype == "HO" || ambertype == "HS" ) {
        atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["Hpolar"];
      }
      else atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["Hnonpolar"];
    }
    else if ( ambertype[0] == 'N' ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["N"];
    else if ( ambertype[0] == 'O' ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["O"];
    else if ( ambertype[0] == 'S' ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["S"];
    else {
      atomIndexToBasicType[atom] = -1;
    }
  }  
  
}

void AtomCategories::assignBasicTypesToAtoms( const PDB *pdb ) {
  cout << " Atom Categories temporarily not allowing PDB input.  Fix this." << endl;
  exit(0);
  //assign basic atom types and identify potential hydroden bond donor atoms
  size_t natoms = pdb->atomLines.size();
  potentialDonors.resize( natoms, 0 );
  atomIndexToBasicType.resize( natoms );
  for ( size_t atom = 0; atom < natoms; atom++ ) {
    string elem = pdb->atomLines[atom].element;
    if ( elem == "C" ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["C"];
    
    else if ( elem == "H" ) {
      //first, set the type to the default Hnonpolar.  Then check to see if
      //it is actually polar, and if so, change type to Hpolar.
      atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["Hnonpolar"];
      if ( (*neighborTable)[atom].size() == 1 ) {
        int donor = (*neighborTable)[atom][0];
        if ( pdb->atomLines[donor].element == "O" || pdb->atomLines[donor].element == "N" ) {
          atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["Hpolar"];
          potentialDonors[donor] = 1; //1 means true
        }
      }
    }
    
    else if ( elem == "N" ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["N"];
    else if ( elem == "O" ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["O"];
    else if ( elem == "S" ) atomIndexToBasicType[atom] = mapBasicTypeNameToTypeIndex["S"];
    else {
      atomIndexToBasicType[atom] = -1;
    }
  }  
  
}

int AtomCategories::getPairType( int atomindex1, int atomindex2 ) const {
  if ( atomindex1 == atomindex2 ||
       neighborTable->isFirstNeighbor( atomindex1, atomindex2 ) ||
       neighborTable->isSecondNeighbor( atomindex1, atomindex2 ) ||
       neighborTable->isThirdNeighbor( atomindex1, atomindex2 ) ) {
    return -1;
  }
    
  int basicType1 = atomIndexToBasicType[ atomindex1 ];
  int basicType2 = atomIndexToBasicType[ atomindex2 ];
  
  if ( basicType1 == basicType_oxygen && basicType2 == basicType_oxygen &&
       ( potentialDonors[atomindex1] || potentialDonors[atomindex2] ) ) {
    return pairID_OH_O;
  }
  
  if ( (basicType1 == basicType_nitrogen && basicType2 == basicType_oxygen && potentialDonors[atomindex1] )
       || (basicType2 == basicType_nitrogen && basicType1 == basicType_oxygen && potentialDonors[atomindex2]) ){
    return pairID_NH_O;
  }

  if ( basicType1 == -1 || basicType2 == -1 ) return -1;
  return basicPairTypeLookup[basicType1][basicType2];
}
