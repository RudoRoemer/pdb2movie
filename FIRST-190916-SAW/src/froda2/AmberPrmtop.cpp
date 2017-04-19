#include "AmberPrmtop.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>


using namespace std;

template <class T>
class AmberPrmtopDataSection {
public:
  std::string sectionstring;
  std::string formatstring;
  int numFieldsPerLine;
  char dataTypeCode;
  int fieldWidth;
  std::vector<T> data;
  
  AmberPrmtopDataSection();
  AmberPrmtopDataSection( std::ifstream& prmtop, const std::string& sectionstring );
  AmberPrmtopDataSection( std::ifstream& prmtop, const char* sectionstring );
  
  ~AmberPrmtopDataSection();
  void getSection( std::ifstream& prmtop, const char* sectionstring );
  void getSection( std::ifstream& prmtop, const std::string& sectionstring );
  void clear();
  std::string str();

  bool checkDataTypeCode( char code );
};

template <class T> inline bool AmberPrmtopDataSection<T>::checkDataTypeCode( char code ) { return false; }
template <> inline bool AmberPrmtopDataSection<std::string>::checkDataTypeCode( char code ) { return toupper(code)=='A'; }
template <> inline bool AmberPrmtopDataSection<double>::checkDataTypeCode( char code ) { return toupper(code)=='E'; }
template <> inline bool AmberPrmtopDataSection<int>::checkDataTypeCode( char code ) { return toupper(code)=='I'; }

template <class T>
AmberPrmtopDataSection<T>::AmberPrmtopDataSection() {
}

template <class T>
AmberPrmtopDataSection<T>::AmberPrmtopDataSection( std::ifstream& prmtop, const std::string& sectionstring ) {
  getSection( prmtop, sectionstring );
}

template <class T>
AmberPrmtopDataSection<T>::AmberPrmtopDataSection( std::ifstream& prmtop, const char* sectionstring ) {
  getSection( prmtop, sectionstring );
}

template <class T>
void AmberPrmtopDataSection<T>::clear() {
  sectionstring = "EMPTY";
  formatstring = "";
  numFieldsPerLine=0;
  dataTypeCode = 0;
  fieldWidth=0;
  data.clear();
}

template <class T>
void AmberPrmtopDataSection<T>::getSection( std::ifstream& prmtop, const char* sectionstring ) {
  getSection( prmtop, std::string(sectionstring) );
}

template <class T>
void AmberPrmtopDataSection<T>::getSection( std::ifstream& prmtop, const std::string& section ) {
  clear();
  std::string currentline="";
  std::string searchsectionstring = "%FLAG " + section;
  int initialfileposition = (int)prmtop.tellg();
  // keep searching forward in file until we find the
  // start of a new data section ( or until we hit EOF )
  bool EOFseenAlready = false;
  while ( true ) {
    // if we hit EOF, begin searching from the beginning of the file
    if ( EOFseenAlready && (int)prmtop.tellg() >= initialfileposition ) {
      std::cout << "Error: can't find Amber7Prmtop section " << section;
      exit(1);
    }
    if ( prmtop.eof() ) {
      prmtop.clear();
      EOFseenAlready = true;
      prmtop.seekg( 0 );
    }
    getline( prmtop, currentline );
    if ( currentline.find( "%FLAG", 0 ) != std::string::npos ) {
      std::stringstream ss;
      std::string token;
      ss << currentline.substr(5, std::string::npos );
      ss >> token;
      if ( token == section ) break;
    }
  }
    
  size_t endOfSectionString = currentline.find_first_of( ' ', 6 );
  sectionstring = currentline.substr( 6, endOfSectionString - 6);
    
  // the next line in the file should be a %FORMAT... line
  getline( prmtop, currentline );
  if ( currentline.find( "%FORMAT", 0, 0 ) == std::string::npos ) {
    std::cout << "Error in prmtop: FORMAT line expected" << std::endl;
    exit(0);
  }
  
  size_t endOfFormatString = currentline.find_first_of( ')', 7 );
  formatstring = currentline.substr( 8, endOfFormatString - 8);
  int i1=0;
  while ( isdigit( formatstring[i1] ) ) {i1++;}
  numFieldsPerLine = atoi( formatstring.substr( 0, i1 ).c_str() );
  dataTypeCode = formatstring[i1++];
  
  if ( !checkDataTypeCode(dataTypeCode) ) {
    std::cout << "Error in prmtop: data type" << std::endl;
    exit(0);
  }
  
  int i2=i1;
  while ( isdigit( formatstring[i2] ) ) {i2++;}
  fieldWidth = atoi( formatstring.substr( i1, i2-i1 ).c_str() );
  
  // the next lines in the file contain the data
  T t;
  while ( !prmtop.eof() && prmtop.peek() != '%' ) { //loop over lines
    getline( prmtop, currentline );
    int offset=0;
    while ( true ) { // loop over fields in the line
      if ( offset + fieldWidth > static_cast<int> (currentline.size()) ) break;
      std::string fieldString = currentline.substr( offset, fieldWidth );
      if ( fieldString.find_first_not_of(" \n\t") == std::string::npos ) break;
      std::stringstream ss;
      ss << fieldString;
      // store a single field
      ss >> t;
      data.push_back( t );
      offset+=fieldWidth;
    } 
  }
}

template <class T>
AmberPrmtopDataSection<T>::~AmberPrmtopDataSection() {}

template <class T>
std::string AmberPrmtopDataSection<T>::str() {
  std::stringstream ss;
  ss << "Found prmtop section: " << sectionstring << std::endl;
  ss << "  Format string: " << formatstring << std::endl;
  ss << numFieldsPerLine << " " << dataTypeCode << " " << fieldWidth << std::endl;
  for ( size_t i=0; i<data.size(); i++ ) {
    ss << data[i] << std::endl;
  }
  return ss.str();
}

AmberPrmtop::AmberPrmtop( string filename ) {
  ifstream prmtop( filename.c_str(), ios::in );
  if (!prmtop) {
    cout << "Could not open Amber7 prmtop file: " << filename << endl;
    exit(1);
  }
  AmberPrmtopDataSection<int> pointersSection( prmtop, "POINTERS" );
  AmberPrmtopDataSection<string> atomNameSection( prmtop, "ATOM_NAME" );
  AmberPrmtopDataSection<double> chargeSection( prmtop, "CHARGE" );
  AmberPrmtopDataSection<int> atomTypeIndexSection( prmtop, "ATOM_TYPE_INDEX" );
  AmberPrmtopDataSection<int> numberExcludedAtomsSection( prmtop, "NUMBER_EXCLUDED_ATOMS" );
  AmberPrmtopDataSection<int> nonbondedParmIndexSection( prmtop, "NONBONDED_PARM_INDEX" );
  AmberPrmtopDataSection<string> residueLabelSection( prmtop, "RESIDUE_LABEL" );  
  AmberPrmtopDataSection<int> residuePointerSection( prmtop, "RESIDUE_POINTER" );  
  AmberPrmtopDataSection<double> bondEquilValueSection( prmtop, "BOND_EQUIL_VALUE" );
  AmberPrmtopDataSection<double> lennardJonesACoefSection( prmtop, "LENNARD_JONES_ACOEF" );
  AmberPrmtopDataSection<double> lennardJonesBCoefSection( prmtop, "LENNARD_JONES_BCOEF" );
  AmberPrmtopDataSection<int> bondsIncHydrogenSection( prmtop, "BONDS_INC_HYDROGEN" );
  AmberPrmtopDataSection<int> bondsWithoutHydrogenSection( prmtop, "BONDS_WITHOUT_HYDROGEN" );
  AmberPrmtopDataSection<int> excludedAtomsListSection( prmtop, "EXCLUDED_ATOMS_LIST" );
  AmberPrmtopDataSection<string> amberAtomTypeSection( prmtop, "AMBER_ATOM_TYPE" );

  // POINTERS section contains a list of various quanitites of interest
  // Here we pick out the ones we want.
  natom = pointersSection.data[0];
  ntypes = pointersSection.data[1]; //number of lennard-jones atom types
  nexc = pointersSection.data[10];
  nres = pointersSection.data[11];
  ifbox = pointersSection.data[27];
  
  //section CHARGE.  List of the charge on each atom.
  charge = chargeSection.data;

  //section ATOM_TYPE_INDEX contains an integer Lennard-Jones atom type index from 1..ntypes
  //(fortran-style indexing).  Here we subtract one from the index, to make it
  //c-style.
  atomTypeIndex = atomTypeIndexSection.data;
  for ( size_t i=0; i<atomTypeIndex.size(); i++ ) {
    atomTypeIndex[i]--;
  }
  
  //section NUMBER_EXCLUDED_ATOMS: for each atom, the number of atoms that are
  //excluded from non-bonded interactions
  numExcludedAtoms = numberExcludedAtomsSection.data;
  
  //section NONBONDED_PARM_INDEX: this list of indices is to be interpreted as
  //a 2-D look-up table.  For two lennard-jones atom types, you can look up the corresponding
  //index into the nonbonded parameter lists.  This lookup table is size
  //ntypes * ntypes.
  //nonbondedParmIndexStorage here gets the 1-D version of the lookup table
  nonbondedParmIndexStorage = nonbondedParmIndexSection.data;

  //nonbondedParmIndex here is set up as the 2-D lookup table.
  //It is to be used like this:  nonbondedParmIndex[type1][type2],
  //where type1 and type2 are an integer between 0 and ntypes-1.
  nonbondedParmIndex.resize(ntypes);
  for ( int i=0; i<ntypes; i++ ) {
    nonbondedParmIndex[i] = &nonbondedParmIndexStorage[i*ntypes];
  }
  //Now we must adjust the fortran-like indices stored in the lookup table
  //to c-style
  for ( size_t i=0; i<nonbondedParmIndexStorage.size(); i++ ) {
    if (nonbondedParmIndexStorage[i] > 0) nonbondedParmIndexStorage[i]--;
  }
  
  //section RESIDUE_LABEL lists the residue string for each residue. 
  residueLabels = residueLabelSection.data;

  //section RESIDUE_POINTER lists the first atom index of each residue.  Here we change
  //the data to C-Style indices.
  residueFirstAtomIndex = residuePointerSection.data;
  for ( int residueIndex = 0; residueIndex < nres; residueIndex++) {
    residueFirstAtomIndex[residueIndex]--;
  }
  
  //Here we build lookupResidueIndexFromAtomIndex, which holds the residue index
  //of each atom.
  lookupResidueIndexFromAtomIndex.resize( natom );
  int residueIndex = 0;
  int nextResiStart = (residueIndex == nres-1) ? natom : residueFirstAtomIndex[residueIndex+1];
  for ( int atomIndex = 0; atomIndex < natom; atomIndex++ ) {
    if ( atomIndex == nextResiStart ) {
      residueIndex++;
      nextResiStart = (residueIndex == nres-1) ? natom : residueFirstAtomIndex[residueIndex+1];
    }
    lookupResidueIndexFromAtomIndex[atomIndex] = residueIndex;
  }
  
  //section BOND_EQUIL_VALUE lists the equilibrium bond lengths,
  //one for each unique bond type.
  bondEquilValue= bondEquilValueSection.data;


  //section LENNARD_JONES_ACOEF lists the coefficient A, one for each
  //nonbonded interaction id.  NOTE this is not listed per atom.  
  //From two atom types, you look up the nonbonded parm index, and use
  //this index to get the corresponding A coefficient.
  lennardJonesACoef= lennardJonesACoefSection.data;
  

  //section LENNARD_JONES_BCOEF, like the A coef above.
  lennardJonesBCoef= lennardJonesBCoefSection.data;
  
  //section BONDS_INC_HYDROGEN.  Data is in triples: 1st integer is atom1,
  //2nd integer is atom2, 3rd integer is the bond data index.  Need to
  //divide atom integers by 3 (and don't add +1) to get the c-Style
  //atom number.  Need to subtract 1 from the bond data index to 
  //convert it to c-style.
  vector<int> bondsIncHydrogen= bondsIncHydrogenSection.data;
  for ( size_t i=0; i<bondsIncHydrogen.size(); i+=3 ) {
    bondsIncHydrogen[i] /= 3;
    bondsIncHydrogen[i+1] /= 3;
    bondsIncHydrogen[i+2] -= 1;
  }
  neighborTable.resize(natom);
  //store bonding data in neighbor table
  for ( size_t i=0; i<bondsIncHydrogen.size(); i+=3 ) {
    int atom1 = bondsIncHydrogen[i];
    int atom2 = bondsIncHydrogen[i+1];
    int bondindex = bondsIncHydrogen[i+2];
    neighborTable[atom1].push_back(pair<int,int>(atom2, bondindex));
    neighborTable[atom2].push_back(pair<int,int>(atom1, bondindex));
  }
    
  //section BONDS_WITHOUT_HYDROGEN.  Data is in triples: 1st integer is atom1,
  //2nd integer is atom2, 3rd integer is the bond data index.  Need to
  //divide atom integers by 3 (and don't add +1) to get the c-Style
  //atom number.  Need to subtract 1 from the bond data index to 
  //convert it to c-style.
  vector<int> bondsWithoutHydrogen= bondsWithoutHydrogenSection.data;
  for ( size_t i=0; i<bondsWithoutHydrogen.size(); i+=3 ) {
    bondsWithoutHydrogen[i] /= 3;
    bondsWithoutHydrogen[i+1] /= 3;
    bondsWithoutHydrogen[i+2] -= 1;
  }
  //store bonding data in neighbor table
  for ( size_t i=0; i<bondsWithoutHydrogen.size(); i+=3 ) {
    int atom1 = bondsWithoutHydrogen[i];
    int atom2 = bondsWithoutHydrogen[i+1];
    int bondindex = bondsWithoutHydrogen[i+2];
    neighborTable[atom1].push_back(pair<int,int>(atom2, bondindex));
    neighborTable[atom2].push_back(pair<int,int>(atom1, bondindex));
  }


  //section EXCLUDED_ATOMS_LIST.  This list is a bit tricky.  For each atom,
  //there is a list of the ids of its excluded atoms.  Each atom's list of
  //excluded atoms is concatenated into one big long list.  To know when
  //one atom's list ends and the next begins, we need the info stored
  //in the NUMBER_EXCLUDED_ATOMS section.
  excludedAtomsList= excludedAtomsListSection.data;
  // first, fix fortran-like indices stored in this array
  for ( size_t i=0; i<excludedAtomsList.size(); i++ ) {
    excludedAtomsList[i]--;
  }
  
  //now, using the number of excluded atoms, we build up the excluded atoms lists
  //for each atom.
  excludedAtoms.resize( natom );
  int istart;
  int iend=0;
  for ( int i=0; i < natom; i++ ) {
    istart = iend;
    iend+=numExcludedAtoms[i];
    for ( int j=istart; j<iend; j++ ) {
      // prmtop file lists a 0 in the EXCLUDED ATOMS LIST to 
      // signify that no atoms are to be excluded.  This can be confusing,
      // because you would think that it meant to exclude atom 0.
      // Here, we check for this special case (however, we have 
      // already subtracted 1 from each index in the list, and so
      // instead of looking for 0, we look for -1).
      if ( excludedAtomsList[j] < 0 ) {
        excludedAtoms[i].insert( excludedAtomsList[j] );
      }
    }
  }

  //section AMBER_ATOM_TYPE.  List of each atom's type (a string).
  //note - this is different from the atom's lennard-jones type.
  amberAtomType = amberAtomTypeSection.data;

  atomName = atomNameSection.data;
  prmtop.close();
}
