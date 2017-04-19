////////////////////////////////////////////////////////////////////////////////
// CovalentBonder.h
//
// This header file defines the CovalentBonder class, as well as some helper
// classes.
////////////////////////////////////////////////////////////////////////////////

#ifndef COVALENTBONDER_H_
#define COVALENTBONDER_H_

#include "MolFramework.h"
#include "Site_Info.h"
#include "global_defs.h"

using namespace std;

typedef string atomname_t;
typedef string resname_t;
typedef string elemname_t;


// The class BondDictionary holds the PDB bonding templates.
// It has member functions for reading the standard dictionary
// and the het dictionary, and for looking up atoms in the dictionary.
// 
// In principle, one can lookup bonding information for an atom as follows.
// If my BondDictionary object is bonddict, I would write
//   bonddict.getRes("GLY")->getAtom("CA")->isNeighbor("CB")
// which would return a bool indicating whether the atom names are neighbors.
// In practice, you don't want to follow a pointer until you first check
// to see whether it is NULL.  NULL indicates that the thing you are looking
// up is not found in the dictionary.  So at each step (say, after getRes, or getAtom),
// you check the pointer before following it.
//////////////////////////////////////////////////////////////////////
class BondDictionary {

public:
  class Res;
  class Atom;
  BondDictionary();
  ~BondDictionary();
  Res *getRes(resname_t);
  void identify_het_residues(string filename, MolFramework *molFramework);
  void read_standard_residue_library(string filename);
private:
  typedef map<resname_t, Res> resmap_t;
  resmap_t resmap;
  bool search_het_group_dictionary( string filename, resname_t residue_name );
  void store_residue_template_data( fstream *file_ptr, resname_t residue_name );
};

//////////////////////////////////////////////////////////////////////
class BondDictionary::Res {
private:
  friend class BondDictionary;
  typedef map<atomname_t, BondDictionary::Atom> atommap_t;
  atommap_t atommap;
public:
  BondDictionary::Atom *getAtom(atomname_t);
  Res();
  ~Res();
};

//////////////////////////////////////////////////////////////////////
class BondDictionary::Atom{
private:
  friend class BondDictionary::Res;
  friend class BondDictionary;
  typedef set<atomname_t> neighborset_t;
  neighborset_t neighborset;  
public:
  Atom();
  ~Atom();
  bool isNeighbor(atomname_t atomname) const;
  size_t numNeighbors() const;
};

inline BondDictionary::BondDictionary() {};
inline BondDictionary::~BondDictionary() {};
inline BondDictionary::Res::Res() {};
inline BondDictionary::Res::~Res() {};
inline BondDictionary::Atom::Atom() {};
inline BondDictionary::Atom::~Atom() {};


// The class BondDistanceTable holds the element-element bond distance cutoffs.
// It has member functions for reading the cutoffs from a file.  It also has
// const_iterators for read-only lookup into the table.
class BondDistanceTable {
private:
  typedef map<string, float> distmap_t;
  distmap_t distmap;
public:
  BondDistanceTable() {};
  BondDistanceTable( string filename );
  ~BondDistanceTable() {};
  void read_file( string filename );
  typedef distmap_t::const_iterator const_iterator;
  const_iterator find( elemname_t elem1, elemname_t elem2 ) const;
  const_iterator begin() const;
  const_iterator end() const;
};

// The CONECT_info class allows quick checking for the presence of
// a CONECT record.  It initializes from the CONECT data in a MolFramework
// object.  In MolFramework, the CONECT data is stored straight from the PDB
// file, so it uses the original atom numbers.  But to be useful, we need to
// search based on FIRST SiteIDs.  The initialization takes care of this
// translation from orig ID to SiteID, and internally stores the CONECT records
// indexed by SiteID.
////////////////////////////////////////////////////////////////////////////////
class CONECT_info {

private:
  vector< set<SiteID> > conect_record_neighborlist;

public:
  CONECT_info() {};
  ~CONECT_info() {};
  void get_from_molFramework( MolFramework* molFramework );
  bool is_listed_in_conect_records(SiteID atom1, SiteID atom2) const;
};

// The CovalentBonder class inspects the atoms in a MolFramework object,
// searching for and adding covalent bonds.  It uses the PDB Bond Dictionary
// files, the covalent bond distance cutoff file, and CONECT records
// to determine bonding.  It checks for certain irregularities as it bonds,
// and issues warnings as appropriate.
//
// To bond atoms, just create a CovalentBonder object, and call the
// autobond member function.  Other public functions are also available.
//
// NOTE: Invoking the constructor will cause the template files and cutoff files
// to be parsed, so you probably only want to instantiate this object once.
class CovalentBonder {
public:
  CovalentBonder(MolFramework *molFramework); //Note: Default constructor has been intentionally omitted
    //to force user to pass in a molFramework object
  ~CovalentBonder();
  void autobond(); //finds bonds and adds bonds to molFramework 
  bool is_locked( SiteID atom1, SiteID atom2) const;
  bool labelled_as_disulfide_bond( SiteID atom_1, SiteID atom_2 ) const;
  bool labelled_as_internucleic_bond( SiteID atom_1, SiteID atom_2 ) const;
  bool labelled_as_peptide_bond( SiteID atom1, SiteID atom2 ) const;
  bool meets_atom_atom_distance_cutoff( SiteID atom1, SiteID atom2 ) ;

private:
  MolFramework *molFramework;
  Site_Info *site_info;
  BondDictionary bonddict; //bonding connectivity from PDB dictionary files
  vector<BondDictionary::Atom*> dict_atom_lookup; //indexed like site_info, this vector
    //contains a pointer to each atom's corresponding entry in the BondDictionary object
  BondDistanceTable bondDistanceTable; //bond distances from bond distance cutoff file
  CONECT_info conect_info; //bonding connectivity from CONECT records
  
  // These variables are for storing any issues that arise during bonding
  set< pair<elemname_t, elemname_t> > missingBondDistances;
  set< pair<resname_t, atomname_t> > missingAtomNames;
  typedef list<pair<SiteID, SiteID > > ListOfSiteIDPairs;
  ListOfSiteIDPairs bondedIntraRes_KnownNeighborsButNotWithinDistance,
                    nonbondedIntraRes_KnownNonNeighborsButWithinDistance,
                    nonBondedInterRes_WithinDistance;
      
  void issue_warnings();
  string stringid(SiteID atom) const;
  int user_interaction(SiteID atom1, SiteID atom2) const;
};

#endif /*COVALENTBONDER_H_*/
