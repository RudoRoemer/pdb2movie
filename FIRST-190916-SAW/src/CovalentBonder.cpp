#include "CovalentBonder.h"
#include "global_defs.h"
#include "generalUtils.h"

#define HARD_DISTANCE_CUTOFF 6.00

extern const Parameters parameters;

size_t BondDictionary::Atom::numNeighbors() const {
  return neighborset.size();
}

bool BondDictionary::Atom::isNeighbor(atomname_t atomname) const {
  return neighborset.find(atomname) != neighborset.end();
}

BondDictionary::Atom *BondDictionary::Res::getAtom(atomname_t atomname)  {
  atommap_t::iterator iatom = atommap.find(atomname);
  return (iatom != atommap.end()) ? &iatom->second : NULL;
}

BondDictionary::Res *BondDictionary::getRes(resname_t resname)  {
  resmap_t::iterator ires = resmap.find(resname);
  return (ires != resmap.end()) ? &ires->second : NULL;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   For each unique residue name in molFramework, check to see if 
//   it is already in this BondDictionary object.
//   If not, search the het group dictionary file, which will add the residue's 
//   template to this BondDictionary object.
// Parameters:
//   filename - the het dictionary filename
//   molFramework - the molFramework object whose het residues are to be added
//     to this BondDictionary object.
////////////////////////////////////////////////////////////////////////////////
void BondDictionary::identify_het_residues(string filename, MolFramework* molFramework){
  
  resname_t* this_resname;
  set<resname_t> residue_set;
  
  // First, create a set of the unique residue names.
  // A set does not hold duplicate entries.
  for( SiteID current_atom = 1; current_atom <= molFramework->total_sites; current_atom++ ){
    this_resname = &(molFramework->site_info[current_atom].residue_name);
    residue_set.insert( *this_resname );
  }
  
  bool found;
  set<resname_t>::iterator ires;
  // For each unique residue name, check to see if it is already in this dictionary
  // object.  If not, search the het dictionary file.
  for( ires = residue_set.begin(); ires != residue_set.end(); ires++) {
    found = resmap.find(*ires) != resmap.end();
    if(!found) 
      search_het_group_dictionary( filename, *ires );
  }
      
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The "standard" templates refer to the 20 standard amino acids, and a 
//   few of the most common HETATM residues, such as water. These common 
//   residues are read and stored in this BondDictionary object. 
// Parameters:
//   filename - the standard template file              
////////////////////////////////////////////////////////////////////////////
void BondDictionary::read_standard_residue_library(string filename){

  string linebuf;

  //    Open the file
  fstream res_library( filename.c_str(), ios::in );
  if( !res_library ){
    cout << " ERROR: Residue Library " << endl << "  " << filename << endl << "was not found." << endl;
    exit(1);
  }
  
  //    Store each residue template  
  while( !res_library.eof() ){
    getline( res_library, linebuf );

    if( !linebuf.find("RESIDUE") )
      store_residue_template_data( &res_library, linebuf.substr(10,3) );
  }

  res_library.close();
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Search the het dictionary file 
//   for a RESIDUE record corresponding to the residue_name passed 
//   to this function as a parameter. If the residue is found, store the 
//   data in this BondDictionary object.
// Parameters:
//   filename - the het dictionary file
//   residue_name - Three letter name of a residue 
// Return Value List:
//   true - If the residue is found in the file.
//   false - if the residue is not found in the file.
////////////////////////////////////////////////////////////////////////////
bool BondDictionary::search_het_group_dictionary( string filename, resname_t residue_name ){

  // 1. Open the HET group dictionary.
  //////////////////////////////////////////////////////////////////////

  fstream het_dictionary( filename.c_str(), ios::in );
  if( !het_dictionary ){
    cout << " ERROR: Residue Library " << endl << "  " << filename << endl << "was not found." << endl;
    exit(1);
  }

  // 2. Search the dictionary file for the residue named "residue_name".
  //    Note: Previous logic would search for residue_name in the
  //    RESIDUE line.  But this would lead to a false match in certain
  //    cases.  For example, when searching for residue ZN, the template
  //    with residue named AZN would be mistakenly matched.
  //    So, this modified logic forces an exact match, ignoring spaces.
  //////////////////////////////////////////////////////////////////////
  string linebuf;
  string thisline_resname;
  string target_resname = residue_name;
  removeWhiteSpacePadding(target_resname);
  while( !het_dictionary.eof() ){
    getline( het_dictionary, linebuf );
    if( linebuf.find("RESIDUE") != 0 ) continue;
    thisline_resname = linebuf.substr(10,3);
    removeWhiteSpacePadding(thisline_resname);
    
    if( thisline_resname == target_resname){ 
      //cout << "Found HET residue |" << residue_name << "|" << endl;
      store_residue_template_data( &het_dictionary, residue_name );
      het_dictionary.close();
      return(true);
    }
  }
  het_dictionary.close();

  // 3. If the residue name is no found in the dictionary, return 0.
  //////////////////////////////////////////////////////////////////////
  return(false);
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Read and store the data of a RESIDUE record in the file passed to this function.
//   The records consist of CONECT data that describe the exact topology of a residue
//   using standard atom names. The CONECT records are stored in a data structure of
//   of type "template_data". When all the CONECT records have been read in, the 
//   data is pushed onto to a vector. 
// Parameters:
//   file_ptr - Pointer to an open input file stream.
//   residue_name - name of the RESIDUE record for which data is being stored.
////////////////////////////////////////////////////////////////////////////
void BondDictionary::store_residue_template_data( fstream *file_ptr, resname_t residue_name ){

  unsigned int neighbor_start = 0;
  size_t end_of_line = 0; // unsigned int here breaks on flexweb (gcc 3.4.3 20050227 (Red Hat 3.4.3-22.1))

  string linebuf;
  atomname_t this_atom_name;
  atomname_t neighbor_name;
  size_t found_str ;
  int inum =0;


  while( linebuf.find("RESIDUE") == string::npos ){
    getline( *file_ptr, linebuf );
    found_str = linebuf.find("RESIDUE");
    if( found_str !=  string::npos ){
       string numb = linebuf.substr(16,4);
       inum = atoi(numb.c_str() );
    }
  
    if( !linebuf.find("CONECT") || inum == 1 ){ // We only need the CONECT data. In case of no CONECT records residue_name should be  this_atom_name.
      
      if (inum != 1){
	this_atom_name = linebuf.substr(11,4);
	removeWhiteSpacePadding( this_atom_name );
      } 
      else {
	this_atom_name = residue_name;
	removeWhiteSpacePadding( this_atom_name );
      }
     
      //Insert this_atom key into dictionary.
      //Even if this_atom has zero neighbors listed, we still want
      //it to be found in the dictionary.
      
      BondDictionary::Atom* dict_atom;
      dict_atom = &resmap[residue_name].atommap[this_atom_name]; // the [] notation adds
      // residue_name and this_atom_name to the corresponding map variables.  By adding these
      // to the maps, default Res and Atom objects are instantiated.
      // This line also assigns to "neighbors" a pointer to the corresponding 
      // Atom object.
      
      neighbor_start = 20; // The atom neighbors start in column 20 of the template file. 
      end_of_line = linebuf.find_first_not_of(" \n", neighbor_start ); // don't store a newline character.
      
      // Read and store the neighbors of "this_atom" until we reach a newline. 
      //////////////////////////////////////////////////////////////////////
      while( end_of_line != string::npos ){ 
        neighbor_name = linebuf.substr( neighbor_start, 4 );
        neighbor_start += 5;
        end_of_line = linebuf.find_first_not_of(" \n", neighbor_start-1 );
        removeWhiteSpacePadding( neighbor_name );

        //  add the neighbor to the dictionary 
        //////////////////////////////////////////////////////////////////////
        dict_atom->neighborset.insert( neighbor_name );
      }      
    }
  }
}
////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////
// Constructor:
//   Initialize the BondDistanceTable object with the bond distance data
//   found in the file passed as an argument.
// Parameters:
//   filename - the bond distances file
////////////////////////////////////////////////////////////////////////////
BondDistanceTable::BondDistanceTable(string filename) {
  read_file( filename );
}

BondDistanceTable::const_iterator BondDistanceTable::find( elemname_t elem1, elemname_t elem2 ) const {
  string this_bond;
  // 1. Create a string listing the two elements in the bond. 
  ////////////////////////////////////////////////////////////////////////////////
  
  if( elem1 > elem2 )
    this_bond = elem2 + elem1;
  else
    this_bond = elem1 + elem2;
  
  // 2. Get the stored distance cutoff for this bond type. 
  ////////////////////////////////////////////////////////////////////////////////
  return distmap.find( this_bond );

}

BondDistanceTable::const_iterator BondDistanceTable::begin() const {
  return distmap.begin();
}

BondDistanceTable::const_iterator BondDistanceTable::end() const {
  return distmap.end();
}

////////////////////////////////////////////////////////////////////////////
// Description:
//   Read and store a table atom-atom distance cutoffs. These cutoffs are an
//   upper bound on the distance when identifying covalent bonds.
////////////////////////////////////////////////////////////////////////////
void BondDistanceTable::read_file(string filename) {

  string linebuf;
  string bond;
  string temp;
  float cutoff = 0;

  // 1. Open the file with the distance cutoff's listed.
  //////////////////////////////////////////////////////////////////////
  
  fstream cutoff_table( filename.c_str(), ios::in );
  if( !cutoff_table ){
    cout << "ERROR: Could not open file distance_cutoffs.txt. " << endl;
    exit(1);
  }
  
  // 2. Read each non-blank line of the distance cutoff data file, and 
  //    store the data in a map that returns a float (the cutoff) when 
  //    indexed by a string (the elements in the bond). 
  //////////////////////////////////////////////////////////////////////
  while( !cutoff_table.eof() ){
    getline( cutoff_table, linebuf );
    if( linebuf.length() ){
      bond   = linebuf.substr( 0,4 );
      temp   = linebuf.substr( 6,5 );
      cutoff = atof( temp.c_str() );

      distmap[bond] = cutoff;    
    }
  }   

  cutoff_table.close();

}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Constructor
//   Fill the BondDictionary objectwith all required standard and het templates.
//   Read data into the BondDistanceTable object.
//   Read data into the MetalList object.
//   Make one pass over the atoms, recording information about each atom,
//   such as which atoms are recognized and which are not, and storing 
//   a pointer to each atom's neighbor data in the templates.
////////////////////////////////////////////////////////////////////////////////
CovalentBonder::CovalentBonder( MolFramework *in_molFramework ) {

  molFramework = in_molFramework;
  site_info = molFramework->site_info;
  
  // initialize the bonddict, bondDistanceTable
  //////////////////////////////////////////////////////////////////////
  string std_res_filename = parameters.path + "/lib/std_dictionary.txt";
  string het_filename= parameters.path + "/lib/het_dictionary.txt";
  string distance_cutoff_filename = parameters.path + "/lib/distance_cutoffs.txt";

  bonddict.read_standard_residue_library(std_res_filename);
  bonddict.identify_het_residues( het_filename, molFramework);
  //bonddict.display(); //For debugging, to view the dictionary
  bondDistanceTable.read_file( distance_cutoff_filename );
  
  // allocate space
  dict_atom_lookup.resize(molFramework->total_sites +1, NULL);
  
  // for each atom, record whether its resname was found,
  // whether its atomname was found, and a pointer to its neighbor template
  // data so that we don't have to keep looking the atom up over and over.
  //////////////////////////////////////////////////////////////////////
  BondDictionary::Res *dict_res;
  BondDictionary::Atom *dict_atom;

  for (SiteID atom = 1; atom <= molFramework->total_sites; atom++) {

    //look up residue in Bond Dictionary
    //////////////////////////////////////////////////////////////////////
    dict_res = bonddict.getRes( site_info[atom].residue_name );
    
    //if residue was found, look up the atom
    //////////////////////////////////////////////////////////////////////
    //? if( site_info[atom].element_name != "H " && 
    //?	dict_res != NULL )
    //KS. 03/18/10 If we found hydrogen in std/het dictionary we know that this hydrogen
    //is bonded 

    if( dict_res != NULL )
      dict_atom = dict_res->getAtom( site_info[atom].atom_name );
    else
      dict_atom = NULL;
      
    // whether the atom was found or not, store the neighbors pointer.  
    // It will be null if the atom was not found.
    //////////////////////////////////////////////////////////////////////
    dict_atom_lookup[atom] = dict_atom;
    
    //record atom names that are not found in the templates.
    //////////////////////////////////////////////////////////////////////
    if( dict_atom == NULL &&
	site_info[atom].element_name != "H " )
      missingAtomNames.insert( pair<resname_t, atomname_t> (site_info[atom].residue_name, site_info[atom].atom_name) ); 
  }
  
  //conect_info.get_from_molFramework( molFramework );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void CONECT_info::get_from_molFramework( MolFramework *molFramework ) {
  
  /*
  conect_record_neighborlist.resize( molFramework->total_sites + 1 );
  
  typedef vector<SiteID> conectlist_t;
  typedef map< SiteID, conectlist_t> conectmap_t;
  conectmap_t *conect_records = &molFramework->conect_records;
  conectmap_t::iterator iconectrecord;
  conectlist_t* list_of_conect_neighbors;
  conectlist_t::iterator ineighbor;
  SiteID atom1 = 0;
  SiteID atom2 = 0;
  
  for ( iconectrecord = conect_records->begin(); iconectrecord != conect_records->end(); iconectrecord++ ){

    atom1 = molFramework->getFIRSTNumber( iconectrecord->first );
    if( atom1 != 0 ){

      list_of_conect_neighbors = &iconectrecord->second;
      for ( ineighbor = list_of_conect_neighbors->begin(); ineighbor != list_of_conect_neighbors->end(); ineighbor++ ) {
	
	atom2 = molFramework->getFIRSTNumber( *ineighbor );
		
	if( atom2 != 0 ){
	  cout << "insert " << atom1 << "-" << atom2 << endl;
	  conect_record_neighborlist[atom1].insert(atom2);
	  conect_record_neighborlist[atom2].insert(atom1);
	}
      }
      cout << endl;
    }
  }
  */
}

////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////
CovalentBonder::~CovalentBonder()
{
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
//   This function iterates over atom pairs and bonds the pairs that should
//   be bonded.  The bonding logic is summarized as follows:
//
//   For each pair of atoms, determine if the pair is INTRA-residue or
//   INTER-residue.
//
//   INTRA-residue bonding is described in the templates.  So, look up 
//   the pair in the templates and bond the atoms if they are listed as
//   neighbors.  If the atoms are found in the templates but are not listed
//   as neighbors, do not bond them.  If at least one of the atoms is not found 
//   in the templates, then we have a problem with atom names, and so we cannot
//   tell from the templates whether the atoms should be bonded.  So, for this
//   case, bond the atoms based on a distance check.  Also, if a CONECT record
//   is found for this intra-residue pair, and if the distance cutoff is met,
//   bond the atoms.
//
//   INTER-residue bonding is not described in the templates, so we cannot simply
//   look up the answer.  But, we also cannot simply bond atoms based on a 
//   distance check, because very often PDB files have atoms placed within
//   bonding distance that really are not supposed to be bonded.  This happens
//   most often with added hydrogens.
//
//   So, to determine inter-residue bonding,
//   we see if the atoms pass a test:  Do the atoms form a peptide bond?
//   Do the atoms form an inter-nucleic acid bond?  Do the atoms form a 
//   disulfide bond?  If so, bond the atoms.  
// 
//   All bonds listed in CONECT record fields before column 32 (this excludes
//   hydrogen bond records) will be added to the network be default, unless
//   the bond has already been identified. No distance check is performed
//   when adding CONECT record bonds. 
////////////////////////////////////////////////////////////////////////////
void CovalentBonder::autobond(){

  vector<SiteID>::iterator iatom2;
  vector<int>::iterator current_grid;
  vector<int> grids_to_check;
  vector<SiteID> *atomsInGrid;
  
  SiteID atom1 = 0;
  SiteID atom2 = 0;
  bool templates_confirm_neighbors;
  bool templates_confirm_nonneighbors;
  bool cannot_discern_bond_from_templates;
  bool passes_distance_check;
  int bars;
     
  for( atom1 = 1; atom1 <= molFramework->total_sites; atom1++ ) {

    if( parameters.verbose >= 2 && 
	!(atom1 % 1000) )
      cout << "  Checking atom " << setw(8) << atom1 << "/" << setw(8) << molFramework->total_sites << endl;
    
    grids_to_check = molFramework->getGrids( atom1 );    
    for ( current_grid = grids_to_check.begin(); current_grid != grids_to_check.end(); current_grid++ ){

      atomsInGrid = &molFramework->coordinate_grid[*current_grid];
      for( iatom2 = atomsInGrid->begin(); iatom2 != atomsInGrid->end(); iatom2++){

        if( *iatom2 <= atom1 ) 
	  continue;

	atom2 = *iatom2;
        
        bars = 0;
        passes_distance_check = meets_atom_atom_distance_cutoff( atom1, atom2 );
	
        //check for intra-residue bonds
	//////////////////////////////////////////////////////////////////////
        if ( !molFramework->isDifferentResidue(atom1, atom2) ) {

	  // HARD_DISTANCE_CUTOFF is defined at the top of this file
	  //////////////////////////////////////////////////////////////////////
          templates_confirm_neighbors = ( dict_atom_lookup[atom1] != NULL && 
					  dict_atom_lookup[atom1]->isNeighbor(site_info[atom2].atom_name) && 
					  dict_atom_lookup[atom2] != NULL );
          templates_confirm_nonneighbors = !templates_confirm_neighbors && dict_atom_lookup[atom1]!=NULL && dict_atom_lookup[atom2]!=NULL;
          cannot_discern_bond_from_templates = !templates_confirm_nonneighbors && !templates_confirm_neighbors;
                    
          if( templates_confirm_neighbors ){
	    if( molFramework->getDistance(atom1, atom2) > HARD_DISTANCE_CUTOFF ){
	      cout << "ERROR: The bond distance between intraresidue atoms " << site_info[atom1].orig_atom_number 
		   << " and " << site_info[atom2].orig_atom_number << " exceeds " << HARD_DISTANCE_CUTOFF << " Angstroms" << endl;
	      cout << "       This bond is defined in the residue template file. Please verify " << endl;
	      exit(1);
	    }
	    if ( is_locked(atom1, atom2) ) 
	      bars = 6;
	    else 
	      bars = 5;
	  }
          else if ( cannot_discern_bond_from_templates && 
		    passes_distance_check ) {
            bars = 5;
          }
          
          // Look for potential problems 
	  //////////////////////////////////////////////////////////////////////
          if ( templates_confirm_neighbors && !passes_distance_check ) {
            //Here, the atoms are known to be neighbors, as confirmed in the templates.
            //However, the bonding distance cutoff is not met.  We have already
            //bonded the atoms, but we should record this irregularity and notify
            //the user.
            bondedIntraRes_KnownNeighborsButNotWithinDistance.push_back( pair<SiteID, SiteID> (atom1, atom2) );
          }
          if ( templates_confirm_nonneighbors
                && passes_distance_check){
            //Here, we know that the atoms are not neighbors
            //So the atoms have not been bonded,
            //yet they are close enough to be bonded.
            //We should record this irregularity and notify the user.
            nonbondedIntraRes_KnownNonNeighborsButWithinDistance.push_back( pair<SiteID, SiteID> (atom1, atom2) );
          }
          
        } //end intra-residue check


        //check for inter-residue bonds
	//////////////////////////////////////////////////////////////////////
        else { 

	  //if( conect_info.is_listed_in_conect_records(atom1, atom2) ){
	  //bars = 5;
	  //}
	  
          if ( passes_distance_check ) {

            if ( labelled_as_peptide_bond(atom1, atom2) ){ 
	      bars = 6; 
	    }
            else if( labelled_as_disulfide_bond(atom1, atom2) ) {
              if ( parameters.lock_disulfide ) 
		bars = 6;
              else 
		bars = 5;
            }
            else if( labelled_as_internucleic_bond(atom1, atom2) ){
	      bars = 5;
	    }
            else {
              //Here, the atoms are in different residues.  We have checked
              //the conditions for inter-residue bonding, but have not found a reason
              //to justify the bond.  However, the atoms are within bonding distance.
              //This typically indicates a severe steric clash, so it is correct
              //not to bond the atoms.  But, we should record this irregularity and
              //notify the user.
              nonBondedInterRes_WithinDistance.push_back( pair<SiteID, SiteID> (atom1, atom2) );
            }
          }

        } //end inter-residue check
          
        
        //if a bond between the pair was detected, write the bond to the structure  
        if( bars ){
	  //cout << "added bond" << endl;
          molFramework->add_to_site_info_array(atom1, atom2, bars);
        }
        
      } //end loop over atoms in this grid
    } // end loop over grids for this atom
  } //end loop over atoms

  // Add all bonds listed in the CONECT records by default, unless the bond
  // already exists in the network, then do nothing. This is to prevent
  // overwriting bonds that may have been added with 6 bars, as the CONECT
  // records provide no information on how many bars to model a bond with.
  //////////////////////////////////////////////////////////////////////
  map< SiteID, vector<SiteID> >::iterator conectRecordsIter = molFramework->conect_records.begin();

  while( conectRecordsIter != molFramework->conect_records.end() ){

    SiteID baseAtom = molFramework->orig_2_FIRST_atom_number[conectRecordsIter->first];

    if( baseAtom != 0 ){
      //cout << "ERROR: Base atom in CONECT record not found." << endl;
      //cout << "       Atom [" << conectRecordsIter->first << "] was listed in a CONECT record" << endl;
      //cout << "       in the input file, but no corresponding atom was found in the file." << endl;
      //exit(1);
      
      vector<SiteID> neighbors = conectRecordsIter->second;

      for( SiteID currentNeighbor = 0; currentNeighbor < neighbors.size(); currentNeighbor++ ){
	
	SiteID neighborAtom = molFramework->orig_2_FIRST_atom_number[ neighbors[currentNeighbor] ];     
	
	if( neighborAtom != 0 ){

	  //cout << "ERROR: Neighbor atom in CONECT record not found." << endl;
	  //cout << "       Atom [" << currentNeighbor << "] was listed in a CONECT record" << endl;
	  //cout << "       in the input file, but no corresponding atom was found in the file." << endl;
	  //exit(1);
	  
	  if( parameters.verbose >= 3 ){
	    cout << "    Adding CONECT between " << baseAtom << " " << neighborAtom << endl;
	  }
	  
	  if( !molFramework->isBonded( baseAtom, neighborAtom ) ){
	    molFramework->add_to_site_info_array( baseAtom, neighborAtom, 5 );
	    //KS. 02/12/2010. Check if we can find the pair of atoms in nonBondedInterRes_WithinDistance.
	    // If this pair is in the list we have to remove it, so program would not give the warning to the users
	    ListOfSiteIDPairs::iterator ipair;
	    ipair = find(nonBondedInterRes_WithinDistance.begin(), nonBondedInterRes_WithinDistance.end(),  pair<SiteID, SiteID> (baseAtom, neighborAtom) );
	    if (ipair != nonBondedInterRes_WithinDistance.end() ){
	       nonBondedInterRes_WithinDistance.erase(ipair);    
	       //   cout << " nonBondedInterRes_WithinDistance erased pair " << baseAtom << " " << neighborAtom <<endl;
	    }

	  }
	}
      }
    }

    conectRecordsIter++;
  }

  //////////////////////////////////////////////////////////////////////  
  issue_warnings();

  if ( parameters.interactive ) {
    ListOfSiteIDPairs::const_iterator ipair = nonBondedInterRes_WithinDistance.begin();
    while( ipair != nonBondedInterRes_WithinDistance.end() ) {
      if ( user_interaction(ipair->first, ipair->second) ) {
        bars=5;
        molFramework->add_to_site_info_array(ipair->first, ipair->second, bars);
      }
      ipair++;
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void CovalentBonder::issue_warnings() {
  cout << endl;
  
  bool anywarnings = false;
  ostringstream o;
  
  if( missingAtomNames.size() != 0) {

    anywarnings = true;

    cout << "WARNING: Some atom names are not found in the PDB standard/het dictionary." << endl;
    
    o << "Some atom names are not found in the PDB standard/het dictionary.  The" << endl;
    o << "dicitonaries contain bonding connectivity for all residues and het groups" << endl;
    o << "in the PDB.  Possible reasons why your atom name is not found are:" << endl; 
    o << "  i) You used what-if to add hydrogens (it changes atom names) " << endl;
    o << "  ii) You have renamed your atoms in some other program" << endl;
    o << "  iii) You have a hand-made het group or residue that is not from the PDB" << endl;
    o << "Be aware that INTRA-residue bonding involving these atoms are distance-based," << endl;
    o << "and bonds are assumed to be rotatable." << endl;
    o << "Also, when FIRST looks for INTER-residue bonds, it looks for particular " << endl;
    o << "atom names.  If these atoms are supposed to form INTER-residue, FIRST may" << endl;
    o << "not be able to auto-detect them." << endl << endl;
    set<pair< resname_t, atomname_t> >::const_iterator ipair;
    for (ipair = missingAtomNames.begin(); ipair != missingAtomNames.end(); ipair++) {
      o << "|" << ipair->first << "|" << ipair->second << "|" << endl;
    }
    o << endl;
    molFramework->add_warning(2, o.str());
    o.str("");
  }
  
  if( missingBondDistances.size() != 0 ) {
    anywarnings = true;
    cout << "WARNING: Some element pairs were not found in the Bond Distance Table." << endl;
    
    o << "Some element pairs were not found in the Bond Distance Table. " << endl;
    o << "Atom pairs involving these elements were left unbonded, unless" << endl;
    o << "they were known to be neighbors from a particular residue template." << endl;
    set< pair<elemname_t, elemname_t> >::const_iterator ipair;
    for (ipair = missingBondDistances.begin(); ipair != missingBondDistances.end(); ipair++) {
      o << "|" << ipair->first << "|" << ipair->second << "|" << endl;
    }
    o << endl;
    molFramework->add_warning(2, o.str());
    o.str("");
  }
  
  if (bondedIntraRes_KnownNeighborsButNotWithinDistance.size() != 0 ) {
    anywarnings = true;
    cout << "WARNING: Some bonded intra-residue pairs do not meet distance cutoff." << endl;

    o << "Some bonded intra-residue pairs do not meet distance cutoff." << endl;
    o << "They are bonded because the standard/het templates list them as" << endl;
    o << "neighbors." << endl;
    ListOfSiteIDPairs::const_iterator ipair = bondedIntraRes_KnownNeighborsButNotWithinDistance.begin();
    for (; ipair != bondedIntraRes_KnownNeighborsButNotWithinDistance.end(); ipair++) {
      o << stringid(ipair->first) << "  to  " << stringid(ipair->second) << endl;
    }
    o << endl;
    molFramework->add_warning(2, o.str());
    o.str("");
  }

  if (nonbondedIntraRes_KnownNonNeighborsButWithinDistance.size() != 0 ) {
    anywarnings = true;
    cout << "WARNING: Some pairs of atoms are within bonding distance but are" << endl;
    cout << "         left unbonded because they are known not to be neighbors." << endl;

    o << "Some pairs of atoms are within bonding distance but are" << endl;
    o << "left unbonded because they are known not to be neighbors." << endl;
    ListOfSiteIDPairs::const_iterator ipair = nonbondedIntraRes_KnownNonNeighborsButWithinDistance.begin();
    for (; ipair != nonbondedIntraRes_KnownNonNeighborsButWithinDistance.end(); ipair++) {
      o << stringid(ipair->first) << "  to  " << stringid(ipair->second) << endl;
    }
    o << endl;
    molFramework->add_warning(2, o.str());
    o.str("");
  }

  if( nonBondedInterRes_WithinDistance.size() != 0 ) {
    anywarnings = true;
    cout << "WARNING: Some pairs of atoms from different residues are within " << endl;
    cout << "         bonding distance, but no justification can be found for " << endl;
    cout << "         bonding them." << endl;

    o << "Some pairs of atoms from different residues are within bonding distance," << endl;
    o << "but no justification can be found for bonding them.  So, FIRST assumes " << endl;
    o << "that these atoms are in a severe steric clash, and that they should" << endl;
    o << "remain unbonded.  The atoms did not pass FIRST's check for" << endl;
    o << "peptide bond, internucleic bond, disulfide bond, or CONECT record." << endl;
    o << "Some of these checks can fail if the atom names are not found in" << endl;
    o << "the standard and or het-group templates." << endl;
    o << "In interactive mode, FIRST prompts you whether to bond these atoms." << endl;
    ListOfSiteIDPairs::const_iterator ipair = nonBondedInterRes_WithinDistance.begin();
    for (; ipair != nonBondedInterRes_WithinDistance.end(); ipair++) {
      o << stringid(ipair->first) << "  to  " << stringid(ipair->second) << endl;
    }
    o << endl;
    molFramework->add_warning(2, o.str());
    o.str("");
  }
  
  if (anywarnings) {
    cout << endl << "See the *_results.txt file for lists of atoms that cause the warnings," << endl;
    cout << "and more information regarding the warnings. " << endl << endl;
    if (parameters.interactive) {
      string input;
      cout << "Press Enter to continue running FIRST..." << endl;
      getline(cin, input);
    }
  }
  
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Check whether atoms are labelled as a disulfide bridge between two cystine 
//   residues. As a validation, the sulfur should be in the thiol form or 
//   disulfide form, but not both.
// Parameters:
//   atom_1 - FIRST_number of an atom. 
//   atom_2 - FIRST_number of an atom.
// Return Value:
//   true - Atoms are both named SG from residue CYS
//   false - Failed criteria.
////////////////////////////////////////////////////////////////////////////////
bool CovalentBonder::labelled_as_disulfide_bond( SiteID atom_1, SiteID atom_2 ) const {

  return ( site_info[atom_1].atom_name == "SG" &&
	   ( site_info[atom_1].residue_name == "CYS" || site_info[atom_1].residue_name == "CYX" ) &&
	   site_info[atom_2].atom_name == "SG" &&
	   ( site_info[atom_2].residue_name == "CYS" || site_info[atom_2].residue_name == "CYX" ) );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The internucleic acid bond connection (O3' to P) is not listed in the standard 
//   residue library templates. Here, we'll check to see if both atoms 
//   have the right atom names.
////////////////////////////////////////////////////////////////////////////////
bool CovalentBonder::labelled_as_internucleic_bond( SiteID atom_1, SiteID atom_2 ) const {

  return(( site_info[atom_1].atom_name == "O3'" && site_info[atom_2].atom_name == "P" ) ||
	 ( site_info[atom_1].atom_name == "P" && site_info[atom_2].atom_name == "O3'" )); 
}


// Description: This relies on conect_record_neighborlist being symmetric,
//   or in other words, a pair can be found by looking it up with atom1 or atom2.
//   Since ordering doesn't matter, it is sufficient to check if atom2 is listed
//   as a neighbor of atom1.
////////////////////////////////////////////////////////////////////////////////
bool CONECT_info::is_listed_in_conect_records(SiteID atom1, SiteID atom2) const {

  if( conect_record_neighborlist[atom1].find(atom2) != conect_record_neighborlist[atom1].end() )
    cout << "     A: found " << atom2 << " in " << atom1 << endl;
  if( conect_record_neighborlist[atom2].find(atom1) != conect_record_neighborlist[atom2].end() )
    cout << "     B: found " << atom1 << " in " << atom2 << endl;

  return conect_record_neighborlist[atom1].find(atom2) != conect_record_neighborlist[atom1].end();
}

//////////////////////////////////////////////////////////////////////////////// 
// Description:
////////////////////////////////////////////////////////////////////////////////
bool CovalentBonder::is_locked( SiteID atom1, SiteID atom2 ) const {

  int mult_atom1 = 0;
  int mult_atom2 = 0;
  
  if( site_info[atom1].element_name > site_info[atom2].element_name ){
    swap( atom1, atom2 );
  }
  
  if (dict_atom_lookup[atom1] == NULL || dict_atom_lookup[atom2] == NULL) return false;
  mult_atom1 = dict_atom_lookup[atom1]->numNeighbors();
  mult_atom2 = dict_atom_lookup[atom2]->numNeighbors();
  
  if( site_info[atom1].element_name == "C " &&
      mult_atom1 <= 3 &&
      site_info[atom2].element_name == "C " &&
      mult_atom2 <= 3 ){
    return( true );
  }
  else if( site_info[atom1].element_name == "C " && mult_atom1 <= 3 &&
     site_info[atom2].element_name == "N " && mult_atom2 <= 3 ){
    //cout << "locking " << site_1 << " " << site_2 << endl;
    return( true );
  }
  
  return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The peptide bond connection is not listed in the standard residue
//   library templates. Here, we'll check to see if both atoms are the 
//   proper distance apart, and have the right atom names.
////////////////////////////////////////////////////////////////////////////////
bool CovalentBonder::labelled_as_peptide_bond( SiteID atom_1, SiteID atom_2 ) const {

  return (( site_info[atom_1].atom_name == "C" && site_info[atom_2].atom_name == "N" ) ||
	  ( site_info[atom_1].atom_name == "N" && site_info[atom_2].atom_name == "C" ));
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool CovalentBonder::meets_atom_atom_distance_cutoff( SiteID atom1, SiteID atom2 ) {
  
  elemname_t elem1 = site_info[atom1].element_name;
  elemname_t elem2 = site_info[atom2].element_name;
  BondDistanceTable::const_iterator iTableEntry = bondDistanceTable.find( elem1, elem2 );  

  if (iTableEntry != bondDistanceTable.end()) {
    return ( molFramework->getDistance(atom1, atom2) <= iTableEntry->second );
  }
  else{
    if (elem1 > elem2) swap(elem1, elem2);
    missingBondDistances.insert( pair<elemname_t, elemname_t> (elem1, elem2) );
    return( false );
  }
}

////////////////////////////////////////////////////////////////////////////////
string CovalentBonder::stringid(SiteID atom) const {
  ostringstream o;
  o << setw(7) << site_info[atom].orig_atom_number
    << setw(5) << molFramework->atomNamePDBFormat(site_info[atom].atom_name, site_info[atom].element_name)
    << setw(4) << site_info[atom].residue_name
    << setw(2) << char(site_info[atom].chain_ID)
    << setw(4) << site_info[atom].seq_number;
 
  return o.str();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
int CovalentBonder::user_interaction( SiteID atom_1, SiteID atom_2 ) const {
  string input;
  int choice = 1;
  clear_screen;
  
  cout << " -- Option - Bond Between Residues -----" << endl << endl;
  cout << "The following atoms from two different residues were found to be within" << endl;
  cout << "covalent bonding distance.  The atoms did not pass FIRST's check for" << endl;
  cout << "peptide bond, internucleic bond, or disulfide bond.  Checks for these" << endl;
  cout << "kinds of bonds can fail if the atom names do not match the standard and" << endl;
  cout << "het templates.  Any other kind of inter-residue bond can only be auto-" << endl;
  cout << "detected if a CONECT record is found in the PDB file.  " << endl << endl;
  cout << "It may be that these atoms are in a severe steric clash, and that they should" << endl;
  cout << "remain unbonded.  " << endl;
  cout.setf( ios::left );

  cout << spacing
       << setw(10) << "atom"
       << setw(10) << "atom"
       << setw(10) << "residue"
       << setw(10) << "sequence"
       << setw(8)  << "chain" << endl;

  cout << spacing
       << setw(10) << "number"
       << setw(10) << "name"
       << setw(10) << "name"
       << setw(10) << "number"
       << setw(8)  << "ID" << endl << endl;

  cout << spacing 
       << setw(10) << site_info[atom_1].orig_atom_number
       << setw(10) << site_info[atom_1].atom_name 
       << setw(10) << site_info[atom_1].residue_name
       << setw(10) << site_info[atom_1].seq_number
       << setw(8)  << char(site_info[atom_1].chain_ID) << endl

       << spacing << "           |" << endl

       << spacing 
       << setw(10) << site_info[atom_2].orig_atom_number
       << setw(10) << site_info[atom_2].atom_name 
       << setw(10) << site_info[atom_2].residue_name
       << setw(10) << site_info[atom_2].seq_number
       << setw(8)  << char(site_info[atom_2].chain_ID)
       << endl << endl

       << spacing << "Distance = " << setprecision(4) << molFramework->getDistance(atom_1, atom_2) << " Angstroms" << endl << endl;
  
  do{
    cout << " 1. CONNECT the atoms (Press Enter to connect)." << endl;
    cout << " 2. DO NOT CONNECT the atoms." << endl << endl;
    cout << "    Choice = ";
    getline(cin, input);
    choice = (int) strtol(input.c_str(), 0, 10);

    if( choice == 1 ||
         did_press_enter(input) ){
      return(1);
    }
    else if( choice == 2 )
      return(0);
    else
      cout << " Please choose 1 or 2." << endl;
    

    cout << "choice = " << choice << endl;
  } while( choice <= 0 || choice > 2 );

  return(1);
}
