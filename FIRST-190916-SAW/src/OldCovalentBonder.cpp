#include "../include/OldCovalentBonder.h"
#include "../include/global_defs.h"
#include "generalUtils.h"

extern const Parameters parameters;

OldCovalentBonder::OldCovalentBonder(MolFramework * in_molFramework)
{
molFramework = in_molFramework;
site_info = molFramework->site_info;
}

OldCovalentBonder::~OldCovalentBonder()
{
}

void OldCovalentBonder::autobond(){
identifyCovalentBonds();
};

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Three checks are completed to identify covalent bonds. The first 
//   check tries to find the current residue in the standard template library. 
//   If the residue is found, the atom is sent to the function find_bonds(), 
//   that will, well, find covalent bonds to this atom. If the residue is
//   not in the standard template library, the HET group dictionary is 
//   opened and searched for the residue name. If the residue is found in 
//   the HET dictionary, the residue template data is read into the 
//   het_group_template, and both the template and the atom are passed to 
//   find_bonds(). 
////////////////////////////////////////////////////////////////////////////////
void OldCovalentBonder::identifyCovalentBonds(){


  bool found = false;
  bool in_standard_template = false;
  bool found_template = false;

  // 1. Read the residues listed in the standard group dictionary into an 
  //    array of structures. ** process_file.cpp **
  ////////////////////////////////////////////////////////////////////////////////
  read_standard_residue_library();
  read_distance_cutoff_table();

  // 2. Connect each atom based on a standard group template, a HET group 
  //    template or a simple distance-based check.
  //////////////////////////////////////////////////////////////////////
  for( unsigned int current_atom = 1; current_atom <= molFramework->total_sites; current_atom++ ){

    if( parameters.verbose >= 2 && 
        !(current_atom%1000) )
      cout << "   Checking atom " << setw(8) << current_atom << "/" << setw(8) << molFramework->total_sites << endl;
    
    // Check to see if we have already failed to find this residue name. Don't 
    // want to keep checking the HET dictionary, it's big.    
    ////////////////////////////////////////////////////////////////////////////////
    if( find(unknown_residue.begin(), unknown_residue.end(), site_info[current_atom].residue_name) == unknown_residue.end() ){
      
      // a1. Check the standard residue dictionary. 
      //////////////////////////////////////////////////////////////////////
      for( unsigned int siteNumber = 0; siteNumber < standard_template.size(); siteNumber++ ){
  
        if( site_info[current_atom].residue_name == standard_template[siteNumber].std_group_name ){
        //cout << "[" << site_info[current_atom].residue_name << "] [" << standard_template[siteNumber].std_group_name << "]" << endl;
          in_standard_template = true;
          found = find_bonds( siteNumber, current_atom );
          break;
        }
      }
      
      // Check the heterogroup dictionary
      //////////////////////////////////////////////////////////////////////
      if( !in_standard_template ){
        found_template = search_het_group_dictionary( site_info[current_atom].residue_name );
  
        if( found_template ){
          found = find_bonds( (standard_template.size()-1), current_atom );
        }
        else
          unknown_residue.push_back( site_info[current_atom].residue_name );
      }
      
    }
    
    // Try a distance only based check for atoms not found in templates. 
    //////////////////////////////////////////////////////////////////////
    if( !found ){
      //cout << " not in template " << current_atom << " " << site_info[current_atom].atom_name << endl;
      int isolated_site = find_bonds( current_atom );
      
      if( parameters.interactive && 
          isolated_site &&
          !(site_info[current_atom].neighbor_list).size() ){
        stringstream warning;
        warning << "No connections to atom " << current_atom << " were found.";
        molFramework->add_warning( 1, warning.str() );
      }
    }
    
    in_standard_template = found = false;
  }

  return;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
//   The "standard" templates refer to the 20 standard amino acids, and a 
//   few of the most common HETATM residues, such as water. These common 
//   residues are read and stored in an array of template structures.                 
////////////////////////////////////////////////////////////////////////////
void OldCovalentBonder::read_standard_residue_library(){

  string linebuf;

  // 1. For each template, store the connections. A template is stored 
  //    in a structure of type "template_data". There are two elements 
  //    in this structure. The name of the group that is being stored, 
  //    and an object of type "map". The map object is accessed by atom 
  //    name, and it returns a list of the atoms connected to it. See 
  //    manual for more details. 
  //////////////////////////////////////////////////////////////////////
  string std_residue_lib = parameters.path + "/lib/std_dictionary.txt";
  fstream res_library( std_residue_lib.c_str(), ios::in );
  if( !res_library ){
    cout << " ERROR: Standard Residue Library not found." << endl;
    exit(1);
  }
    
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
//   Read and store a table atom-atom distance cutoffs. These cutoffs are an
//   upper bound on the distance when identifying covalent bonds. They are 
//   listed in /lib/distance_cutoffs.txt . 
////////////////////////////////////////////////////////////////////////////
void OldCovalentBonder::read_distance_cutoff_table(){

  string linebuf;
  string bond;
  string temp;
  float cutoff = 0;

  // 1. Open the file with the distance cutoff's listed.
  //////////////////////////////////////////////////////////////////////
  string distance_cutoff_table = parameters.path + "/lib/distance_cutoffs.txt";
  fstream cutoff_table( distance_cutoff_table.c_str(), ios::in );
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

      covalent_bond_table.insert( pair<string,float> (bond, cutoff) );    
    }
  }   

  cutoff_table.close();

  // Uncomment the following lines to print the distance cutoff table.
  //////////////////////////////////////////////////////////////////////
  /*
  distance_cutoff = covalent_bond_table.begin();
  while( distance_cutoff != covalent_bond_table.end() ){
    cout << distance_cutoff->first << " " << distance_cutoff->second << endl;
    distance_cutoff++;
  }
  */
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
void OldCovalentBonder::store_residue_template_data( fstream *file_ptr, string residue_name ){

  unsigned int neighbor_start = 0;
  size_t end_of_line = 0; // unsigned int here breaks on flexweb (gcc 3.4.3 20050227 (Red Hat 3.4.3-22.1))

  string linebuf;
  string this_atom;
  string neighbor;
  list<string> temp_list;
  map<string,list<string> > temp_map;


  while( linebuf.find("END") == string::npos ){
    getline( *file_ptr, linebuf );
    
    if( !linebuf.find("CONECT") ){ // We only need the CONECT data.
      
      this_atom = linebuf.substr(11,4);
      removeWhiteSpacePadding( this_atom );

      neighbor_start = 20; // The atom neighbors start in column 20 of the template file. 
      end_of_line = linebuf.find_first_not_of(" \n", neighbor_start ); // don't store a newline character.
      
      // Read and store the neighbors of "this_atom" until we reach a newline. 
      //////////////////////////////////////////////////////////////////////
      while( end_of_line != string::npos ){ 
        neighbor = linebuf.substr( neighbor_start, 4 );
        neighbor_start += 5;
        end_of_line = linebuf.find_first_not_of(" \n", neighbor_start-1 );
        removeWhiteSpacePadding( neighbor );
        temp_list.push_back( neighbor );
      }
      
      // create a map between an atom's CONECT record and the list of atom's
      //  its connected to. 
      //////////////////////////////////////////////////////////////////////
      temp_map.insert( pair<string,list<string> > (this_atom,temp_list));
      temp_list.clear();
    }
  }

  standard_template.push_back( template_data(residue_name, temp_map) );
}
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
// Description:
//   Open the file het_dictionary.txt in the directory "/lib". Search the 
//   dictionary for a RESIDUE record corresponding to the residue_name passed 
//   to this function as a parameter. If the residue is found, store the 
//   data in a new data structure a push the data onto the standard template 
//   vector list. 
// Parameters:
//   residue_name - Three letter name of a residue 
// Return Value List:
//   true - If the residue is found in the file "het_dictionary.txt".
//   false - if the residue is not found in the file "het_dictionary.txt".
////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::search_het_group_dictionary( string residue_name ){

  string linebuf;

  // 1. Open the HET group dictionary.
  //////////////////////////////////////////////////////////////////////
  string het_residue_lib = parameters.path + "/lib/het_dictionary.txt";
  fstream het_dictionary( het_residue_lib.c_str(), ios::in );
  if( !het_dictionary ){
    cout << "ERROR: Could not open file [het_dictionary.txt] " << endl;
    exit(1);
  }

  // 2. Search the dictionary file for the residue named "residue_name".
  //////////////////////////////////////////////////////////////////////
  while( !het_dictionary.eof() ){
    getline( het_dictionary, linebuf );
    
    if( !linebuf.find("RESIDUE") &&
        linebuf.find(residue_name) != string::npos ){ 
      store_residue_template_data( &het_dictionary, residue_name );
      return(true);
    }
  }
  het_dictionary.close();

  // 3. If the residue name is no found in the dictionary, return 0.
  //////////////////////////////////////////////////////////////////////
  return(false);
}
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::find_bonds( int current_template, unsigned int current_atom ){

  int bars = 0;
  //  float distance = 0.0; // FIXME - warning: unused variable 'distance'
  string atom_name;
  vector<int> grids_to_check;

  // 1. Check to see if the current_atom is listed in the residue_template.
  ////////////////////////////////////////////////////////////////////////////////
  atom_name = site_info[current_atom].atom_name;
  atom_connections = (standard_template[current_template].std_group_map).find( atom_name );

  // Debug code
  //////////////////////////////////////////////////////////////////////
  //map<string,list<string> >::iterator test = standard_template[current_template].std_group_map.begin();
  //while( test != standard_template[current_template].std_group_map.end() ) {
  //cout << "[" << atom_name << "] [" << test->first << "]" << endl;
  //test++;
  //}
    
  if( atom_connections != (standard_template[current_template].std_group_map).end() ){ 

    //cout << " this atom " << current_atom << " is in a template " << endl;
    // 2. Check the current_atom's coordintate_grid, and all 28 neighboring grids for 
    //    atoms  with larger  atom numbers.  Next, compare FIRST_number, name, and 
    //    FIRST_chain_ID, in that order (Connected logical AND's will exit the "if" statement
    //    the first time a comparison fails). For nearby atoms that pass the above tests, 
    //    compare their names to those found in the "atom_connections" list. 
    //////////////////////////////////////////////////////////////////////////////////////////

    // a. Create a list of coordinate grids to check. Check each grid. 
    //////////////////////////////////////////////////////////////////////
    grids_to_check = molFramework->getGrids( current_atom );
    current_grid = grids_to_check.begin();
    
    while( current_grid != grids_to_check.end() ){

      // b. Check every atom in each grid. 
      //////////////////////////////////////////////////////////////////////
      neighbor_atom = molFramework->coordinate_grid[*current_grid].begin();
     
      while( neighbor_atom != molFramework->coordinate_grid[*current_grid].end() ){

        // c. Check the distance between current_atom and neighbor_atom. Only check
        //    atoms with larger site numbers, to prevent double counting. 
        //////////////////////////////////////////////////////////////////////
        //cout << current_atom << " " << *neighbor_atom << endl;
        if( *neighbor_atom > current_atom &&
            meets_atom_atom_distance_cutoff(current_atom, *neighbor_atom) ){

          if(is_known_neighbor(current_atom, *neighbor_atom, atom_connections->second, &bars, current_template) )
            molFramework->add_to_site_info_array( current_atom, *neighbor_atom, bars );
          
          else if( is_peptide_bond( current_atom, *neighbor_atom ) )
            molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 6 );
              
          else if( is_disulfide_bond( current_atom, *neighbor_atom ) ){
            if( parameters.lock_disulfide )
              molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 6 );
            else
              molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
          }
      
          else if( is_internucleic_bond( current_atom, *neighbor_atom ) ) //AJR 11.03.05
            molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5);
      
          // If the atom pair does not meet any of the above specific checks, 
          // but it has met the atom-atom distance check, conect them, but tag 
          // it (or send to user in interactive mode).
          //////////////////////////////////////////////////////////////////////
          else{
      
            // Connect hydrogens by default. Valency errors will be picked up in the
            // validate stucture function.
            //////////////////////////////////////////////////////////////////////
            if( site_info[current_atom].element_name == "H " ||
                site_info[*neighbor_atom].element_name == "H " ){
              //cout << "adding hydrogen covalent bond " << current_atom << " " << *neighbor_atom << endl;
              molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
            }
            else{
              if( parameters.interactive ){
                if( molFramework->userInteraction3(current_atom, *neighbor_atom, site_info, molFramework->getDistance(current_atom,*neighbor_atom) ))
                  molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
              }
              else
                molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
            }
          }
        }
        neighbor_atom++;
      }

      current_grid++;
    }
    
    return(1);
  }

  
  // 3. If the current_atom name was not listed in it's residue template, 
  //    send it back to the identifyCovalentBonds() routine for a 
  //    distance-only based neighbor check.
  //////////////////////////////////////////////////////////////////////
  //cout << "atom not found in template list." << endl;
  return(0);
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is an overloaded (underloaded!?) version of the find_bonds()
//   routine. This version is called when the current_atom does not belong
//   to any residue in the standard or HET group dictionaries.
// Parameters:
//   current_atom - FIRST_number of an atom.
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::find_bonds( unsigned int current_atom ){
  
  int found_bond = 0;
  // float distance = 0.0; // FIXME - warning: unused variable 'distance'
  vector<int> grids_to_check;

  // 1. Check the current_atom's coordintate_grid, and all 28 neighboring grids for 
  //    atoms  with larger  atom numbers.  Next, compare FIRST_number, name, and 
  //    FIRST_chain_ID, in that order (Connected logical AND's will exit the "if" statement
  //    the first time a comparison fails). 
  /////////////////////////////////////////////////////////////////////////////////////  
    
  // a. Create a list of coordinate grids to check. Check each grid. 
  ////////////////////////////////////////////////////////////////////////////////
  grids_to_check = molFramework->getGrids( current_atom );
  current_grid = grids_to_check.begin();

  while( current_grid != grids_to_check.end() ){
    
    // b. Check every atom in each grid. 
    //////////////////////////////////////////////////////////////////////
    neighbor_atom = molFramework->coordinate_grid[*current_grid].begin();
    
    while( neighbor_atom != molFramework->coordinate_grid[*current_grid].end() ){

      // c. Check the distance between current_atom and neighbor_atom. Only check
      //    atoms with larger site numbers, to prevent double counting. 
      ////////////////////////////////////////////////////////////////////////////////
      if( *neighbor_atom > current_atom &&
          meets_atom_atom_distance_cutoff(current_atom, *neighbor_atom) ){
        //cout << "checking " << current_atom << "  " << *neighbor_atom << endl;
        found_bond = true;
        // Connect hydrogens by default. Valency errors will be picked up in the
        // validate stucture function.
        //////////////////////////////////////////////////////////////////////
        if( site_info[current_atom].element_name == "H " ||
            site_info[*neighbor_atom].element_name == "H " )
          molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
        else{
          if( parameters.interactive ){
            if( molFramework->userInteraction3(current_atom, *neighbor_atom, site_info, molFramework->getDistance(current_atom,*neighbor_atom)) )
              molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
          }
          else
            molFramework->add_to_site_info_array( current_atom, *neighbor_atom, 5 );
        }
      }
      
      neighbor_atom++;
    }
  
    current_grid++;
  }
  
  // If no bonds to this atom were found, tell the user. 
  ////////////////////////////////////////////////////////////////////////////////
  if( found_bond )
    return(0);
  else
    return(1);
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Check for the occurrance of a disulfide bridge between two cystine 
//   residues. As a validation, the sulfur should be in the thiol form or 
//   disulfide form, but not both.
// Parameters:
//   atom_1 - FIRST_number of an atom. 
//   atom_2 - FIRST_number of an atom.
// Return Value:
//   true - Meets disulfide bond criteria, a disulfide bond exists between 
//     atom_1 and atom_2.
//   false - Failed criteria. No disulfide bond between atom_1 and atom_2.
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::is_disulfide_bond( unsigned int atom_1, unsigned int atom_2 ){

  if( site_info[atom_1].atom_name == "SG"    &&
      site_info[atom_1].residue_name == "CYS"  &&
      site_info[atom_2].atom_name == "SG"        &&
      site_info[atom_2].residue_name == "CYS"      &&
      meets_atom_atom_distance_cutoff( atom_1, atom_2 ) )
    return( true );
  else
    return( false );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::is_known_neighbor( unsigned int current_atom, unsigned int neighbor, 
             list<string> list_of_neighbors, int *bars,
             int current_template ){

  //cout << site_info[current_atom].seq_number << " " << site_info[neighbor].seq_number << endl;
  
  // Simple first pass check that these atoms have the same residue number,
  // name, and chain ID. 
  ////////////////////////////////////////////////////////////////////////////
  if( site_info[current_atom].seq_number      != site_info[neighbor].seq_number ||
      site_info[current_atom].residue_name    != site_info[neighbor].residue_name ||
      site_info[current_atom].FIRST_chain_ID  != site_info[neighbor].FIRST_chain_ID ){
    return(false);
  }

  // Check to see if the neighbor is listed as a known neighbor in the residue template.
  ////////////////////////////////////////////////////////////////////////////
  known_neighbor = list_of_neighbors.begin();
  while( known_neighbor != list_of_neighbors.end() ){

    if( (*known_neighbor).find( site_info[neighbor].atom_name ) != string::npos ) {

      map<string, list<string> >::iterator test_this;
      test_this = (standard_template[current_template].std_group_map).find( *known_neighbor );
      
      if( is_locked(current_atom, list_of_neighbors.size(), 
          neighbor, (test_this->second).size() ) ) 
        *bars = 6;
      else 
        *bars = 5;
      
      return(true);
    }
    
    known_neighbor++;
  }
  
  return(false);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
//   The peptide bond connection is not listed in the standard residue
//   library templates. Here, we'll check to see if both atoms are the 
//   proper distance apart, and have the right atom names.
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::is_peptide_bond( unsigned int atom_1, unsigned int atom_2 ){

  if( site_info[atom_1].atom_name == "C" &&
      site_info[atom_2].atom_name == "N" &&
      ( (site_info[atom_2].seq_number - site_info[atom_1].seq_number) == 1) &&
      site_info[atom_1].FIRST_chain_ID == site_info[atom_2].FIRST_chain_ID &&
      meets_atom_atom_distance_cutoff(atom_1, atom_2) ){
    //cout << "is peptide bond" << endl;
    return(true);
  }
  else if( site_info[atom_1].atom_name == "N" &&
     site_info[atom_2].atom_name == "C" &&
     ( (site_info[atom_1].seq_number - site_info[atom_2].seq_number) == 1) &&
     site_info[atom_1].FIRST_chain_ID == site_info[atom_2].FIRST_chain_ID &&
     meets_atom_atom_distance_cutoff(atom_1, atom_2) ){
    //cout << "is peptide bond" << endl;
    return(true);
  }
  else
    return(false);

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The internucleic acid bond connection (O3* to P) is not listed in the standard 
//   residue library templates. Here, we'll check to see if both atoms are the 
//   proper distance apart, and have the right atom names.
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::is_internucleic_bond( unsigned int atom_1, unsigned int atom_2 ){

  if( site_info[atom_1].atom_name == "O3*" &&
      site_info[atom_2].atom_name == "P" &&
      ( (site_info[atom_2].seq_number - site_info[atom_1].seq_number) == 1) &&
      site_info[atom_1].FIRST_chain_ID == site_info[atom_2].FIRST_chain_ID &&
      meets_atom_atom_distance_cutoff(atom_1, atom_2) ){
    //cout << "is peptide bond" << endl;
    return(true);
  }
  else if( site_info[atom_1].atom_name == "P" &&
     site_info[atom_2].atom_name == "O3*" &&
     ( (site_info[atom_1].seq_number - site_info[atom_2].seq_number) == 1) &&
     site_info[atom_1].FIRST_chain_ID == site_info[atom_2].FIRST_chain_ID &&
     meets_atom_atom_distance_cutoff(atom_1, atom_2) ){
    //cout << "is peptide bond" << endl;
    return(true);
  }
  else
    return(false);

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
int OldCovalentBonder::meets_atom_atom_distance_cutoff( unsigned int atom_1, unsigned int atom_2 ){

  string this_bond;

  // 1. Create a string listing the two elements in the bond. 
  ////////////////////////////////////////////////////////////////////////////////
  if( site_info[atom_1].element_name > site_info[atom_2].element_name )
    this_bond = site_info[atom_2].element_name + site_info[atom_1].element_name;
  else
    this_bond = site_info[atom_1].element_name + site_info[atom_2].element_name;
  
  // 2. Get the stored distance cutoff for this bond type. 
  ////////////////////////////////////////////////////////////////////////////////
  distance_cutoff = covalent_bond_table.find( this_bond );

  if( distance_cutoff != covalent_bond_table.end() ){
    //cout << distance_cutoff->second << " " << getDistance(atom_1,atom_2) << endl;
    if( molFramework->getDistance(atom_1, atom_2) <= (distance_cutoff->second) ){
      return(1);
    }
    else{
      return(0);
    }
  }
  else{
    cout << " Warning: Bond (" << this_bond << ") not found in covalent bond distance table. " 
         << site_info[atom_2].element_name <<endl;
    return(0);
  }
}
////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////// 
// Description:
////////////////////////////////////////////////////////////////////////////////
bool OldCovalentBonder::is_locked( int site_1, int mult_site_1,
             int site_2, int mult_site_2 ){
  
  if( site_info[site_1].element_name > site_info[site_2].element_name ){
    swap( site_1, site_2 );
    swap( mult_site_1, mult_site_2 );
  }

  if( site_info[site_1].element_name == "C " &&
      mult_site_1 <= 3 &&
      site_info[site_2].element_name == "C " &&
      mult_site_2 <= 3 ){
    return( true );
  }
  else if( site_info[site_1].element_name == "C " && mult_site_1 <= 3 &&
     site_info[site_2].element_name == "N " && mult_site_2 <= 3 ){
    //cout << "locking " << site_1 << " " << site_2 << endl;
    return( true );
  }
  
  return(false);
}
////////////////////////////////////////////////////////////////////////////////

