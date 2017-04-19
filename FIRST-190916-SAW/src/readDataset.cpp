#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "readPDB.h"
#include "readDataset.h"

extern const Parameters parameters;
map<int,int> tether_map;

////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////
void readDatasetFile( MolFramework &structure ){

  // Read through pdb file once to count the number of ATOM, HETATM, and 
  // CONECT records. Used for proper array size declarations when constructing
  // a pdb_object. If we find a MODEL record in the file, call the appropriate
  // functions to handle NMR models. 
  ////////////////////////////////////////////////////////////////////////////////
  int total_sites = scanPDB_File( structure );
   
  // Allocate the site_info array in the pdb_object.
  ////////////////////////////////////////////////////////////////////////////////
  structure.initialize( total_sites );

  // Read in the data from the PDB file. Call the overloaded read_pdb() 
  // function if we're reading in a model from a NMR structure.
  ////////////////////////////////////////////////////////////////////////////////
  if( structure.is_nmr ){
    readPDB_Data( structure, structure.model_number );
  }
  else{
    readPDB_Data( structure );
  }

  // Read in the central-force and torsional-force constraints from old
  // *_FIRSTdataset file. 
  ////////////////////////////////////////////////////////////////////////////////
  readCF_AndTF_Constraints( structure );

  // Read hydrogen bond and hydrophobic tether constraints.
  ////////////////////////////////////////////////////////////////////////////////
  readHB_AndPH_Constraints( structure );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   When using the older version of the FIRST input format, with the extension 
//   _FIRSTdataset, read in the central-force (REMARK:CF) and torsional-force 
//   (REMARK:TF) constraints, and store them in the site_info array. 
////////////////////////////////////////////////////////////////////////////////
void readCF_AndTF_Constraints( MolFramework &structure ){

  unsigned int atom_1 = 0;
  unsigned int atom_2 = 0;

  string current_line;

  // 1. Open the dataset file. 
  //////////////////////////////////////////////////////////////////////
  ifstream dataset_file( structure.infile_name.c_str(), ios::in );
  if( !dataset_file ){
    cout << "Error: Could not open dataset file. [in: dataset.cpp-->read_CF_and_TF_constraints]" << endl;
    exit(1);
  }

  while( !dataset_file.eof() ){
    getline( dataset_file, current_line );

    if( current_line.find("REMARK:CF") != string::npos ){
      atom_1 = atoi( (current_line.substr( 9, 6)).c_str() );
      atom_2 = atoi( (current_line.substr(15, 6)).c_str() );
      
      if( atom_1 <= structure.total_sites && atom_2 <= structure.total_sites )
	structure.add_to_site_info_array( atom_1, atom_2, 5 );
      else if( atom_1 <= structure.total_sites && atom_2 > structure.total_sites )
	tether_map.insert( pair<int,int> (atom_2, atom_1) );
    }

    if( current_line.find("REMARK:TF") != string::npos ){
      atom_1 = string_2_int(current_line.substr( 9, 6) );
      atom_2 = string_2_int(current_line.substr(15, 6) );
      
      if( atom_1 <= structure.total_sites && atom_2 <= structure.total_sites ){
	if( find( structure.site_info[atom_1].neighbor_list.begin(), 
		  structure.site_info[atom_1].neighbor_list.end(), 
		  atom_2) != structure.site_info[atom_1].neighbor_list.end() ) {
	  structure.site_info[atom_1].number_of_bars[atom_2] = 6;
	  structure.site_info[atom_2].number_of_bars[atom_1] = 6;
	}
	else{
	  cout << " Warning: Found torsional-force constraint without corresponding central-force constraint." << endl;
	  cout << "          This constaint will be placed in the network with 6 bars." << endl;
	  cout << "          " << atom_1 << " - " << atom_2 << endl;
	  structure.add_to_site_info_array( atom_1, atom_2, 6 );
	}
      }
    }
  }

  dataset_file.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void readHB_AndPH_Constraints( MolFramework &structure ){

  int atom_1 = 0;
  int atom_2 = 0;
  int atom_3 = 0;
  int PH_partner = 0;
  
  float energy = 0.0;
  
  string current_line;

  // 1. Open the dataset file. 
  //////////////////////////////////////////////////////////////////////
  ifstream dataset_file( structure.infile_name.c_str(), ios::in );
  if( !dataset_file ){
    cout << "Error: Could not open dataset file. [in: dataset.cpp-->read_CF_and_TF_constraints]" << endl;
    exit(1);
  }
  
  while( !dataset_file.eof() ){
    getline( dataset_file, current_line );

    if( current_line.find("REMARK:HB") != string::npos ){

      atom_1 = string_2_int( current_line.substr(27, 8) );
      atom_2 = string_2_int( current_line.substr(35, 8) );
      atom_3 = string_2_int( current_line.substr(43, 8) );
      
      // If it's a hydrogen bond remark...
      if( current_line.find("PH hydr phob") == string::npos ){
	energy = string_2_float( current_line.substr(16, 11) );
	structure.hydrogen_bonds.push_back( new_bonds(atom_2, atom_3, parameters.hb_bars, energy) );
      }
      // Else it's a hydrophobic tether...
      else{
	PH_partner = tether_map[atom_1 -1];
	structure.hydrophobic_tethers.push_back( new_bonds(PH_partner, atom_3, parameters.hp_bars) );
      }


    }
  }

  tether_map.clear();
  dataset_file.close();
}
////////////////////////////////////////////////////////////////////////////////
