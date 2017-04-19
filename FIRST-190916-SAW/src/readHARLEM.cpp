#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "readHARLEM.h"

extern Parameters parameters;
//vector<int> seq_numbers;

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Input files created by the software HARLEM contain both the structure information,
//   in PDB file format, but additional information mapping the atom numbers to 
//   HARLEM's "atom name" identifies that are used internally by that program to 
//   uniquely identify every atom. This function will map FIRST's atom numbers to
//   HARLEM's string atom id's. 
////////////////////////////////////////////////////////////////////////////////
void read_harlem_file( MolFramework &structure ){

  ////////////////////////////////////////////////////////////////////////////////
  int total_sites = scan_harlem_file( structure );
   
  ////////////////////////////////////////////////////////////////////////////////
  structure.initialize( total_sites );

  ////////////////////////////////////////////////////////////////////////////////
  read_harlem_data( structure );

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
int scan_harlem_file( MolFramework &structure ){

  bool found_atom_start = false;
  bool found_atom_end   = false;
  int  total_sites = 0;  
  int  counter = 0;
  int  harlem_atom_count = 0;
  string linebuf;
  string test;

  // Open the harlem file
  //////////////////////////////////////////////////////////////////////
  ifstream harlem_file( structure.infile_name.c_str(), ios::in );
  if( !harlem_file ){
    cout << " Error: Could not open file." << endl;
    exit(1);
  }

  while( !harlem_file.eof() ){

    getline( harlem_file, linebuf );

    if( linebuf.find("#ATOMS") != string::npos ){
      found_atom_start = true;
      do{
	getline( harlem_file, linebuf );
	if( linebuf.find("#END ATOMS") != string::npos )
	  found_atom_end = true;
	else
	  total_sites++;
      } while( !found_atom_end );     
    }

    if( linebuf.find("ALLATOMLIST") != string::npos ){
      size_t start = 0;
      size_t end   = 0;
      string harlem_name;

      // Skip the first field
      //////////////////////////////////////////////////////////////////////
      start = linebuf.find_first_not_of(" \t\n");
      end   = linebuf.find_first_of(" \t\n", start);

      // The second field is the number of atom name records. Store this for 
      // error checking.
      //////////////////////////////////////////////////////////////////////
      start = linebuf.find_first_not_of(" \t\n", end);
      end   = linebuf.find_first_of(" \t\n", start);
      harlem_atom_count = atoi( linebuf.substr(start,end-start).c_str() );

      do{
	if( counter != 0 ){
	  getline( harlem_file, linebuf);
	  start = end = 0;
	}

	start = linebuf.find_first_not_of(" \t\n#", end);
	end   = linebuf.find_first_of(" \t\n#", start);
	while( end != string::npos ){
	  harlem_name = linebuf.substr(start, end-start);
	  counter++;
	  //cout << counter << " " << harlem_name << endl;
	  structure.harlem_names.insert( pair<int,string> (counter, harlem_name) );

	  start = linebuf.find_first_not_of(" \t\n#", end);
	  end   = linebuf.find_first_of(" \t\n#", start);
	}

      } while( linebuf.find("###") != string::npos );

    }

  }

  if( harlem_atom_count != total_sites ){
    cout << "WARNING: Total #ATOMS records does not equal total atoms in ALLATOMLIST." << endl;
    //    exit(1);
  }

  if( harlem_atom_count != counter ){
    cout << "WARNING: Total HARLEM atom names read does not equal the total stated by ALLATOMLIST." << endl;
    //    exit(1);
  }

  return( total_sites );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
string get_harlem_element_name( string &atom_name ){
  
  string element = "  ";
  element.at(0) = atom_name.at(0);
  
  return( element );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void read_harlem_data( MolFramework &structure ){

	// bool found_atom_start = false; // FIXME - warning: unused variable 'found_atom_start'
  bool found_atom_end   = false;
  int counter = 0;
  int start = 0;
  int end = 0;
  string linebuf;
  string this_field;
  string temp_string = "    ";
  int old_res_num = 0;
  structure.total_residues = 0;
 
  unsigned int checkChainCounter = 1;
  unsigned int currentOrigChainID = 0;
  int currentFIRSTChainID = 0;

  //////////////////////////////////////////////////////////////////////
  ifstream harlem_file( structure.infile_name.c_str(), ios::in );
  if( !harlem_file ){
    cout << " Error: Could not open file." << endl;
    exit(1);
  }

  //////////////////////////////////////////////////////////////////////
  while( !harlem_file.eof() ){
    
    getline( harlem_file, linebuf );
    if( linebuf.find("#ATOMS") != string::npos ){

      do{
	getline( harlem_file, linebuf );
	if( linebuf.find("#END ATOMS") != string::npos )
	  found_atom_end = true;
	else{
	  counter++;

	  // Unique integer for use internally by FIRST.
	  //////////////////////////////////////////////////////////////////////
	  structure.site_info[counter].FIRST_number = counter;

	  // Store the record name as ATOM for converting to PDB format.
	  //////////////////////////////////////////////////////////////////////
	  structure.site_info[counter].record_name = "ATOM  ";

	  // Atom serial number in orginal file.
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_not_of(" \t\n",0);
	  end   = linebuf.find_first_of(" \n\t", start);
	  this_field = linebuf.substr(start, end-start);
	  structure.site_info[counter].orig_atom_number = this_field.c_str();

	  // Skip the next field
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_not_of(" \t\n",end);
	  end   = linebuf.find_first_of(" \n\t", start);
	  
	  // Map the origianl atom number to the number used by FIRST. They will
	  // only be different if the atom numbers in the input PDB file did not
	  // start with "1". 
	  //////////////////////////////////////////////////////////////////////

	  structure.orig_2_FIRST_atom_number.insert( pair<int,int> ( atoi( this_field.c_str() ), 
								    structure.site_info[counter].FIRST_number) );

	  // Atom name.
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_not_of(" \t\n\"",end);
	  end   = linebuf.find_first_of(" \n\t\"", start);
	  structure.site_info[counter].atom_name = linebuf.substr(start, end-start);

	  //////////////////////////////////////////////////////////////////////
	  // BMH 1.22.04 temporary code. need to get rid of column setting on std and het library templates. 
	  /*
	  this_field = linebuf.substr(start, end-start);

	  if( this_field.size() == 1 )
	    temp_string = " " + this_field + "  ";
	  else if( this_field.size() == 2 )
	    temp_string = " " + this_field + " ";
	  else if( this_field.size() == 3 ) 
	    temp_string = " " + this_field;
	  else
	    temp_string = this_field;
	  */
	  //cout << "(" << this_field << ") (" << temp_string << ")" << endl;
	  //structure.site_info[counter].atom_name = temp_string;

	  // Element name.
	  //////////////////////////////////////////////////////////////////////
	  structure.site_info[counter].element_name = structure.getElementName( structure.site_info[counter].atom_name );
	  	  
	  // X-coordinate in Angstroms.
	  start = linebuf.find_first_not_of(" \t\n\"",end);
	  end   = linebuf.find_first_of(" \n\t", start);
	  this_field = linebuf.substr(start, end-start);
	  //cout << "(X " << this_field << ") ";
	  structure.site_info[counter].coords[X] = atof( this_field.c_str() );
	  //  Compare to min/max in X direction. 
	  if( structure.site_info[counter].coords[X] < structure.min_coords[X] )
	    structure.min_coords[X] = structure.site_info[counter].coords[X];
	  if( structure.site_info[counter].coords[X] > structure.max_coords[X] )
	    structure.max_coords[X] = structure.site_info[counter].coords[X];
	  
	  // Y-coordinate in Angstroms. 
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_not_of(" \t\n",end);
	  end   = linebuf.find_first_of(" \n\t", start);
	  this_field = linebuf.substr(start, end-start);
	  //cout << "(Y " << this_field << ") ";
	  structure.site_info[counter].coords[Y] = atof( this_field.c_str() );
	  //  Compare to min/max in Y direction. 
	  if( structure.site_info[counter].coords[Y] < structure.min_coords[Y] )
	    structure.min_coords[Y] = structure.site_info[counter].coords[Y];
	  if( structure.site_info[counter].coords[Y] > structure.max_coords[Y] )
	    structure.max_coords[Y] = structure.site_info[counter].coords[Y];
	  
	  // Z-coordinate in Angstroms.  
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_not_of(" \t\n",end);
	  end   = linebuf.find_first_of(" \n\t", start);
	  this_field = linebuf.substr(start, end-start);
	  //cout << "(Z " << this_field << ") ";
	  structure.site_info[counter].coords[Z] = atof( this_field.c_str() );
	  //  Compare to min/max in Z direction. 
	  if( structure.site_info[counter].coords[Z] < structure.min_coords[Z] )
	    structure.min_coords[Z] = structure.site_info[counter].coords[Z];
	  if( structure.site_info[counter].coords[Z] > structure.max_coords[Z] )
	    structure.max_coords[Z] = structure.site_info[counter].coords[Z];

	  // Chain identifier in original file;
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_of("\"",end);
	  end   = linebuf.find_first_of(" \n\t\"", start+1);
	  this_field = linebuf.substr(start+1, end-start-1); 
	  //cout << "(C " << this_field << ") " << start << " " << end;
	  structure.site_info[counter].chain_ID = linebuf[start+1];
	  //cout << "new chain ID = " << structure.site_info[counter].chain_ID << endl;
	  
	  // Chain identifier used internally by FIRST.
	  //////////////////////////////////////////////////////////////////////
	    
	  if( structure.site_info[counter].chain_ID == currentOrigChainID &&
	      structure.site_info[counter].FIRST_chain_ID == currentFIRSTChainID ){
	    structure.site_info[counter].FIRST_chain_ID = checkChainCounter;
	  }
	  else{
	    checkChainCounter++;
	    currentOrigChainID = structure.site_info[counter].chain_ID;
	    currentFIRSTChainID = structure.site_info[counter].FIRST_chain_ID;
	    
	    structure.site_info[counter].FIRST_chain_ID = checkChainCounter;    
	  }
	  


	    //	  cout << "new chain ID = " << structure.site_info[counter].chain_ID <<" "<<  structure.site_info[counter].FIRST_chain_ID<< endl;

	  // Residue name. 
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_not_of(" \t\n\"",end);
	  end   = linebuf.find_first_of("# \n\t\"", start);
	  this_field = linebuf.substr(start, end-start);
	  //cout << "(R " << this_field << ") ";
	  structure.site_info[counter].residue_name = this_field;
	  
	  // Map the original chain ID to a unique integer chain ID used by FIRST
	  //////////////////////////////////////////////////////////////////////
	  structure.orig_2_FIRST_chain_ID.insert( pair<int,int> (structure.site_info[counter].chain_ID, 
								 structure.site_info[counter].FIRST_chain_ID) );
	  
	  // Residue sequence number.
	  //////////////////////////////////////////////////////////////////////
	  start = linebuf.find_first_of(" ",end);
	  start = linebuf.find_first_not_of(" \t\n",start);
	  end   = linebuf.find_first_of(" \n\t", start);
	  this_field = linebuf.substr(start, end-start);
	  //cout << "(S " << this_field << ") ";
	  structure.site_info[counter].seq_number = atoi( this_field.c_str() );
	  //seq_numbers.push_back( structure.site_info[counter].seq_number );
          if(structure.site_info[counter].seq_number != old_res_num) // Counter of total residues in the Harlem
	    {
              old_res_num = structure.site_info[counter].seq_number ;
              structure.total_residues++;
	    }
	}

      } while( !found_atom_end );     
    }
  }

}
////////////////////////////////////////////////////////////////////////////////
