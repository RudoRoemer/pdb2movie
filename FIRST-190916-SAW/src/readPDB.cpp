#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "readPDB.h"
#include "hybrid_36_c.h"

int model_number = 0;
int chain_counter = 0;

int total_sites = 0;
int conect_records = 0;

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Read a PDB file into the current MolFramework.
////////////////////////////////////////////////////////////////////////////////
void readPDB_File( MolFramework &structure ){

  if( parameters.verbose )
    cout << "Reading PDB data..." << endl;
	
  total_sites = 0;
  conect_records = 0;
  model_number = 0;
  chain_counter = 0;
  
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
 
  // Check chain id's
  // The the function readATOM_OrHETATM_Line, the internal chain counter, 
  // FIRST_chain_ID, is only incremented on the occurence of a TER record. 
  // To assist with files that may be missing TER records, the following
  // code compares both the original and FIRST chain ID's to a running 
  // count, reassigning FIRST_chain_id's as necessary. 
  //
  // BMH TODO 2006.11.17 - move into function
  //////////////////////////////////////////////////////////////////////
  unsigned int checkChainCounter = 1;
  unsigned int currentOrigChainID = 0;
  int currentFIRSTChainID = 0;
  for( int a = 1; a <= total_sites; a++ ){
    
    // set initial data
    if( a == 1 ){
      currentOrigChainID = structure.site_info[a].chain_ID;
      currentFIRSTChainID = structure.site_info[a].FIRST_chain_ID;
    }
    
    if( structure.site_info[a].chain_ID == currentOrigChainID &&
	structure.site_info[a].FIRST_chain_ID == currentFIRSTChainID ){
      structure.site_info[a].FIRST_chain_ID = checkChainCounter;
    }
    else{
      checkChainCounter++;
      currentOrigChainID = structure.site_info[a].chain_ID;
      currentFIRSTChainID = structure.site_info[a].FIRST_chain_ID;

      structure.site_info[a].FIRST_chain_ID = checkChainCounter;    
    }
  }

  //for( int a = 1; a <= total_sites; a++ ){
    // structure.site_info[a].print();
  //}
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Compare a the model number from a MODEL record in a pdb file to the 
//   model number to be analyzed. If they are the same, return true, indicating
//   the start of the ATOM data corresponding to the NMR model we want to
//   analyze. 
// Parameters:
//   current_line - A string containing a line from the input pdb file.
// Return Type:
//   true - If the model number listed in the string current_line is the model 
//          number we are looking for.
//   false - Otherwise.
////////////////////////////////////////////////////////////////////////////////
bool isCorrectModel( string current_line, int model_number ){
  
  string temp_string;
  temp_string.assign( current_line, 11, 14 );
  
  int this_model = atoi( temp_string.c_str() );

  //cout << "model check " << model_number << " " << this_model << endl;
  if( this_model == model_number )
    return(true);
  else
    return(false);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Here we're reading through the pdb formatted file and counting the 
//   number of ATOM, HETATM, and CONECT records. No data are read in, we 
//   just need to know how many of each record type there are for proper 
//   memory allocation when constructing class objects of the type "molecule".
////////////////////////////////////////////////////////////////////////////////
int scanPDB_FileWithoutPrompting( MolFramework &structure ){
  // TODO - rewrite / replace with another method to initialize a vector<*MolFramework> with the correct number of sites for each model
  
  if( parameters.verbose ) {
    cout << "Scanning PDB file..." << endl;
  }
  
  int hetatm_records = 0;
  int atom_records = 0;
  
  int max_hetatm_records = 0;
  int max_atom_records = 0;
  int max_conect_records = 0;
  
  string current_line;
  
  // 1. Open the pdb or dataset file. 
  //////////////////////////////////////////////////////////////////////
  ifstream pdb_file( structure.infile_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << " Error: Could not open file." << endl;
    exit(1);
  }
  
  // 2. Read through the pdb file once to set the necessary array sizes 
  //    on the class constructor. Reset the file flags. 
  //////////////////////////////////////////////////////////////////////
  while( !pdb_file.eof() ) {
    getline( pdb_file, current_line );
    
    //cout << current_line << endl;

    if( current_line.find("ATOM") == 0 ){
      atom_records++;
      
    } 
    else if( current_line.find("HETATM") == 0 ){
      if( structure.using_dataset ){
	if( current_line.find("BMH") == string::npos )
	  hetatm_records++;
      } 
      else {
	hetatm_records++;
      }
      
    } 
    else if( !current_line.find("CONECT") ) {
      conect_records++;
      
    } 
    else if (!current_line.find("MODEL")) { // reset the current counts for each model
      atom_records = 0;
      hetatm_records = 0;
      conect_records = 0;
      
    } 
    else if( current_line.find("ENDMDL") != string::npos ){
//      cout << "test " << !current_line.find("ENDMDL") << endl;
      // else if( !current_line.find("ENDMDL") ) { // ensure that the max_* values reflect the maximum for each model (if any)
//      cout << "here" << endl;
      if (max_atom_records < atom_records) {
	max_atom_records = atom_records;
      }
      
      if (max_hetatm_records < hetatm_records) {
	max_hetatm_records = hetatm_records;
      }
      
      if (max_conect_records < conect_records) {
	max_conect_records = conect_records;
      }
    }
  }
  
  // ensure that the max_* values reflect the maximum 
  if (max_atom_records < atom_records) {
    max_atom_records = atom_records;
  }
  
  if (max_hetatm_records < hetatm_records) {
    max_hetatm_records = hetatm_records;
  }
  
  if (max_conect_records < conect_records) {
    max_conect_records = conect_records;
  }
  
  pdb_file.close();
  
  atom_records = max_atom_records;
  hetatm_records = max_hetatm_records;
  conect_records = max_conect_records;
  
  if( !atom_records && !hetatm_records ){
    cout << "Error: No ATOM or HETATM records were found in the input file." << endl;
    exit(1);
  }
  
  return( atom_records + hetatm_records );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Here we're reading through the pdb formatted file and counting the 
//   number of ATOM, HETATM, and CONECT records. No data are read in, we 
//   just need to know how many of each record type there are for proper 
//   memory allocation when constructing class objects of the type "molecule".
////////////////////////////////////////////////////////////////////////////////
int scanPDB_File( MolFramework &structure ){
  
  if( parameters.verbose )
    cout << "Scanning PDB file..." << endl;
  
  int hetatm_records = 0;
  int atom_records = 0;
  
  string current_line;
  bool start_reading_file = false;
  
  // 1. Open the pdb or dataset file. 
  //////////////////////////////////////////////////////////////////////
  ifstream pdb_file( structure.infile_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << " Error: Could not open file." << endl;
    exit(1);
  }
  
  // 2. Read through the pdb file once to set the necessary array sizes 
  //    on the class constructor. Reset the file flags. 
  //////////////////////////////////////////////////////////////////////
  while( !pdb_file.eof() ){
    getline( pdb_file, current_line);
    
    if( structure.is_nmr &&
	( current_line.find("MODEL") == 0 ) ){
      if( isCorrectModel(current_line, structure.model_number) )
	start_reading_file = true;
      else
	start_reading_file = false;
    }
    
    else if( !structure.is_nmr )
      start_reading_file = true;
    
    if( start_reading_file ){
      
      while( (current_line.find("ENDMDL") != 0) &&
	     !pdb_file.eof() ){

	if( current_line.find("ATOM") == 0 ){
	  atom_records++;
	}
	
	else if( current_line.find("HETATM") == 0 ){
	  if( structure.using_dataset ){
	    if( current_line.find("BMH") == string::npos )
	      hetatm_records++;
	  }
	  else
	    hetatm_records++;
	}
	
	else if( !current_line.find("MODEL") &&
		 !structure.is_nmr ){
	  structure.model_number = getModelNumber( structure.infile_name );
	  pdb_file.clear();
	  pdb_file.seekg( 0, ios::beg );
	  structure.is_nmr = true;
	  start_reading_file = false;
	  break;
	}
	
	else if( !current_line.find("CONECT") )
	  conect_records++;
	
	getline( pdb_file, current_line);
	
      }

    }
  }
  
  pdb_file.close();
  
  if( !atom_records && !hetatm_records ){
    cout << "Error: No ATOM or HETATM records were found in the input file." << endl;
    exit(1);
  }

  return( atom_records +hetatm_records );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The integer refering to the NMR model number should be in columns 11-14 
//   according the PDB file format V2.2. Identify and return the number of models 
//   stored in the pdb file.
// Return Type:
//   int - Number of the NMR model we are going to analyze.
////////////////////////////////////////////////////////////////////////////////
int countPDB_Models( string infile_name ){
  
  // 1. Open the pdb or dataset file. 
  ////////////////////////////////////////////////////////////////////////////////
  ifstream pdb_file( infile_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << "Error: Could not open file." << endl;
    exit(1);
  }
  
  // 2. Store a list of all the MODEL record lines found in the pdb file.
  //    Prompt the user to select which model they want to use. 
  ////////////////////////////////////////////////////////////////////////////////
  string current_line;
  vector<string> models;
  
  while( !pdb_file.eof() ){
    getline( pdb_file, current_line);
    
    if( !current_line.find("MODEL  ") ){
      models.push_back( current_line );
    }
  }
  pdb_file.close();
  
  return(models.size());
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The integer refering to the NMR model number should be in columns 11-14 
//   according the PDB file format V2.2. Store all the model numbers listed
//   in the input file. Pass the list to the function userInteraction4()
//   in which the user will pick one model to analyze. Return this model
//   number to the calling function. 
// Return Type:
//   int - Number of the NMR model we are going to analyze.
////////////////////////////////////////////////////////////////////////////////
int listModels( string infile_name ){

  // Open the pdb or dataset file. 
  ////////////////////////////////////////////////////////////////////////////////
  ifstream pdb_file( infile_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << "Error: Could not open file." << endl;
    exit(1);
  }
  
  // Store a list of all the MODEL record lines found in the pdb file.
  // Prompt the user to select which model they want to use. 
  ////////////////////////////////////////////////////////////////////////////////
  string current_line;
  vector<string> models;
  
  while( !pdb_file.eof() ){
    getline( pdb_file, current_line);
    
    if( !current_line.find("MODEL  ") ){
      models.push_back( current_line );
    }
  }
  pdb_file.close();
  
  return( userInteraction4(models) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Choose which NMR model to analyze. In interactive mode, the user is
//   given a list of models to choose from. If the user identified which 
//   model to use as a command line option (-n flag), no list is presented,
//   even in interactive mode. In non-interactive mode, the first model 
//   found in the file is used if no model was selected with the -n flag. 
// Return Type:
//   int -  Number of the NMR model we are going to analyze.
////////////////////////////////////////////////////////////////////////////////
int getModelNumber( string infile_name ){
  
  if( parameters.interactive )
    return( listModels(infile_name) );
  else
    return( 1 );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Store the data from the pdb file. No error checking is done at this point, 
//   only data assignments. The structural information is stored in a data 
//   structure of type "pdb_data" The CONECT records are stored in a separate 
//   container for later use in defining the network connectivity. 
////////////////////////////////////////////////////////////////////////////////
void readPDB_Data( MolFramework &structure ){
	
  int counter = 0;
  int rescount =0;
  string current_line;
  
  // 1. Open the pdb or dataset file. 
  ////////////////////////////////////////////////////////////////////////////////
  ifstream pdb_file( structure.infile_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << "Error: Could not open file. [in process_file.cpp--> void read_structural_data_pdb_format() ]" << endl;
    exit(1);
  }
  
  // 2. Read each line of the pdb file, and store necessary information. 
  ////////////////////////////////////////////////////////////////////////////////
  while( !pdb_file.eof() ){
    getline( pdb_file, current_line);
    
    if( !current_line.find("ATOM") ||
	!current_line.find("HETATM") ){
      counter++;
      readATOM_OrHETATM_Line( current_line, structure, counter, rescount);
    }
    else if( !current_line.find("TER") ){
      structure.chainTermini.push_back( counter );
      chain_counter++;
    }
    else if( !current_line.find("CONECT") )
      readCONECT_Record( current_line, structure );

    
  }

  structure.total_residues = structure.unique_res_id.size();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This reads one model from the PDB file is being read in to a MolFramework
// Parameters:
//   structure - reference to MolFramework where the model will be stored
//   pdb_file  - reference to the ifstream positioned at the begining of the model
////////////////////////////////////////////////////////////////////////////////
void readPDB_Model(MolFramework &structure, ifstream &pdb_file) {

  int counter = 0;
  int rescount =0;
  
  chain_counter = 0; // global variable :-(
  string current_line;

  for (getline(pdb_file, current_line); 
       (current_line.find("ENDMDL") && !pdb_file.eof()); 
       getline(pdb_file, current_line)) {
    
    if( !current_line.find("ATOM") ||
        !current_line.find("HETATM") ){
      counter++;
      
      //      readATOM_OrHETATM_Line(current_line, structure, counter);
      readATOM_OrHETATM_Line( current_line, structure, counter, rescount);
    }
    else if( !current_line.find("TER") ){
      structure.chainTermini.push_back( counter );
      chain_counter++;
    }
  }

  structure.total_residues = structure.unique_res_id.size();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is an overloaded version of readPDB_Data(). When the data from the 
//   PDB file is being read in, each model is stored as a separate MolFramework
//   pointer in a vector, structures.
// Parameters:
//   pdb_file_names -  vector<string> of filenames 
////////////////////////////////////////////////////////////////////////////////
void readPDB_Data( vector<MolFramework*> &structures,  vector<string> pdb_file_names){
  
  // Read file_names into structures and store in vectorOfStructures
  for ( vector<string>::iterator filenameIterator=pdb_file_names.begin(); 
	filenameIterator != pdb_file_names.end(); 
	filenameIterator++) {

    string pdb_file_name = (*filenameIterator);
    if (pdb_file_name.length() > 0) {
      readPDB_Data(structures, pdb_file_name); // TODO - extend this to include loading from non pdb files
    }
  }	
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is an overloaded version of readPDB_Data(). When the data from the 
//   PDB file is being read in, each model is stored as a separate MolFramework
//   pointer in a vector, structures.
// Parameters:
//   pdb_file_name - Name of the pdb file
////////////////////////////////////////////////////////////////////////////////
void readPDB_Data( vector<MolFramework*> &structures, string pdb_file_name ){
  // if pdb_file_name has extension '.nam' then it's a list of pdb files; otherwise, it's a single file
  /*  string file_extension;

  int begin_extension = pdb_file_name.find_last_of(".");
	if( begin_extension != string::npos ) {
		file_extension = pdb_file_name.substr(begin_extension);
	}
  
	else {
		file_extension = "";
	}
	
  toLower( file_extension );

  if (file_extension == ".nam") {
	  // TODO - load each line of the .nam file into a vector and call 
    vector<string> pdb_file_names;
    
    ifstream nameFile(pdb_file_name.c_str(),  ios::in);
    
    while (!nameFile.eof()) {
      string linebuf;
      
      getline(nameFile, linebuf);

      pdb_file_names.push_back(linebuf);
    }
    
    nameFile.close();
    
    readPDB_Data(structures, pdb_file_names);
    
    return;
  }*/
  
  string pdb_file_prefix = pdb_file_name;
  int name_length = pdb_file_name.size();
  pdb_file_prefix.erase( name_length-4, 4 ) ;
	
  MolFramework *current_structure = new MolFramework();
  current_structure->infile_name = pdb_file_name;
  current_structure->base_name = pdb_file_prefix;
  
  int total_sites = scanPDB_FileWithoutPrompting( *current_structure );
  current_structure->initialize( total_sites );
  
  // Open the pdb or dataset file. 
  ////////////////////////////////////////////////////////////////////////////////
  ifstream pdb_file( pdb_file_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << "Error: Could not open file. [in process_file.cpp--> void read_structural_data_pdb_format() ]" << endl;
    exit(1);
  }
  
  // Read each line of the pdb file, and store necessary information. 
  ////////////////////////////////////////////////////////////////////////////////
  string current_line;
  
  while( !pdb_file.eof() ){
    getline( pdb_file, current_line);
    
    if( !current_line.find("ATOM") ){
      // TODO - move pdb_file back one line (otherwise will miss one line)
      readPDB_Model(*current_structure, pdb_file);

		current_structure->using_pdb = true;
		//exclude_sites(*current_structure);
		/*set_vdw_radii(*current_structure);
		build_framework(*current_structure);*/
// TODO - verify that we don't need this to build the bond network		analyze_flexibility_and_output(*current_structure);
		
		/*
		 read_data( structure, parameters.infile_name );
		 exclude_sites( structure );
		 set_vdw_radii( structure );
		 build_framework( structure );
		 analyze_flexibility_and_output( structure );*/
      structures.push_back(current_structure);
      
      current_structure = new MolFramework();
      current_structure->infile_name = pdb_file_name;
	  current_structure->base_name = pdb_file_prefix;
      current_structure->initialize( total_sites );
      
    }
    else if( !current_line.find("MODEL") ) {
      readPDB_Model(*current_structure, pdb_file);
      
		current_structure->using_pdb = true;/*
		exclude_sites(*current_structure);
		set_vdw_radii(*current_structure);
		build_framework(*current_structure);*/ 
		// TODO - verify that we don't need this to build the bond network		analyze_flexibility_and_output(*current_structure);
		
      structures.push_back(current_structure);
      
      current_structure = new MolFramework();
      current_structure->infile_name = pdb_file_name;
      current_structure->initialize( total_sites );
      
    } 
    else if( !current_line.find("CONECT") )
      readCONECT_Record( current_line, *current_structure );
  }    
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is an overloaded version of readPDB_Data(). Here, the parameter
//   model_number has been passed. When the data from the PDB file is being read
//   in, only those data corresponding to the given model number are stored. 
// Parameters:
//   model_number - Number of the NMR model to be analyzed.
////////////////////////////////////////////////////////////////////////////////
void readPDB_Data( MolFramework &structure, int model_number ){
	
  int counter = 0;
  int rescount = 0;
  
  // Open the pdb or dataset file. 
  ////////////////////////////////////////////////////////////////////////////////
  ifstream pdb_file( structure.infile_name.c_str(), ios::in );
  if( !pdb_file ){
    cout << "Error: Could not open file. [in process_file.cpp--> void read_structural_data_pdb_format() ]" << endl;
    exit(1);
  }
  
  // Read each line of the pdb file, and store necessary information. 
  ////////////////////////////////////////////////////////////////////////////////
  string current_line;
  string temp_string;
  
  while( !pdb_file.eof() ){
    getline( pdb_file, current_line);
    
    if( current_line.find("MODEL") == 0 ){
      temp_string.assign( current_line, 11, 14 );
      int this_model = atoi( temp_string.c_str() );
      
      if( this_model == model_number ){
				
	while( current_line.find("ENDMDL") ){
	  getline( pdb_file, current_line );
	  
	  if( !current_line.find("ATOM") ||
	      !current_line.find("HETATM") ){
	    counter++;
	    //	    readATOM_OrHETATM_Line( current_line, structure, counter );
	    readATOM_OrHETATM_Line( current_line, structure, counter, rescount);
	  }
	}
      }
    }
    
    else if( !current_line.find("CONECT") )
      readCONECT_Record( current_line, structure );
  }


  structure.total_residues = structure.unique_res_id.size();
  pdb_file.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Store the data from each field defined in a pdb ATOM or HETATM record. Field 
//   boundaries correspond the PDB FILE FORMAT, VERSION 2.2.
// Parameters:
//   current_line - An ATOM or HETATM line from a PDB file.
////////////////////////////////////////////////////////////////////////////////
void readATOM_OrHETATM_Line( string &current_line, MolFramework &structure, int counter, int &residueCounter){

  // Error check for very long lines. These are probably trashed files 
  // from DOS to UNIX transfers, or possibly malicious files aimed at 
  // the webserver. 
  //////////////////////////////////////////////////////////////////////
  if( current_line.size() > 200 ){
    cout << " Error: Line from PDB file " << structure.infile_name << " is longer than 200 characters." << endl;
    exit(1);
  }
  
  // Don't store psuedoatoms from *_FIRSTdataset files. 
  //////////////////////////////////////////////////////////////////////
  if( structure.using_dataset && 
      current_line.find("BMH") != string::npos )
    return;
  
  string temp;
  const char* errmsg;
  int int_atom_number;
  
  //cout << structure.infile_name << " " << counter << " " << current_line << endl;
  
  // Record name .
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].record_name = current_line.substr( 0, 6 );
  
  // Unique integer for use internally by FIRST.
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].FIRST_number = counter;
  
  // Atom serial number in orginal file.
  //////////////////////////////////////////////////////////////////////
  temp.assign( current_line, 6, 5 );
  errmsg = hy36decode(5, temp.c_str(), 5, &int_atom_number);
  if ( errmsg )
    {
      int_atom_number= atoi( temp.c_str() );
      if (parameters.verbose)
	{
	  cout << "WARNING: Atom [" << temp <<"] has non-numeric identifier" << endl;
	}       
    }
  
    structure.site_info[counter].orig_atom_number =  temp;

  
  // Map the original atom number to the number used by FIRST. They will
  // only be different if the atom numbers in the input PDB file did not
  // start with "1". 
  //////////////////////////////////////////////////////////////////////

  structure.orig_2_FIRST_atom_number.insert( pair<int,int> (int_atom_number, 
							    structure.site_info[counter].FIRST_number) );
  
  // Atom name.
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].atom_name = current_line.substr( 12, 4 );

  // Alternate location indicator
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].alt_location = current_line[16];
  
  // Residue name. 
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].residue_name = current_line.substr( 17, 3 );
  
  // Chain identifier in original file;
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].chain_ID = current_line[21];
  
  // Chain identifier used internally by FIRST.
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].FIRST_chain_ID = chain_counter;
  
  // Map the original chain ID to a unique integer chain ID used by FIRST
  //////////////////////////////////////////////////////////////////////
  structure.orig_2_FIRST_chain_ID.insert( pair<int,int> (structure.site_info[counter].chain_ID, 
							 structure.site_info[counter].FIRST_chain_ID) );
  
  // Residue sequence number.
  //////////////////////////////////////////////////////////////////////
  int int_seq_number;
  temp.assign( current_line, 22, 4 );
  errmsg = hy36decode(4, temp.c_str(), 4, &int_seq_number);
  if (errmsg){
    cout << "ERROR: Can not read sequence number [" <<  temp <<"]" <<endl;
      exit(1);
  }
  structure.site_info[counter].seq_number = int_seq_number;
  // if( structure.site_info[counter].record_name == "ATOM  " ) {
  structure.seq_numbers.insert( structure.site_info[counter].seq_number );
  //cout<<"setting: "<< counter << " "<<atoi( temp.c_str() )<< "\t " << structure.site_info[counter].residue_name<<endl;
  //}

  // Code for insertion of residues.
  //////////////////////////////////////////////////////////////////////
  structure.site_info[counter].insert_code = current_line[26];

  // Store a unique identifier for each residue, giving us a total count for
  // the number of residues. 
  ///////////////////////////////////////////////////////////////////////
  stringstream resnum_insertion_code_chain_ID;   
  resnum_insertion_code_chain_ID << structure.site_info[counter].seq_number << ";"
				 << structure.site_info[counter].insert_code << ";"
				 << structure.site_info[counter].chain_ID;

  // AJR 03.30.06 change second term to be the raw residue count;
  if( structure.unique_res_id.find( resnum_insertion_code_chain_ID.str() ) == structure.unique_res_id.end() ) { 
    structure.unique_res_id.insert( pair<string,int> (resnum_insertion_code_chain_ID.str(), residueCounter++) );
    //cout<<residueCounter<<" "<<counter<<"  "<<resnum_insertion_code_chain_ID.str()<<endl;
  }

  // X-coordinate in Angstroms.
  //////////////////////////////////////////////////////////////////////
  temp.assign( current_line, 30, 8 );
  structure.site_info[counter].coords[X] = atof( temp.c_str() );
  //  Compare to min/max in X direction. 
  if( structure.site_info[counter].coords[X] < structure.min_coords[X] )
    structure.min_coords[X] = structure.site_info[counter].coords[X];
  if( structure.site_info[counter].coords[X] > structure.max_coords[X] )
    structure.max_coords[X] = structure.site_info[counter].coords[X];
  
  // Y-coordinate in Angstroms. 
  //////////////////////////////////////////////////////////////////////
  temp.assign( current_line, 38, 8 );
  structure.site_info[counter].coords[Y] = atof( temp.c_str() );
  //  Compare to min/max in Y direction. 
  if( structure.site_info[counter].coords[Y] < structure.min_coords[Y] )
    structure.min_coords[Y] = structure.site_info[counter].coords[Y];
  if( structure.site_info[counter].coords[Y] > structure.max_coords[Y] )
    structure.max_coords[Y] = structure.site_info[counter].coords[Y];
  
  // Z-coordinate in Angstroms.  
  //////////////////////////////////////////////////////////////////////
  temp.assign( current_line, 46, 8 );
  structure.site_info[counter].coords[Z] = atof( temp.c_str() );
  //  Compare to min/max in Z direction. 
  if( structure.site_info[counter].coords[Z] < structure.min_coords[Z] )
    structure.min_coords[Z] = structure.site_info[counter].coords[Z];
  if( structure.site_info[counter].coords[Z] > structure.max_coords[Z] )
    structure.max_coords[Z] = structure.site_info[counter].coords[Z];
  
  // Occupancy
  //////////////////////////////////////////////////////////////////////
  if( current_line.size() >= 59 ){
    temp.assign( current_line, 54, 6 );
    structure.site_info[counter].occupancy = atof( temp.c_str() );
  }
  
  // Temperature factor.
  //////////////////////////////////////////////////////////////////////
  if( current_line.size() >= 65 ){
    temp.assign( current_line, 60, 6 );
    structure.site_info[counter].temp_factor = atof( temp.c_str() );
  }
  
  // Seg ID, columns 73-76
  //////////////////////////////////////////////////////////////////////
  if( current_line.size() >= 76 ) {
    temp.assign( current_line, 72, 4 );
    //seg_id field is left-justified according to pdb format,
    //so here we remove spaces from right-hand side
    temp.erase( temp.find_last_not_of(" ") + 1); // remove spaces
    structure.site_info[counter].seg_id = temp;
  }
      
  // Element name, columns 77-78
  //////////////////////////////////////////////////////////////////////
  if( current_line.size() >= 78 ){
    temp.assign( current_line, 76, 2 );
    //check if the element characters are valid (alphabetic or spaces are valid,
    //but two spaces are not valid)
    if( ( isalpha(temp[0]) || temp[0] == ' ' )
            && ( isalpha(temp[1]) || temp[1] == ' ' )
            && !(temp[0] == ' ' && temp[1] == ' ') ) { 
            // if the element characters were right-justified, make them left-justified
            if( temp[0] == ' ' ) swap( temp[0], temp[1] );
            // now copy the element name into storage
            structure.site_info[counter].element_name = temp;
            }
  }

  // Try to determine the element type if none was listed. 
  //////////////////////////////////////////////////////////////////////
  if( structure.site_info[counter].element_name.size() == 0 ||
      structure.site_info[counter].element_name[0] == ' ' )
    structure.site_info[counter].element_name = structure.getElementNamePDBFormat( structure.site_info[counter].atom_name );
  
  // Charge, columns 79-80
  if( current_line.size() >= 80 ){
    temp.assign( current_line, 78, 2 );
    structure.site_info[counter].charge = atof( temp.c_str() );
  }

  // FIRST_group_ID, columns 82-87
  // these columns are not defined by the PDB, but if there
  // is data in these columns and the -use_group_id option
  // is on, we read into the FIRST_group_ID field (integers only)
  //////////////////////////////////////////////////////////////////////
  if( parameters.use_group_id && current_line.size() >=87 ){
    temp.assign( current_line, 81, 6); 
    structure.site_info[counter].FIRST_group_ID = atoi( temp.c_str() );
  }

  // Remove padding from the atom name string.
  //////////////////////////////////////////////////////////////////////
  removeWhiteSpacePadding( structure.site_info[counter].atom_name );

  //cout << setw(8) << "[" << structure.site_info[counter].atom_name << "]" 
  //     << setw(8) << "[" << structure.site_info[counter].element_name << "]" << endl;

  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Read and store the CONECT records from a PDB formatted file.
//   The CovalentBonder class will use these to help determine bonding.
//   NOTE: The original atom numbers from the CONECT records are stored,
//   not the corresponding FIRST atom numbers.
////////////////////////////////////////////////////////////////////////////////
void readCONECT_Record(const string &current_line, MolFramework &structure ){

  //cout << current_line << endl;

  SiteID base_atom = atoi( (current_line.substr(6, 5)).c_str() );  
  SiteID neighbor = 0;

  int start = 11;

  //cout << "Storing CONECT " << base_atom << endl;
  // Read and store the neighbors of the base atom until we reach a newline. 
  ////////////////////////////////////////////////////////////////////////////////
  while( start < 31 &&
	 current_line.find_first_of("123456789", start) != string::npos ){

    neighbor = atoi( (current_line.substr(start, 5)).c_str() );
    //cout << "   [start:"<< start << ", end:" << end << "] [" << neighbor << "]" <<endl;

    if( neighbor != 0 ){
      vector<SiteID>::iterator knownConectPartner = find( structure.conect_records[base_atom].begin(), 
							  structure.conect_records[base_atom].end(), 
							  neighbor );
      if( knownConectPartner == structure.conect_records[base_atom].end() ){
	structure.conect_records[base_atom].push_back( neighbor );
      }
    }
    start += 5;
  }
  //cout << endl;
}
////////////////////////////////////////////////////////////////////////////////

