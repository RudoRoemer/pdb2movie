#include "global_defs.h"
#include "Parameters.h"
#include "generalUtils.h"
#include "interact.h"
#include "readPDB.h"
#include "readDataset.h"
#include "readMOL2.h"
#include "readAXYZ.h"
#include "readHARLEM.h"
#include "readBBG.h"
#include "MolFramework.h"
#include "PebbleGame.h"
#include "flexweb.h"
#include "output.h"
#include "Color.h"
#include "Cutoff.h"
#include "RigidCluster.h"
#include "RigidClusterFactory.h"
#include "RigidClusterAnalysis.h"
#include "hybrid_36_c.h"

#include <cstring>

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Print a logo screen displaying the name of the program, and version number.
////////////////////////////////////////////////////////////////////////////////
void printLogo(){

  cout << endl;
  cout << " ---------------------------------------------------------- " << endl;
  cout << "                           FIRST                            " << endl;
  cout << "     Floppy Inclusion and Rigid Substructure Topography     " << endl;
  cout << "                        VERSION " << FIRST_VERSION_ID << endl << endl;
  cout << "                           FRODA                            " << endl;
  cout << "      Framework Rigidity Optimised Dynamic Algorithm        " << endl;
  cout << "                        VERSION " << FRODA_VERSION_ID << endl << endl;
  cout << "                            and                     " << endl << endl;
  cout << "                           TIMME                            " << endl;
  cout << "  Tool for Identifying Mobility in Macromolecular Ensembles " << endl;
  cout << "                        VERSION " << TIMME_VERSION_ID << endl << endl;
  cout << "                Arizona State University            " << endl << endl;
  cout << " Copyright (c) 2004-2006 Brandon Hespenheide, Stephen Wells " << endl;       
  cout << "     Scott Menor and Dan Farrell. All Rights Reserved.      " << endl;       
  cout << " ---------------------------------------------------------- " << endl << endl;
  cout << " Random number generator: mt19937ar.cpp                     " << endl;
  cout << " Copyright (c) 1997 - 2002, Makoto Matsumoto and            " << endl;
  cout << "                            Takuji Nishimura, " << endl;
  cout << " All rights reserved.                                       " << endl;
  cout << " Copyright (c) 2005, Mutsuo Saito,                          " << endl;
  cout << " All rights reserved.                                       " << endl;
  cout << " ---------------------------------------------------------- " << endl << endl;
  cout << "Hybrid-36 PDB serial and sequence numbers: hybrid_36_c.cpp  " << endl;
  cout << "cctbx Copyright (c) 2006, The Regents of the University of  " << endl;
  cout << "California, through Lawrence Berkeley National Laboratory   " << endl;
  cout << "(subject to receipt of any required approvals from the U.S. " << endl;
  cout << "Dept. of Energy).  All rights reserved.                                            " << endl;
  cout << " ---------------------------------------------------------- " << endl << endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Determine what type of file has been read in. Call the appropriate routine
//   to read the data based on the input file type. 
////////////////////////////////////////////////////////////////////////////////
void readData( MolFramework &structure, string infile_name ){

  string file_extension;
  structure.infile_name = infile_name;
      
  size_t begin_extension = structure.infile_name.find_last_of(".");
  if( begin_extension != string::npos )
    file_extension = structure.infile_name.substr(begin_extension);
  else 
    file_extension = "";

  toLower( file_extension );

  // Find the path to file.
  //////////////////////////////////////////////////////////////////////
  size_t last_slash = structure.infile_name.find_last_of("/\\");
  structure.path_name.assign( structure.infile_name, 0, last_slash+1 );
  structure.base_name = structure.infile_name.substr( last_slash+1 );
  parameters.path_to_working_dir = structure.path_name;

  //cout << "infile name " << structure.infile_name << endl;
  //cout << "path name " << structure.path_name << endl;
  //cout << "base_name " << structure.base_name << endl;

  // If the user defined a NMR model to use in the command line, set the 
  // structure to be an NMR file.
  //////////////////////////////////////////////////////////////////////
  if( parameters.use_model_number != 0 ){
    structure.is_nmr = true;
    structure.model_number = parameters.use_model_number;
  }
  
  // PDB input file
  //////////////////////////////////////////////////////////////////////
  if( file_extension == ".pdb" || 
      file_extension == ".ent" ) {
    int name_length = (structure.base_name).size();
    structure.base_name.erase( name_length-4, 4 ) ;
    structure.ctf_file_name = structure.base_name + ".ctf";
    structure.using_pdb = true;
    readPDB_File( structure );
  }
  // old FIRSTdataset input file
  //////////////////////////////////////////////////////////////////////
  else if( structure.infile_name.find("_FIRSTdataset",0) != string::npos ){
    int name_length = (structure.base_name).size();
    structure.base_name.erase( name_length-13, 13 ) ;
    structure.ctf_file_name = structure.base_name + ".ctf";  
    structure.using_dataset = true;
    readDatasetFile( structure );
  }
  // Sybyl (c) Mol2 input file
  //////////////////////////////////////////////////////////////////////
  else if( structure.infile_name.find(".mol2",0) != string::npos ){
    cout << " Error: Mol2 file format compatibility not completed." << endl;
    exit(1);
    int name_length = (structure.base_name).size();
    structure.base_name.erase( name_length-4, 4 ) ;
    structure.ctf_file_name = structure.base_name + ".ctf";  
    structure.using_mol2 = true;
    readMOL2_File( structure );
  }
  // Generic ATOM - X - Y - Z input file
  //////////////////////////////////////////////////////////////////////
  else if( file_extension == ".axyz" ){
    int name_length = structure.base_name.size();
    structure.base_name.erase( name_length-5, 5 ) ;
    structure.using_axyz = true;    
    readAXYZ_File( structure );
  }
  // Generic body-bar graph file input
  //////////////////////////////////////////////////////////////////////
  else if( file_extension == ".bbg" ){
    int name_length = structure.base_name.size();
    structure.base_name.erase( name_length-4, 4 ) ;
    structure.using_bbg = true;    
    readBBG_File( structure );
  }
  // HARLEM input file
  //////////////////////////////////////////////////////////////////////
  else if( file_extension == ".hlm" ){
    int name_length = structure.base_name.size();
    structure.base_name.erase( name_length-4, 4 ) ;
    structure.using_harlem = true;    
    read_harlem_file( structure );
  }
  else{
    cout << " Unrecognized input file format." << endl;
    structure.unrecognized_file_format = true;
  }
 
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Currently, each line of the input file that corresponds to an atom or site
//   is stored in the site info array. This function will go through and remove
//   specific sites based on user input and/or default behaviour. Examples of
//   when sites would be removed are the existence of multiple side-chain 
//   conformations in the input file, or attempting to mimic the loss of a 
//   specific group or residue.
//
//   This function contains two loops over all the sites in the input structure.
//   The removed sites will not be present in the output. 
////////////////////////////////////////////////////////////////////////////////
void excludeSites( MolFramework &structure, bool shrink_and_remap_only ){

  if( !shrink_and_remap_only ){

    // Look for alternative side chain lables
    //////////////////////////////////////////////////////////////////////
    excludeAltSideChainLocation( structure );
  }
  
  // Loop through the site_info array. If we find a site that is labeled
  // excluded, increment the "gap" counter, and swap the excluded site 
  // with the next in the array. Continue swapping element N with element 
  // N+gap, incrementing gap every time we exclude a site. This should
  // push all the excluded sites to the end of the array. Remember to 
  // check the last site in the array explicitly after the loop, as the
  // logic doesn't allow for swapping the next element past the end of 
  // the array. 
  // 
  // After shrinking the site_info array, reset the number of total_sites,
  // and remap the internally used FIRST number with the original atom
  // number, should one exist. 
  // BMH need to put in explicit copy constructor in Site_Info.
  //////////////////////////////////////////////////////////////////////
  int gap = 0;
  // int excluded = 0; // TODO - removeme?
  Site_Info temp_info;
  
  for(unsigned int siteNumber = 1; siteNumber < (structure.total_sites -gap); siteNumber++ ){

    gap = 0;

    if( structure.site_info[siteNumber].excluded ){
      while( structure.site_info[siteNumber +gap].excluded &&
	     siteNumber+gap < structure.total_sites ) 
	gap++;

      temp_info = structure.site_info[siteNumber];
      structure.site_info[siteNumber] = structure.site_info[siteNumber +gap];
      structure.site_info[siteNumber +gap] = temp_info;
    }
  }

  //for( int a = 1; a <= structure.total_sites; a++ ){
  //cout << setw(5);
  //structure.site_info[a].print();
  //}

  if( structure.site_info[structure.total_sites].excluded )
    structure.total_sites -= --gap;
  else{
    temp_info = structure.site_info[structure.total_sites];
    structure.site_info[structure.total_sites] = structure.site_info[structure.total_sites -gap];
    structure.site_info[structure.total_sites -gap] = temp_info;
    structure.total_sites -= gap;
  }

  remapOrigToFIRSTNumbers( structure );

  //for( int a = 1; a <= structure.total_sites+gap; a++ ){
  //cout << setw(5);
  //structure.site_info[a].print();
  //}
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: Retains only one unique location for each side chain and main 
//              chain atom in the pdbfile by excluding alternate atom positions 
//              from the structure.
////////////////////////////////////////////////////////////////////////////////
void excludeAltSideChainLocation( MolFramework &structure ){

  // note should remove alt mainchain locations as well -- default is to keep 
  // 1st location (A) for each atom with alt_location != ' ';
  for( unsigned int siteNumber = 1; siteNumber < structure.total_sites; siteNumber++ ){
    if( structure.site_info[siteNumber].alt_location != " ") {
      if(structure.site_info[siteNumber].alt_location > "siteNumber") structure.site_info[siteNumber].excluded = true;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This routine is called from within the "exclude_sites()" function. If sites
//   have been excluded from the network, all the FIRST numbers (which are used 
//   internally to keepm track of each site) are reordered from 1 to N. This 
//   function will remap the original atom numbers (found in the input file) to 
//   the newly assigned, internally-used, FIRST numbers.
////////////////////////////////////////////////////////////////////////////////
void remapOrigToFIRSTNumbers( MolFramework &structure ){

  structure.orig_2_FIRST_atom_number.clear();
  string temp1;
  const char* errmsg;
  int int_atom_number;
  int slen;

  for( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ){
 
    // Unique integer for use internally by FIRST.
    //////////////////////////////////////////////////////////////////////
    structure.site_info[siteNumber].FIRST_number = siteNumber;
    
    // Map the origianl atom number to the number used by FIRST. They will
    // only be different if the atom numbers in the input PDB file did not
    // start with "1". 
    //////////////////////////////////////////////////////////////////////
    temp1= structure.site_info[siteNumber].orig_atom_number;
    slen = strlen ( temp1.c_str());
    errmsg = hy36decode(5, temp1.c_str(), temp1.size(), &int_atom_number);
    if ( errmsg )
      {
	if (isNumber( temp1 )){
	  int_atom_number = atoi( temp1.c_str() );
	  structure.orig_2_FIRST_atom_number.insert( pair<int,int> (int_atom_number, siteNumber) );
    //cout << "new mapping " << structure.site_info[siteNumber].orig_atom_number << " - " << siteNumber << endl;
	}
	  else{
	    cout << "ERROR: Atom [" << temp1 <<"] has invalid number" << endl;
	    exit(1);
	  }
      }
    else {
      structure.orig_2_FIRST_atom_number.insert( pair<int,int> (int_atom_number, siteNumber) );
    }
  }

  structure.reset_site_info_labels();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void setVdwRadii( MolFramework &structure ){

  for( unsigned int siteNumber = 1; siteNumber < structure.total_sites; siteNumber++ ){

    if     ( structure.site_info[siteNumber].element_name == "H " ) structure.site_info[siteNumber].vdw_radius = 1.00;
    else if( structure.site_info[siteNumber].element_name == "C " ) structure.site_info[siteNumber].vdw_radius = 1.70;
    else if( structure.site_info[siteNumber].element_name == "N " ) structure.site_info[siteNumber].vdw_radius = 1.55;
    else if( structure.site_info[siteNumber].element_name == "O " ) structure.site_info[siteNumber].vdw_radius = 1.40;
    else if( structure.site_info[siteNumber].element_name == "S " ) structure.site_info[siteNumber].vdw_radius = 1.80;
    else if( structure.site_info[siteNumber].element_name == "P " ) structure.site_info[siteNumber].vdw_radius = 1.80;
    else if( structure.site_info[siteNumber].element_name == "MG" ) structure.site_info[siteNumber].vdw_radius = 1.50;
    else if( structure.site_info[siteNumber].element_name == "SI" ) structure.site_info[siteNumber].vdw_radius = 2.10;
    else if( structure.site_info[siteNumber].element_name == "MN" ) structure.site_info[siteNumber].vdw_radius = 1.40;
    else structure.site_info[siteNumber].vdw_radius = 1.50; // Default radii for unrecognized atom types
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//     
void buildFramework( MolFramework &structure ){
////////////////////////////////////////////////////////////////////////////////

  // These functions are always called for any input type OTHER than 
  // body-bar graphs. All non body-bar graph file types should represent
  // some type of chemical structure. 
  //////////////////////////////////////////////////////////////////////
  if( !structure.using_bbg ){
    structure.assignAtomsTo3DGrid();
    structure.identifyCovalentBonds();
    structure.countNeighbors();
    structure.assignHydrogenBondStatus();

    if (parameters.findExposed) {
      if( parameters.verbose ){
	cout << "Calling findExposed" << endl;
      }
      findExposedAtoms ( structure );
      if( parameters.verbose ){
	cout << "Finished findExposed" << endl;
      }
    }
  }

  // Additional processing of the input information will depend on the input
  // file type used. Each allowed input file type sets a unique flag, generally
  // in the form "using_*", where * describes the input file type. For each 
  // recognized type, perform the necesary operations to create the body-bar
  // network.
  //////////////////////////////////////////////////////////////////////
  if( (structure.using_pdb || structure.using_harlem) && !parameters.run_timme ){
    structure.identifyHydrogenBonds();
    structure.identify_stacked_rings();
    structure.identifyHydrophobicTethers();
    structure.addUserDefinedConstraints();
  }
  else if( structure.using_dataset ){
  }
  else if( structure.using_axyz ){
    structure.identifyHydrogenBonds();
    structure.identifyHydrophobicTethers();
    structure.addUserDefinedConstraints();
  }

  // If we're running in interactive mode, give the user some options to check
  // the file before heading to the pebble game. Also, perform some structure 
  // validation tests. 
  //////////////////////////////////////////////////////////////////////
  structure.queryNetwork();

  // Run some simple error checks on the network.
  //////////////////////////////////////////////////////////////////////
  if( !structure.using_bbg )
    structure.validateStructure();

  // Output bond lists. Have to output the covalent bonds list before the
  // prepare_network_for_analysis, as all the hbonds and hphobes will be 
  // added to the site_info array then. 
  //////////////////////////////////////////////////////////////////////
  if( parameters.phout )
    structure.outputHydrophobicTetherList();

  if( parameters.hbout )
    structure.outputHydrogenBondList();

  if( parameters.srout )
    structure.outputStackedRingsList();
  
  if( parameters.covout )
    structure.outputCovalentBondList();
  
  // Add noncovalent bonds to the neighbor list in site_info. Also, 
  // identify which bonds, if any, can be diluted from the network. Add these to a 
  // vector named "dilution_list". Also, allow the user to add their own "ranking" 
  // to each bond, which will be used if we are running a bond dilution.
  //////////////////////////////////////////////////////////////////////
  if (!parameters.run_timme) {
    structure.prepareFrameworkForAnalysis();
  }

  if( parameters.flexweb ) {
    output_details_flexweb();
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void analyzeFlexibilityAndOutput( MolFramework &structure ){

  if( parameters.bond_dilution ){
    switch( parameters.bond_dilution ){
    case 1: hydrogenBondDilutionEnergyOnly( structure );
      break;
    case 2: hydrogenBondDilutionAssemblyCodeEnergyOnly( structure );
      break;
    }

    // BMH hack to get FRODA to run after a bond dilution.
    parameters.bond_dilution = false;
  }

  PebbleGame pebble_game( &structure );

  pebble_game.runPebbleGame( RUN_DECOMPS );
  
  structure.pruneSideChains();
  structure.pruneDanglingEnds();
  
  structure.computeMeanCoordination( true );
  //structure.saveTopology( structure.dilution_list.size() );
  
  if( structure.using_axyz ){
    structure.outputRCD_MOL2FORMAT();
    structure.outputPyMolScriptRCD_MOL2Format();
  }
  else if( structure.using_harlem ){
    structure.outputHARLEMFormat();
    structure.outputRCD_PDBFormat( "_RCD.pdb" );
    structure.outputPyMolScriptRCD_PDBFormat();
  }
  else if( structure.using_pdb ||
	   structure.using_dataset ) {
    structure.outputRCD_PDBFormat( "_RCD.pdb" );
    structure.outputPyMolScriptRCD_PDBFormat();
  }
  
  structure.outputRawData();
  structure.outputBondData();
  structure.outputText();


  // Create output specific for FlexWeb
  //////////////////////////////////////////////////////////////////////
  if( parameters.flexweb ){
    if (parameters.runFRODA) {
      output_FRODA_TIMME_Jmol_script(structure);
    } 
    else {
      output_Jmol_script( structure );
    }
    output_RCD_XML( structure );
  }
  
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void hydrogenBondDilutionAssemblyCodeEnergyOnly( MolFramework &structure ){

  cout << "Not yet implemented" << endl;
  return;
  /*
  RigidClusterDilution hydrogenBondDilution;
  hydrogenBondDilution.setDisassemble( true );
   
  PebbleGame pebble_game( &structure );
  
  structure.pruneSideChains();
  structure.pruneDanglingEnds();
  pebble_game.runPebbleGame( NO_DECOMPS );
  structure.computeMeanCoordination( true );
  structure.reset_site_info_labels();
  structure.reset_pruned_label();
  pebble_game.resetPebbleGame();

  pebble_game.runPebbleGame( RUN_DECOMPS );
  structure.computeMeanCoordination();
  structure.reset_site_info_labels();
  pebble_game.resetPebbleGame();
  
  hydrogenBondDilution.addStep();
  hydrogenBondDilution.addPieces( structure.RC_atom_list );

  for( unsigned int current_bond = 0; current_bond < structure.dilution_list.size(); current_bond++ ){

    if( parameters.verbose )
      cout << "Breaking bond " << structure.dilution_list.size() -current_bond << endl;
    
    structure.removeBond( current_bond );
    
    structure.pruneSideChains();
    structure.pruneDanglingEnds();

    pebble_game.runPebbleGame( NO_DECOMPS );
    structure.computeMeanCoordination( true );


    structure.reset_site_info_labels();
    pebble_game.resetPebbleGame();
    structure.reset_pruned_label();
    
    pebble_game.runPebbleGame( RUN_DECOMPS );
    structure.computeMeanCoordination();
    
    hydrogenBondDilution.addStep();
    hydrogenBondDilution.addPieces( structure.RC_atom_list );
        
    structure.reset_site_info_labels();
    pebble_game.resetPebbleGame();
  }

  // Add all the bonds back to the network in case we're planning to  do
  // a FRODA run. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int current_bond = 0; current_bond < structure.dilution_list.size(); current_bond++ ){
    structure.addBond( current_bond );
  }  

  MainChain1D_AssemblyGraphics testGraphics( hydrogenBondDilution );
  testGraphics.loadGraphicsFileShapes("PS");
  testGraphics.mapStructureToAssemblyCoordinates( structure );
  testGraphics.printAssemblyToFile("testoutput");
  */
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void hydrogenBondDilutionEnergyOnly( MolFramework &structure ){

  // Remove any previous temporary files.
  //////////////////////////////////////////////////////////////////////
  system( "\\rm -f hbdil_temp.pdb decomp_list" );
  
  PebbleGame pebble_game( &structure );
  
  structure.pruneSideChains();
  structure.pruneDanglingEnds();
  pebble_game.runPebbleGame( NO_DECOMPS );
  structure.computeMeanCoordination( true );
  structure.reset_site_info_labels();
  structure.reset_pruned_label();
  pebble_game.resetPebbleGame();

  pebble_game.runPebbleGame( RUN_DECOMPS );
  structure.computeMeanCoordination();
  structure.saveTopology( structure.dilution_list.size() );
  structure.reset_site_info_labels();
  pebble_game.resetPebbleGame();
  
  multimap<Cutoff, set<SiteID> > mapFromCutoffToJoinedSites;
  unsigned int totalJoinedPairs = 0;
  
  //////////////////////////////////////////////////////////////////////
  for( unsigned int current_bond = 0; current_bond < structure.dilution_list.size(); current_bond++ ){
    
    if( parameters.verbose )
      cout << "Breaking bond " << structure.dilution_list.size() -current_bond << endl;
    
    structure.removeBond( current_bond );
    
    structure.pruneSideChains();
    structure.pruneDanglingEnds();

    pebble_game.runPebbleGame( NO_DECOMPS );
    structure.computeMeanCoordination( true );

    structure.reset_site_info_labels();
    pebble_game.resetPebbleGame();
    structure.reset_pruned_label();
    
    pebble_game.runPebbleGame( RUN_DECOMPS );
    structure.computeMeanCoordination();
    //compareTopologies();
    structure.saveTopology( current_bond );
    
    if( parameters.diluteWithRigidClusterFactory ) {

      // BMH why do we need a new copy? 
      //////////////////////////////////////////////////////////////////////
      vector< vector<int> > RC_atom_list = structure.RC_atom_list;
      set< set<SiteID> > rigidClusterDecomposition;
      
      // add all rigid clusters containing more than one site to the current 
      // rigidClusterDecomposition
      //////////////////////////////////////////////////////////////////////
      for( vector< vector<int> >::iterator RC_atom_list_iterator = RC_atom_list.begin();
	   RC_atom_list_iterator != RC_atom_list.end();
	   RC_atom_list_iterator++ ){

        set<SiteID> cluster( RC_atom_list_iterator->begin(), 
			     RC_atom_list_iterator->end() );
        
        if (cluster.size() > 1) {
          rigidClusterDecomposition.insert(cluster);
	}
      }
      
      Cutoff cutoff = structure.dilution_list[current_bond].energy;
      //Cutoff averageR = structure.mean_coordination;
      
      for( set< set<SiteID> >::iterator clusterIterator = rigidClusterDecomposition.begin();
           clusterIterator != rigidClusterDecomposition.end();
           clusterIterator++ ){          

	// BMH do we need to create another variable here? 
	set<SiteID> joinedSites = *clusterIterator;
	totalJoinedPairs += joinedSites.size();
	mapFromCutoffToJoinedSites.insert( pair<Cutoff, set<SiteID> > (cutoff, joinedSites));
      }
    }
    
    structure.reset_site_info_labels();
    pebble_game.resetPebbleGame();
  }

  // Add all the bonds back to the network in case we're planning to  do
  // a FRODA run. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int current_bond = 0; current_bond < structure.dilution_list.size(); current_bond++ ){
    structure.addBond( current_bond );
  }  

  if (!parameters.diluteWithRigidClusterFactory) {

    
    structure.MakeTmpFileForHBDilute();
    string command = parameters.path;
    command += "/hbdilute/hbdilute " + structure.path_name + "decomp_list b " + structure.path_name + "hbdil_temp.pdb";
    system( command.c_str() );
    //command = "\\rm " + structure.path_name + "hbdil_temp.pdb " + structure.path_name + "decomp_list";
    //system( command.c_str() );
    
    if( structure.path_name.size() ){
      command = "\\mv -f " + structure.base_name + ".ps " + structure.path_name + structure.base_name + ".ps";
      system( command.c_str() );
    }
  }
  
  if (parameters.diluteWithRigidClusterFactory) {

    RigidClusterFactory rigidClusterFactory;
    
    rigidClusterFactory.decomposeIntoRigidClustersOLD(mapFromCutoffToJoinedSites);
    
    RigidClusterHierarchy *rigidClusterHierarchy = rigidClusterFactory.getRigidClusterHierarchy();
    
    RigidClusterAnalysis rigidClusterAnalysis(rigidClusterHierarchy);
    rigidClusterAnalysis.setStructure(&structure);
    rigidClusterAnalysis.setXAxisType(parameters.dilutionPlotXAxisType);
    
    rigidClusterAnalysis.assignRigidClusterLabels(); // FIXME - is there a way to do this while joining sites?
    
    rigidClusterAnalysis.computeRigidClusterColoring();
    
    string dilutionFilename = structure.path_name + structure.base_name + "_dilution";
    rigidClusterAnalysis.saveDilution(dilutionFilename); //output the Stripe Plot line
    
    string dilutionPdbFilename = structure.path_name + structure.base_name + "_dilution.pdb";
    ofstream dilutionPdbFile(dilutionPdbFilename.c_str());

    map<SiteID, float> emptyMap;	

    rigidClusterAnalysis.savePDBmodel(&structure,
                                      emptyMap,
                                      emptyMap,
                                      dilutionPdbFile);
    dilutionPdbFile.close();
    
    cout << "FIXME - implement new rigid cluster decomposition" << endl;
    structure.reset_site_info_labels();
    pebble_game.resetPebbleGame();
  }  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Convert a string to all upper-case letters.
// Parameters:
//   str: Input string.
////////////////////////////////////////////////////////////////////////////////
string toUpper( string &str ){

  for( unsigned int a = 0; a < str.length(); a++ ){
    if( str[a] >= 'a' && str[a] <= 'z' )
      str[a] -= 32;
  }
  
  return(str);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Convert a string to all lower-case letters.
// Parameters:
//   str: Input string.
////////////////////////////////////////////////////////////////////////////////
string toLower( string &str ){

  for( unsigned int a = 0; a < str.length(); a++ ){
    if( str[a] >= 'A' && str[a] <= 'Z' )
      str[a] += 32;
  }
  
  return(str);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: Generic string tokenizer. The function takes as arguments the
//   string to be parsed and an optional string containing the delimiters. By 
//   default, the function will assume white-space characters as delimiters. The
//   function returns a vector of strings containing the tokens parsed from the
//   input string. NOTE: If your input string contains newline characters, these
//   will be treated as a delimiter. The string will continue past intervening 
//   newlines until it reaches the end of the string, as determined by string::npos.
////////////////////////////////////////////////////////////////////////////////
vector<string> tokenize_string( string &input, string delimiters ){

  // Add the newline character to the list of tokens, just in case it 
  // wasn't included.
  //////////////////////////////////////////////////////////////////////
  if( delimiters.find("\n") == string::npos )
    delimiters += "\n";

  vector<string> tokens;

  size_t start = input.find_first_not_of(delimiters,0);
  size_t end   = input.find_first_of(delimiters,start);

  if( start == string::npos )
    return( tokens );

  while( end != string::npos ){
    tokens.push_back( input.substr(start, end-start) );
    start = input.find_first_not_of( delimiters, end );
    end   = input.find_first_of( delimiters, start );
  }

  if( start != string::npos )
    tokens.push_back( input.substr(start, end-start) );

  // Testing
  //////////////////////////////////////////////////////////////////////
  //cout << "initial string ["<< input << "]" << endl;
  //cout << tokens.size() << " tokens:" << endl;
  //for( int a = 0; a < tokens.size(); a++ )
  //cout << "\t" << " " << tokens[a] << endl;

  return( tokens );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: Remove leading and trailing white-space characters from a 
//   string. 
////////////////////////////////////////////////////////////////////////////////
void removeWhiteSpacePadding( string &str ){

  size_t start = str.find_first_not_of(" \t");
  size_t end   = str.find_last_not_of(" \t\n");
  if ( start != str.npos && end != str.npos ) {
  //cout << str << " [" << str.substr(start, end-start+1) << "] " << start << " " << end << endl;
    str = str.substr( start, end-start+1 );
  }
  else str = "";
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool sortOnSize( const vector<unsigned int> &a, const vector<unsigned int> &b ){
  return( a.size() >  b.size() );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string hr(){
  return("# --------------------------------------------------------------------");
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string next_rcd_color_as_name( bool reset ){

  static int new_color = -1;

  if( reset ){
    new_color = -1;
    return("");
  }

  new_color++;
    
  return( Color::getPymolColor(new_color) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string next_jmol_color_as_name( bool reset ){

  static int new_color = -1;
	
  if( reset )
    new_color = -1;

  new_color++;
    
  return( Color::getJmolColor(new_color) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
string next_svg_color_as_name( bool reset ){
	
	static int new_color = -1;
	
	if( reset )
		new_color = -1;
	
	new_color++;
    
	return( Color::getSVGColor(new_color));
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//
////////////////////////////////////////////////////////////////////////////////
bool isComment( string &linebuf ){

  removeWhiteSpacePadding( linebuf );
  if( linebuf.at(0) == '#' )
    return( true );

  return( false );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This functions assumes that the input structures are from the Yale Morph 
//   Server. These files have been "aligned" by renaming and renumbering the 
//   residue name and sequence number. This function will loop over the first
//   structure, and search for a similar atom_name/res_name/res_number in the
//   second structure. If a pair is found, both atoms are renumber incrementally
//   from -1 backwards, to -N. Then each positive atom number is relabeled from +(N+1).
//   Finally, the sign of the negative atom numbers is flipped. The result is two
//   structures that have atoms 1-N targeted (the -fancy command-line flag must be
//   on). 
////////////////////////////////////////////////////////////////////////////////
void renumber_structures_from_morph_server( MolFramework &structure, MolFramework &target ){

  int counter = -1;
  int aray_a[structure.total_sites];
  int aray_b[structure.total_sites];

  if( parameters.verbose )
    cout << " Renumbering atom numbers for FRODA targeted dynamics..." << endl;

  if( !parameters.using_target ){
    cout << " Error: The morphing function requires a target structure." << endl;
    cout << "        Please use the -target<file_name> command-line option." << endl;
    exit(1);
  }

  if( !parameters.fancy_target ){
    cout << " Error: The morphing function requires fancy targeting." << endl;
    cout << "        Please use the -fancy command-line option." << endl;
    exit(1);
  }

  // Find atoms to target
  //////////////////////////////////////////////////////////////////////
  for( unsigned int a = 1; a < structure.total_sites; a++ ){
    for( unsigned int b = 1; b < target.total_sites; b++ ){

      if( !target.site_info[b].checked &&
	  structure.site_info[a].FIRST_group_ID == target.site_info[b].FIRST_group_ID &&
	  structure.site_info[a].chain_ID == target.site_info[b].chain_ID &&
	  structure.site_info[a].seq_number == target.site_info[b].seq_number &&
	  structure.site_info[a].residue_name == target.site_info[b].residue_name &&
	  structure.site_info[a].atom_name == target.site_info[b].atom_name ){
	//cout << "found pair: " <<  structure.site_info[a].orig_atom_number 
	//  << " " << target.site_info[b].orig_atom_number << endl;
	//cout << "found pair: " <<  structure.site_info[a].FIRST_number 
	//     << " " << target.site_info[b].FIRST_number << endl;
 	/*cout << structure.site_info[a].FIRST_group_ID << " " << structure.site_info[a].chain_ID 
	     << " " << structure.site_info[a].seq_number << " " <<  structure.site_info[a].residue_name
	     << " " << structure.site_info[a].atom_name << " -- " << target.site_info[b].FIRST_group_ID
	     << " " << target.site_info[b].chain_ID << " " << target.site_info[b].seq_number
	     << " " << target.site_info[b].residue_name << " " << target.site_info[b].atom_name << endl;*/
	
	//	structure.site_info[a].orig_atom_number =  counter;
	aray_a[a] =  counter;
	//	target.site_info[b].orig_atom_number =  counter--;
	aray_b[b] =  counter--;
	target.site_info[b].checked = true;
	break;	
      }
    }
  }

  target.resetCheckedLabel();
  std::stringstream ss;

  // Relabel targeted and nontargeted atoms.
  //////////////////////////////////////////////////////////////////////
  counter = abs( counter );
  //char* temp_str;

  for( unsigned int a = 1; a < structure.total_sites; a++ ){
    if( aray_a[a] > 0)
      aray_a[a] = counter++;
    //    if( structure.site_info[a].orig_atom_number > 0 )
    //  structure.site_info[a].orig_atom_number =   counter++;
    else
      // structure.site_info[a].orig_atom_number =  abs(structure.site_info[a].orig_atom_number);
      aray_a[a] = abs( aray_a[a] );

    // itoa( aray_a[a], structure.site_info[a].orig_atom_number, 10 );

    ss<< aray_a[a];
    structure.site_info[a].orig_atom_number =  ss.str();
    //structure.site_info[a].orig_atom_number = (string) temp_str;
  }

  for( unsigned int b = 1; b < target.total_sites; b++ ){
    if( aray_b[b] > 0)
      aray_b[b] = counter++;

    // if( target.site_info[b].orig_atom_number > 0 )
    //  target.site_info[b].orig_atom_number =  counter++;
    else
      //  target.site_info[b].orig_atom_number =  abs(target.site_info[b].orig_atom_number);
      aray_b[b] = abs( aray_b[b] );
    // itoa( aray_b[b], structure.site_info[b].orig_atom_number, 10 );
    ss<< aray_b[b];
    structure.site_info[b].orig_atom_number =  ss.str();


  }

}
////////////////////////////////////////////////////////////////////////////////
 
////////////////////////////////////////////////////////////////////////////////
// Description:
//   Searches the space using a dot grid and labels those atoms that are exposed
//   to solvent as isExposed
////////////////////////////////////////////////////////////////////////////////
void findExposedAtoms( MolFramework &structure ){

  //cerr << "Entering findExposed" << endl;

  vector<Dot> dot; //points used for sampling
  double rangeX; //spread of coordinates
  double rangeY;
  double rangeZ;
  double pad = 5.0; //margin around actual atom positions, Angstroms
  double spacing = 0.4; //dot separation
  int nX, nY, nZ; // maximum index of dots
  int nTotal; //size of 1-D dot array
  double solventRadius = 1.4; //radius of water molecule

  int countRange = (int) (4.0/spacing); //check indices within 4 angstroms

  //work out the limits for the dot grid
  rangeX = structure.max_coords[X] - structure.min_coords[X] + 2*pad;
  rangeY = structure.max_coords[Y] - structure.min_coords[Y] + 2*pad;
  rangeZ = structure.max_coords[Z] - structure.min_coords[Z] + 2*pad;

  nX  = (int) ( rangeX/spacing);
  nY  = (int) ( rangeY/spacing);
  nZ  = (int) ( rangeZ/spacing);
  nTotal = nX*nY*nZ;

  //cerr << "Making dot array of size " << nTotal << endl;
  dot.resize(nTotal);
  //cerr << "Made dot array" << endl;

  //plot some positions for the dot grid
  for (int countX = 0; countX < nX; countX++){
    for (int countY = 0; countY < nY; countY++){
      for (int countZ = 0; countZ < nZ; countZ++){
        int index =(countZ*nX*nY) + (countY*nX) + countX;
        dot.at(index).position.x = ( structure.min_coords[X] - pad ) + (countX*spacing);
        dot.at(index).position.y = ( structure.min_coords[Y] - pad ) + (countY*spacing);
        dot.at(index).position.z = ( structure.min_coords[Z] - pad ) + (countZ*spacing);
        dot.at(index).nearestAtom = 0;
        dot.at(index).distance2 = 0;
        dot.at(index).going = false;
        dot.at(index).gone = false;
      }
    }
  }
  //cerr << "placed all dots" << endl;

  for ( unsigned int atomA = 1; atomA <= structure.total_sites; atomA ++){
    //cerr << "Working on atom " << atomA << endl;
    double cut2 = structure.site_info[atomA].vdw_radius + solventRadius;
    //cerr << "Relevant distance " << cut2 << endl;
    cut2 *= cut2; //square of cutoff distance

    Vector atomPos;
    atomPos.x = structure.site_info[atomA].coords[0];
    atomPos.y = structure.site_info[atomA].coords[1];
    atomPos.z = structure.site_info[atomA].coords[2];
    //cerr << "Made it to after atomPos" << endl;

    int atomX = (int) ( ( atomPos.x - structure.min_coords[X] + pad ) /spacing );
    int atomY = (int )( ( atomPos.y - structure.min_coords[Y] + pad ) /spacing);
    int atomZ = (int )( ( atomPos.z - structure.min_coords[Z] + pad ) /spacing);
    //int atomIndex = (atomZ*nX*nY) + (atomY*nX) + atomX;

    //cerr << "Made it to start of count loop" << endl;
    for ( int theX = atomX - countRange; theX < atomX + countRange; theX++){
      for ( int theY = atomY - countRange; theY < atomY + countRange; theY++){
        for ( int theZ = atomZ - countRange; theZ < atomZ + countRange; theZ++){
          int theIndex =(theZ*nX*nY) + (theY*nX) + theX;
          Vector delta = dot.at(theIndex).position;
          delta.x -= atomPos.x;
          delta.y -= atomPos.y;
          delta.z -= atomPos.z;
          double distance2 = dotProduct(delta, delta);
          if ( dot.at(theIndex).nearestAtom == 0 ){
            if ( distance2 < cut2 ){
              dot.at(theIndex).nearestAtom = atomA;
              dot.at(theIndex).distance2 = distance2;
              //cerr << "Dot " << theIndex << " owns atom " << atomA << endl;
            }
          }
          else if ( distance2 < dot.at(theIndex).distance2 ){
            dot.at(theIndex).nearestAtom = atomA;
            dot.at(theIndex).distance2 = distance2;
            //cerr << "Dot " << theIndex << " owns atom " << atomA << endl;
          }

        }
      }
    }


  }

  //now every point either belongs to an atom or not
  //remove all the points with no atom (solvent points)

  //cerr << "Stripping out solvent points." << endl;
  for ( int index = 0; index < nTotal; index++){
    if ( dot.at(index).nearestAtom == 0){
      dot.at(index).gone = true;
    }
  }

  //cerr << "Stripping out interior points." << endl;
  for ( int index = 0; index < nTotal; index++){
    if ( dot.at(index).gone ) continue;
    //this strips off all the outer solvent
    //check +/- x,y,z
    int indexToCheck;
    int nGone = 0;
    //first the x direction
    indexToCheck = index +1;
    if ( dot.at(indexToCheck).gone ) {
      nGone++;
    }
    indexToCheck = index -1;
    if ( dot.at(indexToCheck).gone ) {
      nGone++;
    }
    //now the +/- y direction
    indexToCheck = index +nX;
    if ( dot.at(indexToCheck).gone ) {
      nGone++;
    }
    indexToCheck = index -nX;
    if ( dot.at(indexToCheck).gone ) {
      nGone++;
    }
    //and now the +/- z direction
    indexToCheck = index +(nX*nY);
    if ( dot.at(indexToCheck).gone ) {
      nGone++;
    }
    indexToCheck = index -(nX*nY);
    if ( dot.at(indexToCheck).gone ) {
      nGone++;
    }
    if ( nGone == 0) {
      dot.at(index).going = true;
    }
  }

  //any dots that are not Gone or Going are surface dots
  //label them, and their atoms

  //cerr << "Identifying surface atoms." << endl;
  for ( int index = 0; index < nTotal; index++){
    if ( dot.at(index).gone || dot.at(index).going ) {
      continue;
    }
    unsigned int atomA = dot.at(index).nearestAtom;
    structure.site_info[atomA].isExposed = true;
    //cerr << "Atom " << atomA << " is exposed." << endl;

    for ( unsigned int nay = 0; nay < structure.site_info[atomA].neighbor_list.size(); nay++) {
      unsigned int atomB = structure.site_info[atomA].neighbor_list.at(nay);
      if ( structure.site_info[atomB].isExposed ) continue; 
      Vector delta = dot.at(index).position;
      delta.x -= structure.site_info[ atomB].coords[0];
      delta.y -= structure.site_info[ atomB].coords[1];
      delta.z -= structure.site_info[ atomB].coords[2];

      double cut2 = structure.site_info[atomB].vdw_radius + solventRadius;
      cut2 *= cut2; //square of cutoff distance
      if (  dotProduct( delta, delta ) < cut2 ) {
        structure.site_info[ atomB].isExposed = true;
        //cerr << "Atom " << atomB << " is exposed." << endl;
      }
    }

  }

  int nexposed=0;
  int ninterior=0;
  for ( unsigned int atomA = 1; atomA <= structure.total_sites; atomA ++){
    if (structure.site_info[ atomA ].isExposed ){
      nexposed++;
    }
    else ninterior++;
  }
  cerr << "Exposed: " << nexposed << ", interior " << ninterior << endl;

}
 

bool isNumber( string &str ){
  int i;
  bool res= true;
  for (i=0; i< (int) strlen(str.c_str() ); i++)
    {
      if( !isdigit(str[i]) )
	{
	  break;
	  return false;
	}
	 }
    return res;
  }

