#include "global_defs.h"
#include "Parameters.h"
#include "generalUtils.h"
#include "flexweb.h"
#include "MolFramework.h"
#include "PebbleGame.h"
#include "Froda.h"
#include "TIMME.h"
#include "RigidUnitSystem.h"
#include "Froda2Builder_FromMolFramework.h"
#include "SymReplica.h"
#include "SymmetryMatrices.h"

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Main.cpp
//
//   Brandon Hespenheide, Ph.D,
//   Center for Biological Physics
//   Arizona State University
//   Biophysics Theory Group
//   PO Box 871504
//   Tempe, AZ 85287
//   bmh@asu.edu
//
//   Copyright (c) 2004-2006 Brandon Michael Hespenheide, Stephen Wells,
//                           Scott Menor and Dan Farrell
////////////////////////////////////////////////////////////////////////////////

// Global variables.
////////////////////////////////////////////////////////////////////////////////
Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
int main( int argc, char **argv ) {

  if( argc == 1 ) 
    parameters.printUsage(); 
  
  // Process the command-line arguments. Creates necessary file names, 
  // and checks for discrepencies in the list of command-line arguments.
  // Parameters that are set by using command-line arguments override those found 
  // in the parameters file FIRST.par, which override the default values set
  // when the parameters data structure is initialized in ** global_defs.h **
  //////////////////////////////////////////////////////////////////////
  parameters.readCommandLineOptions( argc, argv );

  if( parameters.flexweb ){
    setTaskNames();
  }
  
  // Print the program logo.
  //////////////////////////////////////////////////////////////////////
  if( parameters.verbose )
    printLogo(); 

  // Symmetry replicas
  //////////////////////////////////////////////////////////////////////
  vector <Matrix> myMatrices;
  if (parameters.useSymmetry) {
    // create replicas
    SymReplica myReplicas(parameters.infile_name);
    // use replicas for subsequent FIRST calculations
    parameters.infile_name.insert(parameters.infile_name.find(".pdb",0),"_sym");
    myReplicas.getMatrices(myMatrices); 
  }
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // FIRST routines
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  MolFramework structure;
  
  if( parameters.run_first ){
    
    parameters.mapFromTaskNameToStatus["Build molecular framework"] = "Running";
    outputStatus();
    readData( structure, parameters.infile_name );
    
    // output the PID when running on flexweb.
    if( parameters.flexweb ) 
      output_PID( structure.path_name );

    excludeSites( structure );
    setVdwRadii( structure );
    buildFramework( structure );
    
    parameters.mapFromTaskNameToStatus["Build molecular framework"] = "Complete";
    parameters.mapFromTaskNameToStatus["Analyze flexibility"] = "Running";
    outputStatus();
    if ( !parameters.inputGhosts) {
      analyzeFlexibilityAndOutput( structure );
    }
    else {
      cerr << "Skipping pebble game, using input ghost file g.in" << endl;
    }
    parameters.mapFromTaskNameToStatus["Analyze flexibility"] = "Complete";     
    outputStatus();
  }
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // TIMME routines
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  if( parameters.run_timme ){
    
    parameters.includeHydrogenbonds = false;
    parameters.includeHydrophobicTethers = false;
    
    parameters.mapFromTaskNameToStatus["Build molecular framework"] = "Complete";
    outputStatus();
    
    TIMME timme; 
    
    timme.rigidClusterDecomposition();
    
    timme.saveOutput();
  }
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // FRODA routines
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  if( parameters.runFRODA ){

    parameters.mapFromTaskNameToStatus["Initialize FRODA dynamics"] = "Running";
    outputStatus();
    
    MolFramework target;
    MolFramework restart;
    
    if( parameters.using_target ){
      
      if( parameters.verbose ) cout << " Loading target structure." << endl;
      readData( target, parameters.target_name );
      setVdwRadii( target );
    }
    
    if( parameters.using_restart ){
      
      if( parameters.verbose ) cout << " Loading restart structure." << endl;
      readData( restart, parameters.restart_name );
      setVdwRadii( restart );
      
    }

    structure.freeGrid();
    Froda *froda = new Froda();
    if (parameters.useSymmetry) {
      froda->setMatrices(&myMatrices);
    }
    froda->runFRODA( structure, target, restart );
    delete froda;
  }
  
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  // MORPHing rountine
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  if( parameters.runMorph ){
    
    readData( structure, parameters.infile_name );
    setVdwRadii( structure );
    buildFramework( structure );
    if ( !parameters.inputGhosts) {
      analyzeFlexibilityAndOutput( structure );
    }
    else {
      cerr << "Skipping pebble game, using input ghost file g.in" << endl;
    }
    
    MolFramework target;
    MolFramework restart;
    
    if( parameters.using_target ){
      if( parameters.verbose ) cout << " Loading target structure." << endl;
      readData( target, parameters.target_name );
      setVdwRadii( target );
    }
    
    // The following function require that the input structure are from 
    // the Yale Morph server (ala Gerstein's group). They will match up 
    // residue names and numbers, and then renumber atom numbers 
    // appropriately for targeting in FRODA.
    //////////////////////////////////////////////////////////////////////
    renumber_structures_from_morph_server( structure, target );
    
    Froda *froda = new Froda();
    froda->runFRODA( structure, target, restart );
    delete froda;
  }
  
  /* 
  //If we want to run froda2 completely outside of froda,
  // Here's where we could do it.
  if( parameters.froda2 ) {
    //preprocessing on the MolFramework to get rid
    //of the grid, and to prep the neighbor table
    structure.freeGrid();
    structure.stripNonCovalentNeighbors();
    
    Froda2Builder_FromMolFramework builder( structure, parameters );  
    RigidUnitSystem *rigidUnitSystem = builder.getRigidUnitSystem();
    
    //here we would create the PerturbRelaxCycle object,
    //using the PerturbRelaxCycleBuilder class.
    //Then we would run the simulation.
  }
  */
  
  
  return(0);
}
////////////////////////////////////////////////////////////////////////////////
