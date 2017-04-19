#include "Parameters.h"
#include <cstring>

////////////////////////////////////////////////////////////////////////////////
Parameters::Parameters(){
    
  verbose = 1;
  interactive = true;
  run_number = 0;
  energy_cutoff = -1.0;
  hphobe_fxn = 3;
  PH_cutoff = 0.5;
  PH_nucleic_cutoff = 0.15;
  exclude_water = 0;
  resolution_factor = 0.0;
  grid_length = 4.5;
  maxdist = 2.5;
  output_level = 1;
  
  numberOfFRODAwebOverlaySnapshots = 20;
  
  cutoff_SB_hyd_accpt_dist = 4.0;
  cutoff_SB_donor_accpt_dist = 5.0;
  cutoff_SB_donor_hyd_accpt_angle = 100.0;
  cutoff_SB_hyd_accpt_base_angle = 80.0;
  cutoff_HB_hyd_accpt_dist = 4.0;
  cutoff_HB_donor_accpt_dist = 5.0;
  cutoff_HB_donor_hyd_accpt_angle = 100.0;
  cutoff_SR_normal_plane_angle = 30.0;
  cutoff_SR_center_normal_angle = 40.0;
  cutoff_SR_distance = 5.5;
  cutoff_SRnucleic_distance = 3.55;
 
  rama = 0;
  hp_bars = 2;
  rs_bars = 3;
  hb_bars = 5;
  ud_bars = 5;
  cov_bars = 5;
  
  energyFxnA = false;
  run_query = false;
  bond_dilution = 0;
  lock_disulfide = true;
  read_user_defined_constraint_file = false;
  use_model_number = 0;
  alt_location_label = "A";
  output_bbg = false;
  use_unpruned_mean_coord_in_stripy_plot = false;
  correction_for_2_4 = true;
  
  found_alt_side_chain_loc = false;
  run_first = true;
  runFRODA = false;
  run_timme = false;
  publicationStyle = false;
	covalentTIMME = true;
  timmeScaleFactorA = 1.0; // by default - don't scale 
  
  xPeriod = 0.0f;
  yPeriod = 0.0f;
  zPeriod=  0.0f;
  
  useSpectralColoring = false;
  saveSigma = false;
  diluteWithRigidClusterFactory = false;
  timme_backbone = false;
	timme_complement = false;
  timme_nThNeighborForFlexibility = 3;
  timme_nThNeighborForPlasticity = 3;
  timme_bodyClusterSize = 1;
  
  timme_neighborCutoffDistance = 10.0; // based on comparison between O(n^2) vs O(n) algorithms for MD simulation of ADK
  timme_outputPrefix = "timme";
  timme_useAssembly = false;
	
  use_first_numbering = false;
  
  includeHydrogenbonds = true;
  includeHydrophobicTethers = true;
  
  keep_missing_bbg_vertices = false;
  
  min_output_cluster_size = 20;
  do_not_validate = false;
  using_target = false;
  using_restart = false;
  
  step_size = 0.1;
  output_freq = 100;
  total_conformations = 1000;
  
  stripeThickness = 7.;
  separationBetweenStripes = 0.;
  minimumRigidClusterSize = 10; // adjusted default based on ADK 1.0ns MD simulations
  useLinearTimme = true;
  outputNeighborList = false;
  useLinearCutoffScale = true;

  body = false; 
  dihedral = false;
  dtol = 0.5;
  dstep = 0.01;
  prob = 1.0;
  
  maxFitCycles = 100;
  bodyResponseEvery =31;
  atleast =2;
  
  fancy_target=0;
  ghost_debug=0;
  every_cycle=0;
  
  suppress_mc = 1.0;
  ghostTol = 0.125;
  ph_radius = 0.5;
  vdw_overlap = 0.85;
  
  mobileRC1 = false;
  centri = false;
  centri_maxr = 1E10;
  localCenter = false;
  localCenters.clear();
  
  make_raw_data_file = true;
  hbin = false;
  hbout = false;
  phin = false;
  phout = false;
  covin = false;
  covout = false;
  srin = false;
  srout = false;
  outputGhosts = false;
  inputGhosts = false;
  all_5_bars = false;
  energy_set_on_command_line = false;
  flexweb = false;
  runMorph = false;
  morphFrames = 0;
  skip_identify_hydrogen_bonds = false;
  skip_identify_hydrophobic_tethers = false;
  skip_identify_aromatics = false;
  
  polarGeometry = true;
  polarHRadius = 0.5;
  
  cylindrical = false; 
  cyl_low_z = 0.0;
  cyl_high_z= 0.0;
  
  local = false;
  lstep = 0.001;
  
  vdwNotice = 0.9;
  doAnneal = false;
  annealScale = 0.01;
  constantRate = false;
  acceptRate = 0.5;
  doChill = false;
  chillSteps = 1000;
  
  useSeed = false;
  seed = 0;
  
  propto = false;
  reportByRSMD = false;
  RMSDSpacing = 0.5;
  
  use_group_id = false;
  scoreMC = false;
  frodaMCScale = 1.0;

  doGetArea = false; 

  doGetPhobicRg = false;
 
  mcRg = false;
  RgWeight = 0;
  minimumRg = 0; 
  
  mcPhobicRg = false;
  phobicRgWeight = 0;
  minimumPhobicRg = 0; 
  
  mcSASA = false;
  SASAWeight = 0;
  minimumSASA = 0; 
  
  mcPhobicSASA = false;
  phobicSASAWeight = 0;
  minimumPhobicSASA = 0; 
  
  mcPolarSASA = false;
  polarSASAWeight = 0;
  
  doHBScoring = false;
  HBWellDepth = 8.0;
  SBWellDepth = 10.0;

  nosterics = false;
  report_resrmsd = false;
  
  chatty = false;
  
  bondnew = true;
  
  useCM = false;
  useCMforward = false;
  useCMreverse = false;
  CMforward = 1.0;
  CMreverse = 0.1;
  CMlower = 10;
  CMupper = 20;
  CMmax = 0.01;

  persist = false;
  nPersist = 10;
  
  findExposed = false; //check for surface exposure?
  exposedBonds = 2; //include exposed bonds?
  
  dilutionPlotXAxisType = perResidue;
  useHatbox = false;
  hatboxZa = -100;
  hatboxZb = 100;
  hatboxXY = 100;
  
  usePairs = false;
  stopWhenAllPairConstraintsSatisfied = false;
  
  useElasticVector = false;
  useInternalBias = false;
  
  useTalksto = false;
  changeWeights = false;
  
  // electron-density options
  useTheoMap = false;
  useEZDMap = false;
  trimMap = false;
  useED = false;
  edTol = 1.0;
  edResFac = 0.0;
  edNoise = 0.0;
  trimMapFactor = 1;
  // SAXS options
  useSAXS = false;
  saxsTol = 1.0;
  saxsFile = "";
  
  //Froda2 options
  froda2 = false;
  froda2Hybrid = false;
  froda2UseGhostTol = false;
  froda2Momentum = false;
  froda2Tol = 0.01;
  froda2RepulsionType = "froda";
  froda2AmberPrmtopFile = "";
  froda2AmberTrajOutput = false;

  stopAtEffectiveTimeLimit = false;
  effectiveTimeLimit = 10; // 10 ps
  effectiveTimeFreq = 1;   // 1 ps
}

////////////////////////////////////////////////////////////////////////////////
bool Parameters::isParameter( char *arg ){

  if( !strncmp( arg, "-", 1) )
    return true;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
bool Parameters::checkNonBoolParameter( int argNumber, int argc, char **argv, string myArg ){
 
  string currentArg;
  currentArg = (argv[argNumber])+1;

  if( currentArg != myArg )
    return false;

  if( argNumber+1 >= argc ){
    cout << " Error: A Missing value for command-line argument [" << myArg << "]" << endl;
    exit(1);
  }
  //else if( isParameter( argv[argNumber+1] ) ){
  //  cout << " Error: B Missing value for command-line argument [" << myArg << "]" << endl;
  //  exit(1);
  //}

  return true;
}

////////////////////////////////////////////////////////////////////////////////
bool Parameters::checkBoolParameter( int argNumber, int argc, char **argv, string myArg ){
 
  string currentArg;
  currentArg = (argv[argNumber])+1;
  
  if( currentArg != myArg )
    return false;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Read each white-space delimited string that was entered on the command line. 
//   If the string matches one of the known parameters, store the information. If 
//   an unknown flag is entered, exit. Each command-line option should be listed 
//   separately preceeded by a hyphen.
// Parameters:
//   argc - Total number of white-space delimited command-line arguments passed to
//          FIRST.
//   argv - character array of command line fields.
//   parameters - A data structure for storing the parameters used by FIRST.
////////////////////////////////////////////////////////////////////////////////
void Parameters::readCommandLineOptions( int argc, char **argv ){

  // Set the path to the FIRST distribution root. Required for finding 
  // the library files. If the user compiles the program, FIRST_ROOT is
  // set in the Makefile. In the "exe_only" version, FIRST_ROOT is NOT
  // set, and this program will look for an environment variable named
  // FIRST_ROOT. If no such variable exists, the program will exit. 
  //////////////////////////////////////////////////////////////////////

  #ifdef FIRST_ROOT
    path = FIRST_ROOT;
  #endif

  #ifndef FIRST_ROOT
    char *envPath;
    envPath = getenv("FIRST_ROOT");
    if( envPath != NULL )
      path = (string) envPath;
  #endif

  cout << "Reading command line options." << endl;
  bool fileFound = false;
  time_t current_time;
  time(&current_time);
  start_time = ctime(&current_time);

  command_line = argv[0];
  string space = " ";

  // Store the command line arguments. 
  //////////////////////////////////////////////////////////////////////
  for( int a = 1; a < argc; a++ ) {
    command_line += space + argv[a];
  }

  // Read all the command line arguments. 
  //////////////////////////////////////////////////////////////////////
  for( int a = 1; a < argc; a++ ) {

    // The options are split into sections, one for each of the three major
    // operations performed by FIRST, FRODA and TIMME. 
    //////////////////////////////////////////////////////////////////////
    if( isParameter(argv[a]) ){

      //////////////////////////////////////////////////////////////////////
      if( checkBoolParameter( a, argc, argv, "all5bars" ) ){
	all_5_bars = true;
      }

      /////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "anneal" ) ||
               checkNonBoolParameter( a, argc, argv, "targetMC" ) ){
	a++;
	doAnneal = true;
        if (argv[a][strlen(argv[a])-1]=='%') {
          constantRate = true;
          argv[a][strlen(argv[a])-1] = '\0';
          acceptRate = atof(argv[a])/100.0;
        }
        else {
          annealScale = atof( argv[a] );
        }
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter(a, argc, argv, "chill" ) ) {
        a++;
        doChill = true;
        chillSteps = atoi(argv[a]);
      }
      
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "atleast" ) ){
	a++;
	atleast = atoi( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////      
      else if( checkBoolParameter( a, argc, argv, "bbg_1sites" ) ){
	keep_missing_bbg_vertices = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "nosterics" ) ){
	nosterics = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "body" ) ){
	body = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "bondnew") ) {
        bondnew = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "bondold") ) {
        bondnew = false;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "bresponse" ) ){
	a++;
	bodyResponseEvery = atoi( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "centri_maxr") ){
	a++;
	centri_maxr = atof( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "centri") ){
	centri = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "changeweights" ) ){
	changeWeights = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "chatty" ) ){
	chatty = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "cyl" ) ){
	centri = true;
	cylindrical = true;
      } 
	
      //////////////////////////////////////////////////////////////////////		
      else if ( checkBoolParameter( a, argc, argv, "dust") ) {
	covalentTIMME = false;
	timme_bodyClusterSize = 0;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "lowz") ){
	a++;
	cyl_low_z = atof( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "highz" ) ){
	a++;
	cyl_high_z = atof( argv[a] );
      }

      // Print a list of covalent bonds
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "covout" ) ){
	covout = true;
      }

      // Read a list of covalent bonds
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "covin" ) ){
	covin = true;
      }

      // Alter the value of PH_cutoff used in the Hphobe tether fxn. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "c" ) ) {
	a++;
	PH_cutoff = atof( argv[a] );
	if( PH_cutoff < 0 ){
	  cout << " Error: The value for command-line argument " 
	       << argv[a-1] << " should be greater than 0.0." << endl;
	  exit(1);
	}
      }

      // apply conjugate motion, forward direction, to enhance motion
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "CMforward" ) ){
	a++;
	CMforward = atof( argv[a] );
	useCMforward = true;
	useCM = true;
      }

      // apply conjugate motion in reverse direction, to recover from jams
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "CMreverse" ) ){
	a++;
	CMreverse = atof( argv[a] );
	useCMreverse = true;
	useCM = true;
      }

      // conjugate-motion lower limit - enhance motion when fitting is good 
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "CMlower" ) ){
	a++;
	CMlower = atoi( argv[a] );
      }

      // conjugate-motion upper limit - reverse motion when fitting is bad  
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "CMupper" ) ){
	a++;
	CMupper = atoi( argv[a] );
      }

      // conjugate-motion upper limit - reverse motion when fitting is bad  
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "CMmax" ) ){
	a++;
	CMmax = atof( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "detailed" ) ){
	every_cycle = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "dihedral" ) ){
	dihedral = true;
      }

      // Option to use bond dilution.
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "dil" ) ){
	a++;
	bond_dilution = atoi( argv[a] );

	if( bond_dilution < 0 ||
	    bond_dilution > 2 ){
	  cout << " Error: Invalid bond dilution (-d) scheme selected." << endl;
	  exit(1);
	}
      }

      //
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "dstep" ) ){
	a++;
	dstep = atof( argv[a] );
      }
      
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "dtol" ) ){
	a++;
	dtol = atof( argv[a] );
      }
      
      // Change energy cutoff of H-bonds. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "E" ) ){
	a++;
	energy_cutoff = atof( argv[a] );
	energy_set_on_command_line = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "skiphbonds" ) ) {
	includeHydrogenbonds = false;
      } 		

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "skipphobes" ) ) {
	includeHydrophobicTethers = false;
      }

      //search for solvent-exposed atoms
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "findexposed" ) ) {
	findExposed = true;
      }

      //include hbonds and phobic tethers with no more than N exposed atoms
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "exposedbonds" ) ){
	a++;
	exposedBonds = atoi( argv[a] );
	findExposed = true;
      }
      
			//////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "complement" ) ){
				timme_complement = true;
      }
			
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "backbone" ) ){
	timme_backbone = true;
      }

      else if (checkBoolParameter(a, argc, argv, "numbers")) {
	use_first_numbering = true;
      }
			
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "neighbor" ) ) {
	a++;
	timme_nThNeighborForFlexibility = atoi( argv[a] );
	timme_nThNeighborForPlasticity  = atoi( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////			
      else if ( checkNonBoolParameter( a, argc, argv, "cutoffR") ) {
	a++;
	timme_neighborCutoffDistance = atof(argv[a]);
      } 
      
      //////////////////////////////////////////////////////////////////////
      else if (checkNonBoolParameter( a, argc, argv, "unitSize") ) {
	a++;
	timme_bodyClusterSize = atoi(argv[a]);
      }

      //////////////////////////////////////////////////////////////////////
      else if ( checkBoolParameter( a, argc, argv, "useAssembly") ) {
	timme_useAssembly = true;
      } 
				
      // Access the FRODA routines. 
      ////////////////////////////////////////////////////////////////////////////////	
      else if( checkBoolParameter( a, argc, argv, "FRODA" ) ) {
	runFRODA = true;
      }

      // perform area estimates in FRODA
      ////////////////////////////////////////////////////////////////////////////////	
      else if( checkBoolParameter( a, argc, argv, "getarea" ) ) {
	doGetArea = true;
      }

      // perform phobic-radius-of-gyration estimates in FRODA
      ////////////////////////////////////////////////////////////////////////////////	
      else if( checkBoolParameter( a, argc, argv, "getphobicrg" ) ) {
	doGetPhobicRg = true;
      }

      // Set scale for folding monte carlo. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scoreMC" ) ){
	scoreMC = true;
	a++;
	frodaMCScale = atof( argv[a] );
      }

      //Score hydrogen bonds using Mayo potential
      else if ( checkBoolParameter( a, argc, argv, "scoreHB" ) ){
        doHBScoring = true;
      }

      // Set radius of gyration scoring. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scoreRg" ) ) {
	a++;
	RgWeight = atof( argv[a] );
        mcRg = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scoreRgmin" ) ){
	a++;
	minimumRg = atof( argv[a] );
      }

      // Set phobic radius of gyration scoring. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scorephobicRg" ) ) {
	a++;
	phobicRgWeight = atof( argv[a] );
        mcPhobicRg = true;
	doGetPhobicRg = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scorephobicRgmin" ) ){
	a++;
	minimumPhobicRg = atof( argv[a] );
      }

      // Set solvent-exposed surface area (SASA) scoring. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scoreSASA" ) ) {
	a++;
	SASAWeight = atof( argv[a] );
        mcSASA = true;
	doGetArea = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scoreSASAmin" ) ){
	a++;
	minimumSASA = atof( argv[a] );
      }

      // Set phobic solvent-exposed surface area (SASA) scoring. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scorephobicSASA" ) ) {
	a++;
	phobicSASAWeight = atof( argv[a] );
        mcPhobicSASA = true;
	doGetArea = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scorephobicSASAmin" ) ){
	a++;
	minimumPhobicSASA = atof( argv[a] );
      }

      // Set polar solvent-exposed surface area (SASA) scoring. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "scorepolarSASA" ) ) {
	a++;
	polarSASAWeight = atof( argv[a] );
        mcPolarSASA = true;
	doGetArea = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "freq" ) ){
	a++;
	output_freq = atoi( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "fancy" ) ) {
	fancy_target = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "localcenter") ){
	localCenter = true;
        a++;
        localCenters.push_back( atoi( argv[a] ) );
        cerr << "Taking local center: " << localCenters.at( localCenters.size() -1 );
      }


      // Turn on local targeting
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "local" ) ) {
	local = true; 
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "lstep") ){
	a++;
	lstep = atof( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "flexweb" ) ) {
	flexweb = true;
      } 

      else if ( checkBoolParameter( a, argc, argv, "publication" ) ) {
        publicationStyle = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "gdebug" ) ){
	ghost_debug = true;
      }

      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "gtol" ) ){
	a++;
	ghostTol = atof( argv[a] );
      }

      // Select a different Hphobic tether identification fxn. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "H" ) ){
	a++;
	hphobe_fxn = atoi( argv[a] );
	
	if( hphobe_fxn == 1 ){
	  PH_cutoff = 0.25;
	  for( int ii = 1; ii < argc; ii++ ) { // If H is set to 1 check if c parameter exists. By default PH_cutoff=0.25 if H is set to 1
	    if( checkNonBoolParameter( ii, argc, argv, "c" ) ) {
	      ii++;
	      PH_cutoff = atof( argv[ii] );
	      if( PH_cutoff < 0 ){
		cout << " Error: The value for command-line argument " 
		     << argv[ii-1] << " should be greater than 0.0." << endl;
		exit(1);
	      }
	    }
	  }
	}

      }
     

      // Use hatbox constraint ( for membrane simulation).
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "usehatbox" ) ){
	useHatbox = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "hatboxza" ) ){
	a++;
	hatboxZa = atof( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "hatboxzb" ) ){
	a++;
	hatboxZb = atof( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "hatboxxy" ) ){
	a++;
	hatboxXY = atof( argv[a] );
      }

#ifndef BASE_ONLY
      // Options for dealing with electron density maps
      /////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "edres" ) ){
	a++;
        double resolution = atof( argv[a] );
	edResFac = 0.114040345081 * resolution * resolution;
        useED = true;
        // If you're curious about where that number came from, see 
        // Craig's notebook for 08/10/06 - 08/11/06
        // and 08/07/07 - 08/09/07
      }
	  
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "theomap" ) ){
	a++;
	theoMapFile = argv[a];
        useTheoMap = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "ezdmap" ) ){
	a++;
	ezdMapFile = argv[a];
        useEZDMap = true;
      } 

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "trimmap" ) ){
        a++;
        trimMapFactor = atoi(argv[a]);
	trimMap = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "edtol" ) ){
	a++;
	edTol = atof( argv[a] );
        // CJ TODO: this parameter doesn't actually work
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "ednoise" ) ) {
	a++;
        edNoise = atof( argv[a] );
      }

      // SAXS options
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "saxstarget" ) ) {
	a++;
        useSAXS = true;
        saxsFile = argv[a];
      }
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "saxsout" ) ) {
        useSAXS = true;
        saxsFile = "";
      }       
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "saxstol" ) ) {
	a++;
        saxsTol = atof( argv[a] );
      }   

      // Symmetry-enforcement options
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "sym" ) ) {
        useSymmetry = true;
      }  
#endif
 
      // Option to write hbonds to file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "hbout" ) ){
	hbout = true;
      }

      // Option to write ghosts to file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "gout" ) ){
	outputGhosts = true;
      }

      // Option to read ghosts in file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "gin" ) ){
	inputGhosts = true;
      }

      // Option to read hbonds from file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "hbin" ) ){
	hbin = true;
      }
      
      // Print usage
      //////////////////////////////////////////////////////////////////////
      else if( !strncmp( argv[a], "-help", 5 ) ){
	printUsage();
	exit(0);
      }

      // Read in location of FIRST root directory. This will override all
      // previously set "path" variables.
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "L" ) ){
	a++;
	path = argv[a];	

	if( path.size() == 0 ){
	  cout << " Error: No path to FIRST distribution was set." << endl;
	  cout << "        The -L flag was empty." << endl << endl;
	  exit(1);
	}
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "maxfitc" ) ){
	a++;
	int cycles = atoi( argv[a] );
	if( cycles > 500 ){
	  cout << " The maximum value for the command-line option -maxfitc is 500." << endl;
	  exit(1);
	  //cout << " The value option " << argv[a-1] << " has been reset to -maxfitc 500.";
	  //maxFitCycles = 500;
	}
	else{
	  maxFitCycles = cycles;
	}
      }

      // Set the minimum rigid cluster size to color
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "minRC" ) ){
	a++;
	min_output_cluster_size = atoi( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "mobRC1" ) ){
	mobileRC1 = 1;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "mode" ) ){
	useElasticVector = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if ( checkBoolParameter( a, argc, argv, "modei" ) ) {
	useElasticVector = true;
	useInternalBias = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "morph" ) ){
	a++;
	morphFrames = atoi( argv[a] );
	reportByRSMD = true;
	runMorph = true;
	run_first = false;
	runFRODA = false;
      }

      // When reading an NMR pdb file, explicitly choose the model number to
      // analyze. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "nmrmdl" ) ){
	a++;
	use_model_number = atoi( argv[a] );
      }                                                           

      // Do not factor in terms to make <r> more like glasses.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "nofactorR" ) ){
	correction_for_2_4 = false;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "nohbonds" ) ){
	skip_identify_hydrogen_bonds = true;
      }

      //
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "nohphobes" ) ){
	skip_identify_hydrophobic_tethers = true;
      }

      //
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "noaromatics" ) ) {
	skip_identify_aromatics = true;
      }

      // Use parameter defaults for all values, except those altered 
      // using command-line flags. User selection screens are turned off.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "non" ) ){
	interactive = false;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "noticevdw" ) ){
	a++;
	vdwNotice = atof( argv[a] );
      }

      // Do not run the structure validation routine if this flag is set.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "noval" ) ){
	do_not_validate = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "energyFxnA" ) ){
	energyFxnA = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "outbbg" ) ){
	output_bbg = true;
      }
      
      // Set the level of detail of textual output.
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "o" ) ){
	a++;
	output_level = atoi( argv[a] );
      }

      // Option to read pair associations from file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "pair" ) ){
	  usePairs = true;
      }

      else if( checkBoolParameter( a, argc, argv, "pairEarlyStop" ) ){
        stopWhenAllPairConstraintsSatisfied = true;
      }
      // 
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "ph_tol" ) ){
	a++;
	ph_radius = atof( argv[a] );
      }

      // Option to read hphobic tethers from file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "phout" ) ){
	phout = true;
      }

      // Option to read hphobic tethers from file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "phin" ) ){
	phin = true;
      }

      //options here for polar and phobic interactions
      
      // Option to deactivate polar geometry .
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "polarGeometry_off" ) ){
	polarGeometry = false;
      }

      // Option to set polar hydrogen radius.
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "polar_h_rad" ) ){
	a++;
	polarHRadius = atof( argv[a] );
      }

      // 
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "prob" ) ){
	a++;
	prob = atof( argv[a] );
      }

      ////////////////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "propto" ) ){
	propto = true;
      }

      // Open the query_network menus even in noninteractive mode.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "q" ) ){
	run_query = true;
      }

      // 
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "rama" ) ){
	rama = 1;
      }

      // Do NOT create the raw data file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "raw" ) ){
	make_raw_data_file = false;
      }

      // Report by rmsd not frequency
      else if( checkNonBoolParameter( a, argc, argv, "rep_by_rmsd" ) ){
	a++;
	RMSDSpacing = atof( argv[a] );
        reportByRSMD = true;
      }

      // Report rmsd by residue in output files
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "resrmsd" ) ){
	report_resrmsd = true;
      }

      // Read the filename of a restart structure for FRODA
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "restart" ) ){
	a++;
	restart_name = argv[a];	
	using_restart = true;

	if( restart_name.size() == 0 ){
	  cout << " Error: No restart file name given." << endl << endl;
	  cout << " USAGE: FIRST <original_input_file> -FRODA -restart<restart_file_name>." << endl << endl;
	  cout << "        The original input file is required in order to reproduce the rigid" << endl;
	  cout << "        cluster decomposition that is unique to this run." << endl;
	  exit(1);
	}
      }

      //////////////////////////////////////////////////////////////////////      
      else if( checkNonBoolParameter( a, argc, argv, "seed" ) ){
	a++;
	seed = atoi( argv[a] );
        useSeed = true;
      }

      // Option to output stacked ring interactions to file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "srout" ) ){
	srout = true;
      }

      // Option to read stacked ring interactions from file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "srin" ) ){
	srin = true;
      }   
      // 
      ////////////////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "suppress" ) ){
	a++;
	suppress_mc = atof( argv[a] );
      }

      // Option to allow rotation about disulfide bonds.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "S" ) ){
	lock_disulfide = false;
      }

      // Option to move single residue at a time.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "singlemove" ) ){
	  makeSingleMove = true;
      }

      // step size for FRODA moves 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "step" ) ){
	a++;
	step_size = atof( argv[a] );
      }

      // persist FRODA moves for several iterations
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "persist" ) ){
        persist = true;
	a++;
	nPersist = atoi( argv[a] );
      }

      // Read the filename of a target structure for directed dynamics in FRODA
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "target" ) ){
	a++;
	target_name = argv[a];
	using_target = true;
      }

      // set the minimumRigidClusterSize for dilution plots
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "minimumRigidClusterSize" ) ){
	a++;
        minimumRigidClusterSize = atoi( argv[a] );
      }

#ifndef BASE_ONLY      
      // Run the TIMME routines on a set of structures. 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "TIMME" ) ){
	run_timme = true;
	run_first = false;
      }

      // Color dilutions and decompositions with the new "spectrum" coloring scheme
      //////////////////////////////////////////////////////////////////////
      else if ( checkBoolParameter(a, argc, argv, "spectrum")) {
        useSpectralColoring = true;
      }

      //
      //////////////////////////////////////////////////////////////////////
      else if ( checkNonBoolParameter( a, argc, argv, "xPeriod") ) {
        a ++;
        xPeriod = atof(argv[a]);
      }

      //
      //////////////////////////////////////////////////////////////////////
      else if ( checkNonBoolParameter( a, argc, argv, "yPeriod") ) {
        a ++;
        yPeriod = atof(argv[a]);
      }
      
      //
      //////////////////////////////////////////////////////////////////////
      else if ( checkNonBoolParameter( a, argc, argv, "zPeriod") ) {
        a ++;
        zPeriod = atof(argv[a]);
      }

      // Run the TIMME routines on a set of structures. 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "saveSigma" ) ){
        saveSigma = true;
      }

      else if ( checkNonBoolParameter( a, argc, argv, "scaleFactor") ) {
        a ++;
        timmeScaleFactorA = atof(argv[a]);
      }
      
      // Use the (slow) exhaustive O(n^2) TIMME algorithm. 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "exhaustive" ) ){
        useLinearTimme = false;
      }

      // Use the nonlinear (per-change) cutoff scaling for dilution plots
      //////////////////////////////////////////////////////////////////////
      else if ( checkBoolParameter( a, argc, argv, "oldCutoffScale") ) {
        useLinearCutoffScale = false;
      }

      // Save the neighbor list used in TIMME analysis
      //////////////////////////////////////////////////////////////////////
      else if ( checkBoolParameter( a, argc, argv, "outputNeighborList") ) {
        outputNeighborList = true;
      }
      
      //perSite, perResidue, perRigidLabel
#endif      
      // dilution plot per-site 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "perSite" ) ){
        dilutionPlotXAxisType = perSite;
      }
      
      // dilution plot per-residue 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "perResidue" ) ){
        dilutionPlotXAxisType = perResidue;
      }
      
      // dilution plot per-residue 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "perRigidLabel" ) ){
        dilutionPlotXAxisType = perRigidLabel;
      }
      
      // Change stripeThickness. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "stripeThickness" ) ) { // FIXME - pick a better flag
	a++;
        stripeThickness = atof( argv[a] );
      }
      
      // Change separationBetweenStripes. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "separationBetweenStripes" ) ) {  // FIXME - pick a better flag
	a++;
        separationBetweenStripes = atof( argv[a] );
      }
      
      // TODO - remove this when development is complete or pick a better flag
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "NEWRCD") ) { 
        diluteWithRigidClusterFactory  = true;
      }
 
      // Change energy cutoff of H-bonds. 
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "totconf" ) ){
	a++;
	total_conformations = atoi( argv[a] );
      }

      // Read a list of user-defined constraints from a file.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "ud" ) ){
	read_user_defined_constraint_file = true;
      }

      // Use the unpruned value for <r> in the stripy plots. 
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "unprunedR" ) ){
	use_unpruned_mean_coord_in_stripy_plot = true;
      }

      // Enables reading/writing of the FIRST_group_ID field,
      // a FIRST-defined extra field in pdb files
      //////////////////////////////////////////////////////////////////////      
      else if( checkBoolParameter( a, argc, argv, "use_group_id" ) ){
        use_group_id = true;
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "vdw" ) ){
	a++;
	vdw_overlap = atof( argv[a]  );
      }

      // Option to turn off run-time messages.
      //////////////////////////////////////////////////////////////////////
      else if( checkNonBoolParameter( a, argc, argv, "v" ) ){ 
	a++;
	verbose = atoi( argv[a] );
      }

      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "froda2" ) ){
        froda2 = true;
      }
      else if( checkBoolParameter( a, argc, argv, "froda2Hybrid" ) ){
        froda2Hybrid = true;
      }
      else if( checkBoolParameter( a, argc, argv, "froda2UseGhostTol" ) ){
        froda2UseGhostTol = true;
      }
      else if( checkBoolParameter( a, argc, argv, "froda2Momentum" ) ){
        froda2Momentum = true;
      }
      else if( checkNonBoolParameter( a, argc, argv, "froda2Tol" ) ){
        a++;
        froda2Tol = atof( argv[a]  );
      }
      else if( checkNonBoolParameter( a, argc, argv, "froda2RepulsionType" ) ){
        a++;
        froda2RepulsionType = argv[a];
      }
      else if( checkNonBoolParameter( a, argc, argv, "froda2AmberPrmtopFile" ) ){
        a++;
        froda2AmberPrmtopFile = argv[a];
      }
      else if( checkBoolParameter( a, argc, argv, "froda2AmberTrajOutput" ) ){
        froda2AmberTrajOutput = true;
      }
      else if( checkNonBoolParameter( a, argc, argv, "effectiveTimeLimit" ) ){
        a++;
        effectiveTimeLimit = atof( argv[a]  );
        stopAtEffectiveTimeLimit = true;
      }
      else if( checkNonBoolParameter( a, argc, argv, "effectiveTimeFreq" ) ){
        a++;
        effectiveTimeFreq = atof( argv[a]  );
        stopAtEffectiveTimeLimit = true; 
      }
      
      // Option to turn off run-time messages.
      //////////////////////////////////////////////////////////////////////
      else if( checkBoolParameter( a, argc, argv, "version" ) ){

	cout << endl;
	cout << " Floppy Inclusions and Rigid Substructure Topography (FIRST) " 
	     << FIRST_VERSION_ID << endl;
	cout << " Framework Rigidity Optimized Dynamics Algorithm (FRODA) " 
	     << FRODA_VERSION_ID << endl;
	cout << " Tool for Identifying Mobility in Macromolecular Ensembles (TIMME) version " 
	     << TIMME_VERSION_ID << endl << endl;
	cout << " Copyright (c) 2004 - 2007 Brandon Hespenheide, Stephen Wells, Scott Menor." << endl << endl;
	cout << " Use the -help flag for more options." << endl << endl;
	exit(0);
      }

      // The current flag was not recognized.
      //////////////////////////////////////////////////////////////////////
      else{
	cout << " Error: Command line option " << argv[a] << " not recognized." << endl;
	exit(1);
      }
    }

    // Handle the arguments that don't begin with a "-" character. These 
    // should only be file names. 
    //////////////////////////////////////////////////////////////////////
    else{
      if( !fileFound ){
	infile_name = argv[a];
	fileFound = true;
      }
      else{
	cout << " Error: Too many input file names found. Error in command line." << endl;
	exit(1);
      }
    }
  }

  // Check to see if a file name was given.
  //////////////////////////////////////////////////////////////////////
  if( infile_name.size() == 0 ){
    cout << " Error: No valid file name was given." << endl;
    exit(1);
  }
    
  // Check to see if the given file name exists. 
  //////////////////////////////////////////////////////////////////////
  fstream test( infile_name.c_str() );
  if( !test ){
    cout << " Error: Could not find the file: [" << infile_name << "] in this directory." << endl;
    exit(1);
  }

  // Check that a path to the root directory has been set.
  //////////////////////////////////////////////////////////////////////
  if( path.size() == 0 ){
    cout << endl;
    cout << "ERROR: No path to the FIRST directory was found." << endl << endl;
    cout << "       FIRST_ROOT = [" << path << "]" << endl << endl;
    cout << "       FIRST needs to know the location of the directory where the external" << endl;
    cout << "       text files are found. If this is a binary distribution, you may:" << endl;
    cout << "       1. Set an envirnoment variable named FIRST_ROOT that points to the " << endl;
    cout << "          root directory of the FIRST distribution (FIRST/) " << endl;
    cout << "       OR " << endl;
    cout << "       2. You may set FIRST_ROOT on the command line by using the" << endl;
    cout << "          command-line flag -L <path-to-FIRST-root-directory>." << endl << endl;

    exit(1);
  }   

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Print the command-line options for the program.
////////////////////////////////////////////////////////////////////////////////
void Parameters::printUsage(){

  cout << endl 
       << " USAGE: FIRST [optional FLAGS] <filename> -target <file_name> -restart <filename>" << endl << endl;
  cout << " [optional FLAGS]: F = real number; I = integer. Do not include () in the argument." << endl << endl;
  cout << "    EXAMPLE: -c 0.53 sets the cutoff distance for identifying hydrophobic " << endl
       << "                    tethers to 0.53 Angstroms." << endl << endl;

  cout << " FIRST specific options (includes general options):" << endl;
  cout << " ----------------------------------------------------------------------" << endl;
  cout << " -bbg_1sites  Used when reading body-bar graph files. Forces the inclusion of DOF" << endl;
  cout << "                 for site numbers that are skipped in the file." << endl;
  cout << " -c (F) ----  Set the R_cutoff value to (F) Angstroms for identifying hydrophobic" << endl;
  cout << "                 tethers. Default varies by function (See -H flag)." << endl;
  cout << " -covout ---  Output a list of covalent bonds. The file will be named \"cov.out\"" << endl;
  cout << " -dil (I)---  Select which hydrogen bond dilution scheme to use. Default is -dil1." << endl;
  cout << " -E (F) ----  Set the energy cutoff to \"F\" for selecting which hydrogen bonds to include." << endl;
  cout << "                 The default energy cutoff is -1.0 kcal/mol." << endl;
  cout << " -FRODA ----  Use the program FRODA to explore the conformational space of the protein." << endl;
  cout << " -H (I)-----  Select the function to identify hydrophobic tethers." << endl;
  cout << "                 Ranges from 0(no hydrophobics) - 3(most restrictive function)." << endl
       << "                 The default is -H3." << endl;
  cout << " -hbin  ---   Read hydrogen bonds from a file. Check manual for file format." << endl;
  cout << " -hbout  --   Output a list of all hydrogen bonds identified by FIRST." << endl;
  cout << " -help   --   Print this usage screen." << endl;
  cout << " -minRC (I)   Set the minimum rigid cluster size for display." << endl; 
  cout << " -noaromatics Run without identifying stacked aromatic interactions." << endl;
  cout << " -nohbonds    Run without identifying hydrogen bonds." << endl;
  cout << " -nohphobes   Run without identifying hydrophobic tethers." <<endl;
  cout << " -non   ---   Run in noninteractive mode." << endl;
	cout << " -numbers -   Use inferred FIRST numbering for input / output files rather than " << endl;
	cout << "                 PDB site IDs. " << endl;
  cout << " -o (I) ---   Level of output sent to the *_results.txt file." << endl;
  cout << "                 The range is 0(no output detail) - 3(most detail). Default = 1." << endl;
  cout << " -phin  ---   Read hydrophobic tethers from a file. Check manual for file format." << endl;
  cout << " -phout  --   Output a list of hydrophobic tethers. List will depend on which hydrophobic" << endl;
  cout << "                 tether function was used (see -H flag)." << endl;
  cout << " -q   -----   Run the \"query network\" menu system even in noninteractive mode." << endl;
  cout << " -srin  ---   Read stacked aromatic rings from a file. Check manual for file format." << endl;
  cout << " -srout  --   Output a list of all stacked aromatic rings identified by FIRST." << endl;
  cout << " -S   -----   Lock disulfide bonds against rotation." << endl;
  cout << " -TIMME  --   Perform TIMME analysis on a set of structures. The structures should be" << endl;
  cout << "                 stored in a multi-model pdb file." << endl;
  cout << " -unprunedr   Do not prune side chains and dangling ends before computing the mean coordination" << endl;
  cout << " -v(I) ----   Verbosity level. Default is 1. Range is 0 - 3." << endl;
  cout << " -version     Print version information." << endl;
  cout << endl;
  
  cout << " FRODA specific options:" << endl;
  cout << " -anneal (F) - Use Monte-Carlo targetting with scale (F)." << endl;
  cout << "                  Accept probability during targetting is (delta-RMSD)/(scale)." << endl;
  cout << " -body  ----  During Monte Carlo step, make random moves of whol clusters, not individual atoms." << endl;
  cout << " -bresponse (I) - Body-response is a mechanism which helps the fitting of ghosts to atoms." << endl;
  cout << "                     It is applied every (I) cycles of fitting. Default value is 31." << endl;
  cout << " -centri ---  Run FRODA in \"centrifuge\" mode. Will bias each atom away from the geometric center" << endl;
  cout << "              of the molecule (with -mobRC1 option) or from the rigid core (default)." << endl;
  cout << " -changeweights -- uses an experimental algorithm which gives greater 'weight' to some atoms." << endl;
  cout << "             Atoms with serious mismatches are 'weighted' more heavily during fitting." << endl;
  cout << "             May improve stability, but use with caution." << endl;
  cout << " -CMforward (F) - Apply forward momentum to the motion of bodies." << endl;
  cout << "                     (F) is the fraction of momentum to apply; default is 0.2." << endl;
  cout << " -CMlower (I) - Momentum is applied when the system is fitting easily." << endl;
  cout << "                   This is when the system fits in less than (I) cycles. Default is 10." << endl;
  cout << " -CMmax (F) -- Maximum size of a motion due to momentum, in Angstroms. Default is 0.05." << endl;
  cout << " -cyl -- Activate \"cylinder centrifuge\", which biases atoms away from the z axis." << endl;
  cout << " -dihedral -- perturb structure by altering dihedral angles." << endl;
  cout << " -dstep (F) - During directed runs, this is the directed step size, in Angstroms." << endl;
  cout << "              Directed steps are performed during non-anneal targetting and \"centrifuge\" modes.:" << endl;
  cout << " -dtol (F) --  During directed runs, this is an acceptable RMSD, in Angstroms, between the" << endl;
  cout << "                current conformer and the target state. signalling the end of the run." << endl;
  cout << " -fancy  --  Switch on 'fancy' targetting, using PDB ID numbers, for partial targets." << endl;
  cout << " -freq (I) --  Every (I)th conformer is output as a pdb file. Default = 100." << endl;
  cout << " -getarea  -- FRODA will calculate the solvent-exposed surface area of the conformer." << endl;
  cout << "              This is an approximation based on a algorithm by Brooks and Guvench '03." << endl;
  cout << "              FRODA reports a total SASA, also phobic SASA and polar SASA." << endl;
  cout << " -getphobicrg -- calculates a phobic radius of gyration for every FRODA conformer." << endl;
  cout << "               This is calculated using the beta carbons of the following residues:" << endl;
  cout << "               ALA VAL LEU ILE MET PHE TRP TYR CYS  and is reported as REMARK Phobic Rg." << endl;
  cout << " -gin      -- read ghost information from a file g.in; do not run pebble game." << endl;
  cout << " -gout     -- output ghost information to file g.out." << endl;
  cout << " -gtol (F) --  Tolerance (Angstroms) when fitting atoms to ghost templates." << endl;
  cout << "              Default value is 0.125 Angstroms." << endl;
  cout << " -highz (F) -- The cylinder-centrifuge mode acts on atoms with Z coordinate less than this." << endl;
  cout << " -local -- Use \"local targetting\", which moves atoms towards their neighbors in the target structure." << endl;
  cout << " -lowz (F) -- The cylinder-centrifuge mode acts on atoms with Z coordinate greater than this." << endl;
  cout << " -lstep (F) -- Directed step size for \"local targetting\" moves." << endl;
  cout << " -maxfitc (I) -- Perform up to (I) cycles of fitting before a conformer search fails." << endl;
  cout << " -mobRC1 -- Makes rigid cluster 1 (the core) mobile like other clusters." << endl;
  cout << " -nosterics -- Deactivate steric contact routine (atoms are points)." << endl;
  cout << " -noticevdw (F) -- Steric contact begins at (R) times the sum of atomic radii." << endl;
  cout << "                  Default value is 0.9, calibrated to a 6-12 potential." << endl;
  cout << " -pair  -- read a list of paired atoms from the file pair.in." << endl;
  cout << "           Format is atom1 atom2 desiredDistance <label>; label is optional." << endl;
  cout << "           atom1 and atom2 will be biased towards each other during perturbations." << endl;
  cout << "           Label A will apply bias to only atom1; label B to atom 2." << endl; 
  cout << "           The bias size is given using -dstep." << endl;
  cout << " -persist (I) -- the random biases persist for (I) iterations before changing." << endl;
  cout << " -ph_tol (F) -- Distance allowance (Angstroms) on hydrophobic tethers. Atoms move freely within " << endl;
  cout << "               (F) Angstroms of contact; longer separations are forbidden. Default 0.5." << endl;
  cout << " -polar_h_rad (F) -- Radius (Angstroms) for a polar (OH,NH) hydrogen atom. Default is 0.25." << endl;
  cout << " -polarGeometry_off -- Steric routine takes no account of polar effects. Not recommended." << endl;
  cout << " -prob (F) --  Probability for an atom or body to be moved during a Monte Carlo step." << endl;
  cout << "                Values not in range 0-1 default to 1 (everything moves)" << endl;
  cout << " -propto -- Directed steps in targetting are proportional to distance to target." << endl;
  cout << " -rama --  Report information on backbond Phi/Psi angles as a REMARK in the pdb files." << endl;
  cout << "           Variable angles are sorted into the forbidden, generous, allowed and core regions." << endl;
  cout << " -resrmsd -- Report RMSD by residue for each pdb file in a resrmsd.txt file." << endl;
  cout << " -restart <filename> - Begin the dynamics from a previously saved conformation. The" << endl;
  cout << "                   name of the file should follow the flag." << endl;
  cout << " -scoreHB -- include strength of hydrogen bonds in Monte Carlo score." << endl;
  cout << "             Strength is calculated using Mayo potential as for FIRST." << endl;
  cout << " -scoreMC (F)  - activates Scoring Monte Carlo, with pseudotemperature (F)." << endl;
  cout << "     Conformers are scored and a Metropolis criterion used to accept/reject them." << endl;
  cout << "     The convention is that lower scores are better." << endl;
  cout << " -scoreRg (F) - score conformers based on radius of gyration." << endl;
  cout << " -scoreRgmin (F) - sets a minimum threshold for radius of gyration scoring." << endl;
  cout << "  The Rg score is 0 if Rg < scoreRgmin; else it is given by:" << endl;
  cout << "                  scoreRg * ( Rg - scoreRgmin )." << endl;
  cout << " -scorephobicRg (F) - score conformers based on phobic radius of gyration." << endl;
  cout << " -scorephobicRgmin (F) - sets a minimum threshold for phobic radius of gyration scoring." << endl;
  cout << "  The phobicRg score is 0 if phobicRg < scorephobicRgmin; else it is given by:" << endl;
  cout << "                  scorephobicRg * ( phobicRg - scorephobicRgmin )." << endl;
  cout << " -seed (I) -- Initialises the random number generator using this seed." << endl;
  cout << "             If no seed is given, the system time at start of execution is used." << endl;
  cout << " -singlemove -- only one atom or cluster will be perturbed at a time." << endl;
  cout << " -step (F) -- Stepsize for random motion, in Angstroms." << endl;
  cout << " -target <filename> - Perform a targeted dynamics FRODA run from the input structure" << endl;
  cout << "                   to this structure. No spaces between the flag and the filename." << endl;
  cout << " -targetMC (F) - equivalent to -anneal, activates Monte Carlo targeting with scale F." << endl;
  cout << " -totconf (I)  Total number of conformations to produce. Default = 1000. (NOTE: Not all" << endl;
  cout << "                   conformations are written to file. Writing to file is set by the -freq flag." << endl;
  cout << " -vdw (F) -- steric clashes beyond this fraction of the sum of radii are forbidden." << endl;
  cout << "            Default is 0.85." << endl;
  cout << endl;

  cout << " TIMME specific options:" << endl;   // TODO - advertise the following switches 

  cout << " ----------------------------------------------------------------------" << endl;
  cout << " -neighbor (I) Use n'th nearest neighbors to calculate empirical bond network flexibility" << endl;
  cout << "              for an ensemble. Default is 2." << endl;
	cout << " -dust        Treat particles as a dust for TIMME analysis" << endl;
	cout << " -exhaustive  Perform an exhaustive O(n^2) TIMME analysis" << endl;
	cout << " -cutoffR     The maximum separation between pairs of sites included in a linear TIMME run" << endl;
	cout << " -saveSigma   output the sigma matrix (from a TIMME run) in Matrix Market exchange format" << endl;

		

/*cout << "   -backbone- Perform TIMME analysis on backbone atoms only (faster but incomplete; " << endl;
  cout << "              default is to perform analysis on all atoms)" << endl;
*/

  cout << endl;

  exit(0);
}
////////////////////////////////////////////////////////////////////////////////
