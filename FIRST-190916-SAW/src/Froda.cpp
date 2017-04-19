#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "Froda.h"
#include "mt19937ar.h"
#include "flexweb.h"
#include "Timer.h"
#include "Sterics.h"
#include "EDStructure.h"
#include "SAXS.h"
#include "GenericMap.h"
#include "AreaProximity.h"

#include "Froda2Builder_FromMolFramework.h"
#include "RigidUnitSystem.h"
#include "PerturbRelaxCycle.h"
#include "ConstraintEnforcingPotential.h"
#include "Mismatch.h"
#include "Timer.h"
#include "PhaseSpacePathLengthIntegrator.h"
#include "FitPerturber.h"
#include "MomentumPerturber.h"
#include "SymmetricPerturber.h"
#include "GlobalMotionRemover.h"
#include "RandomCenterPerturber.h"
#include "RandomRotorPerturber.h"
#include "MinimizeSystem.h"
#include "OutputCurrentStatus.h"
#include "SymmetryMatrices.h"

#include <limits>

#include <signal.h>

extern Parameters parameters;

const Vector Froda::NULL_VEC = Vector(0,0,0);
const double Froda::RAD_TO_DEG = 180.0/PI;
const double Froda::TWOPI = 2*PI;

////////////////////////////////////////////////////////////////////////////////
Froda::Froda() {
  
  whenTargetFound=0;
  debug_vdw = false; debug_rigid=false; debug_connect=false; debug_phobe=false; debug_hbond=false; debug_ghost=false;
  chatty_fitting =false; // report every cycle
  isBad=false;
  isRestart =0; isTargeted=0; foundTarget=0; // for targeted dynamics routines
  moveProb = 1.0;
  
  doAnneal= false; // for simulated annealing routine
  
  cylindrical = false; // for cylindrical centrifuge for membrane proteins
  
  localTargeting = false; // target locally
  
  /* tolerances during fitting */
  ghostTol = 0.125; // mismatch to ghost positions, Angstroms
  phobeTol = 0.5; // tolerance on hydrophobic variation
  epsilon = 0.0001; // very small vector length, equivalent to zero
  
  verysmall = 1E-5;
  
  regridEvery =2, regridCounter=0; // update vdW coarse gridding.
  
  maxFailsInRow = 9; // stop when jammed
  
  didIFail = false;
  failScore = 0;
  failLimit = 0;
  
  centrifugal=0; // if yes, bias the random steps to unfold the protein
    
  frameNumber = 0;
  
  cacheEvery = 3;
  maxRotorCycles = 100;
  fancy_targetting =0; // check original pdb numbers for partial target files
  
  nConfsFound =0;
  
  smallgrad=1E-5; large_rotor_2=0.5;
  
  n_myHbonds=0; n_myPhobes=0;
  
  //these variable are used to count VDW overlaps and sort them by severity;
  //e.g. an overlap by 0.16A goes into the bin (overlap/overlapStep)
  overlapStep = 0.1;
  
  //per-residue list: get all the calphas

  nFailedCycles = 0;
  
  checkRama =0;
  
  n_accepted=0; n_rejected=0;
 
  sterics = NULL;
  areaProximity = NULL;

  checkingWeights = false;
  running_allmain_MSD = 0.0;
  froda2Hybrid = parameters.froda2Hybrid;

}

////////////////////////////////////////////////////////////////////////////////
Froda::~Froda() {

  delete allSiteSelector;
  delete backboneSiteSelector;
  delete sidechainSiteSelector;
  delete targetedSiteSelector;
  delete targetedBackboneSiteSelector;
  if (parameters.useSymmetry) {
    delete symMat;
  }
}



////////////////////////////////////////////////////////////////////////////////
// Description: HERE BEGINNETH THE FRODA
///////////////////////////////////////////////////////////////////////////////
void Froda::runFRODA( MolFramework &structure, MolFramework &target, MolFramework &restart ){

  ofstream problemLog( "problem.log" );

  // setup site selectors
  allSiteSelector = new AllSiteSelector();
  backboneSiteSelector = new BackboneSiteSelector(&structure);
  sidechainSiteSelector = new SidechainSiteSelector(&structure);

  //cerr << setprecision(20);
  
  targetedSiteSelector = new TargetedSiteSelector(this);
  
  CompositeSiteSelector * targetedBackboneSiteSelector = new CompositeSiteSelector();
  targetedBackboneSiteSelector->addSiteSelector(backboneSiteSelector);
  targetedBackboneSiteSelector->addSiteSelector(targetedSiteSelector);
  this->targetedBackboneSiteSelector = targetedBackboneSiteSelector;

  preferredSiteSelector = targetedBackboneSiteSelector; // FIXME - make this a command-line option
  
  if( parameters.verbose ){
    cout << "Initializing FRODA... " << endl;
  }

  ofstream *runningRMSDfile = NULL;

  if (parameters.flexweb) {
    if (structure.total_sites > 99999)
      {
	cout << "Total number of atoms exceeds 100,000. FRODA stopped... " << endl;
	problemLog << "Total number of atoms exceeds 100,000. FRODA stopped..." << endl;
	exit(1);
      }
    runningRMSDfile = new ofstream();
    string runningRMSDfileName = structure.path_name;
    runningRMSDfileName += "RunningRMSD.txt";

    std::cout << runningRMSDfileName << std::endl;
    runningRMSDfile->open(runningRMSDfileName.c_str());

    *runningRMSDfile << "Iteration # \t Current \t Cumulative" << endl;

    *runningRMSDfile << fixed << setw(11) << 0 << " \t " << setw(7) << 0 << " \t " << setw(10) << 0 << endl;
  }
  
  // capture kill signal(s) so we can clean up gracefully  
  // TODO - factor these out into Main (or someplace more appropriate)
  signal(SIGINT, &Froda::killFRODA);
  signal(SIGTERM, &Froda::killFRODA);
  signal(SIGQUIT, &Froda::killFRODA);
 
  if (parameters.flexweb) {
    string file_name = structure.path_name;
    file_name += ".pdbSnapshotList";
    ofstream currentSnapshotFile( file_name.c_str() );
    currentSnapshotFile << structure.base_name << "_RCD.pdb";
    currentSnapshotFile.close();
  }
  

  obtainParameters( structure ); //read values from parameter list into local variables

  if (isRestart && verbose) cout << "RESTART RUN." << endl;
  if (isTargeted && verbose) cout << "TARGETED RUN." << endl;
  if (centrifugal && verbose) cout << "CENTRIFUGAL RUN." << endl;

  bool usesterics = !nosterics;

  //initialise the random number generator
  if ( !useSeed) {
    timeval *currentTime = new timeval();
    gettimeofday(currentTime, NULL);

    seed = currentTime->tv_usec;
  }

  //SAW 2 DAN: presumably we only need to initialise the RNG once, we can both use it.

  init_genrand(static_cast<unsigned long>(seed));

  if( parameters.verbose >= 2 ){
    cout << "  Using random seed: " << seed << endl;
  }

  //establish the central point of the protein for reference
  myCentralPoint.x = (structure.max_coords[X] + structure.min_coords[X])*0.5;
  myCentralPoint.y = (structure.max_coords[Y] + structure.min_coords[Y])*0.5;
  myCentralPoint.z = (structure.max_coords[Z] + structure.min_coords[Z])*0.5;

  nTotalSites = structure.total_sites;

  initialiseMainArrays(structure, restart);

  //SAW 2 DAN: somewhere around here, I guess we initialise your objects
  //Note that initialisePhobic and stripNonCovalentNeighbors modify the structure neighbor_list

  if ( !parameters.inputGhosts ) {
    initialisePhobic( structure );
    stripNonCovalentNeighbors( structure );
    getGhostMembershipFromRCD( structure );
    buildGhostGeometry( structure );
  }
  else {
    getGhostMembershipFromFile( structure );
    initialisePhobic( structure );
    stripNonCovalentNeighbors( structure );
    buildGhostGeometry( structure );
    getRCDFromGhosts( structure );
    structure.outputRCD_PDBFormat( "_RCD.pdb" );
    structure.outputPyMolScriptRCD_PDBFormat();
    structure.outputText();
  }
  nTotalClusters = structure.total_clusters;

  //SAW 2 DAN: by this point the FRODA ghosts are made.
  //For your engine, I guess we use the same information to initialise your ghosts.
  
  const double s_over_t_Carbon_Angstroms_per_ps = 7.894;
  PhaseSpacePathLengthIntegrator *pathLengthIntegrator = NULL;
  if ( parameters.stopAtEffectiveTimeLimit ) {
    // here we must use the FIRST style indexing, where index 0
    // is empty.
    size_t N = structure.total_sites + 1;
    vector<char> boolMask( N );
    boolMask[0] = 0;
    for ( size_t i = 1; i < N; i++ ) {
      boolMask[i] = structure.site_info[i].element_name == "C ";
    }
    pathLengthIntegrator = new PhaseSpacePathLengthIntegrator( &currentPos );
    pathLengthIntegrator->activateMasking( boolMask );
  }
  
  ////////Froda2 initialization
  vector<Vec3> froda2_MeanPositions;
  RigidUnitSystem *rigidUnitSystem = NULL;
  ConstraintEnforcingPotential *constraintEnforcingPotential = NULL;
  MinimizeSystem *minim = NULL;
  PerturbRelaxCycle *cycle = NULL;
  if ( froda2Hybrid ) {
    cout << "\n************" << endl;
    cout << " Using hybrid Froda/Froda2 " << endl;
    
    // turn off Froda's sterics.
    // (they are still activated inside Froda2) 
    usesterics = false;
    nosterics = true;
    
    froda2_MeanPositions.resize( currentPos.size() - 1 );
    
    Froda2Builder_FromMolFramework builder( structure, parameters, symMat );  
    rigidUnitSystem = builder.getRigidUnitSystem();
    constraintEnforcingPotential = builder.getConstraintEnforcingPotential();

    minim = new MinimizeSystem( rigidUnitSystem, constraintEnforcingPotential );

    //tolerance condition
    if ( parameters.froda2UseGhostTol ) {
      minim->setToleranceCondition( "mismatch", ghostTol );
      cout << " tolerance condition: mismatch\n";
      cout << " tolerance value: " << ghostTol << " Angstroms" << endl;
    }
    else {
      minim->setToleranceCondition( "maxPreconditionedGradComponent", parameters.froda2Tol );
      cout << " tolerance condition: estimated distance of each degree of freedom \n";
      cout << "                      from its nearby minimum\n";
      cout << " tolerance value: " << parameters.froda2Tol << " Angstroms" << endl;
    }
    
    //minimizer
    int NminimizationSteps = 200;
    cout << " maximum number of conj grad minimization steps: " << NminimizationSteps << endl;
    minim->setNminimizationSteps( NminimizationSteps );
    
    cycle = new PerturbRelaxCycle( minim );

    // Here we check to see if any perturbations have been requested
    // that do not already exist in Froda2.  If so, then we will
    // do the perturbation of atoms in Froda, and use Froda2's 
    // FitPerturber to translate this perturbation into a Froda2
    // rigid unit perturbation.
    // NOTE: It is important that the FitPerturber be set
    // before any of the other perturbers.  Otherwise
    // the effect of the FitPerturber would be to undo
    // any previously applied perturbations.
    if ( isTargeted && !doAnneal || 
         isTargeted && localTargeting ||
         centrifugal || 
         usePairs || 
         useElasticVector ) {
      FitPerturber *fitpert = new FitPerturber( rigidUnitSystem );
      fitpert->setMeanPointsTarget( &froda2_MeanPositions );
      cycle->addCommand_perturb( fitpert, &FitPerturber::perturb );
      cout << " enabled Fit to Froda Perturbation " << endl;
    }
    
    // Here, we check if any perturbations have been requested
    // that DO exist in Froda2.  If so, we will enable the perturbation
    // in Froda2.
    
    //momentum perturbation
    if ( useCM ) {
      cout << "Froda Conjugate Momentum is not compatible with Froda2Hybrid.\n";
      cout << "Do not use any Froda CM options with Froda2Hybrid.\n";
      cout << "You can activate the momentum perturbation using -froda2Momentum." << endl;
      exit(0);
    }
    if ( parameters.froda2Momentum ) {
      GlobalMotionRemover *globalMotionRemover = new GlobalMotionRemover( rigidUnitSystem );
      cycle->addCommand_cycleStart( globalMotionRemover, &GlobalMotionRemover::setCurrentPointsAsTarget );
      cycle->addCommand_cycleEnd( globalMotionRemover, &GlobalMotionRemover::fitCurrentPointsToTarget );    
      MomentumPerturber *mom = new MomentumPerturber( rigidUnitSystem );
      cycle->addCommand_perturb( mom, &MomentumPerturber::perturb );
      cycle->addCommand_cycleStart( mom, &MomentumPerturber::setQ1 );
      cycle->addCommand_cycleEnd( mom, &MomentumPerturber::determineDeltaQ );    
      cout << " enabled Momentum Perturbation" << endl;
    }
    
    //random perturbation
    if ( parameters.step_size > numeric_limits<double>::epsilon() && !parameters.useSymmetry ) {
      RandomCenterPerturber *pertC = new RandomCenterPerturber( rigidUnitSystem, parameters.step_size );
      cycle->addCommand_perturb( pertC, &RandomCenterPerturber::perturb );
      cout << " enabled random perturbation of rigid unit centers: \n" << 
              "    " << parameters.step_size << " Angstroms" << endl;
      RandomRotorPerturber *pertR = new RandomRotorPerturber( rigidUnitSystem, parameters.step_size  );
      cycle->addCommand_perturb( pertR, &RandomRotorPerturber::perturb );
      cout << " enabled random rotation of rigid unit centers: \n" << 
              "    " << parameters.step_size << " Angstroms on outermost atoms" << endl;
      //There is no Froda option that disables random perturbations.
      //So, to disable Froda's random perturbations, I (Dan) have put Froda's
      //call of the random perturbation routine in an "if" statement.
    }
    // symmetric perturbation
    if (parameters.useSymmetry && parameters.step_size > numeric_limits<double>::epsilon()) {
      SymmetricPerturber *pertSym = new SymmetricPerturber(rigidUnitSystem,symMat,parameters.step_size);
      cycle->addCommand_perturb(pertSym, &SymmetricPerturber::perturb);
      cout << " enabled symmetric perturbation:\n" << 
              "    " << parameters.step_size << " Angstroms" << endl;
    }
    //output enable
    OutputCurrentStatus *outputCurrentStatus = 
      new OutputCurrentStatus( rigidUnitSystem, cycle, minim, constraintEnforcingPotential );
    outputCurrentStatus->setPeriod( 10 );
    
    if ( parameters.stopAtEffectiveTimeLimit ) {
      cycle->addCommand_cycleEnd( pathLengthIntegrator, &PhaseSpacePathLengthIntegrator::updatePath );
      outputCurrentStatus->enablePathLength( pathLengthIntegrator, 1.0/s_over_t_Carbon_Angstroms_per_ps );
    }
    
    cycle->addCommand_cycleEnd( outputCurrentStatus, &OutputCurrentStatus::display );
    
    cout << "************\n" << endl;
  }
  ////////End Froda2 initialization  

  defineAtomsAsMainOrSide( structure);

  //report some atom counts
  if( parameters.verbose >= 2 ){ 
    cout << "  Total Main-Chain Atoms: " << nMainchainAtoms << " (" << nMobileMainchainAtoms << " are mobile)" << endl;
    cout << "  Total Side-Chain Atoms: " << nSidechainAtoms << " (" << nMobileSidechainAtoms << " are mobile)" << endl;

    cout << "  " << nMobileAtoms << " mobile sites within ";
    if (mobileRC1){
      cout << myNClusters << " mobile rigid clusters." << endl ;
    }
    else{
      cout << myNClusters-1 << " mobile rigid clusters." << endl ;
    }
  }
  
  defineClustersAsMainOrSide();

  getGeneralRamaPlot();
  getProlineRamaPlot();
  getGlycineRamaPlot();
  makeBackboneList( structure ); //note that this routine labels the alpha carbons
  makeRamaList( structure );

  //cerr << "Left Rama routines OK." << endl;

  //SAW 25April06 Removed last traces of "talksto"
  //are we using mode biasing?
  if ( useElasticVector ) {
    makeResidueList( structure );
    makeResidueBasis();

    readElasticNetworkMode( structure );
    formBiasVectors();
  }


  //cerr << "setting up arrays for RMSD by residue." << endl;
  setUpRMSDByResidueArrays();

  // Dan made this Sterics initialization conditional on the
  // already existing usesterics flag.  Seems to me that
  // if usesterics is off, we don't need to create the object  
  if ( usesterics ) {
    sterics = new Sterics( this, &structure );
  }
  //SAW 2 DAN: do your objects use a different kind of sterics object, or can we share this one?


  identifyPolarAtoms(structure);
  identifyPhobicAtoms(structure);

  if( doGetPhobicRg ) { 
    findPhobicBetaCarbons( structure );
  }
  
  findDonorHydrogens( structure );
  findHBondingAtoms( structure );
 
  findNeighbors(structure);

  // BMH commented out. Need better memory management.
  if ( doGetArea ) {
    areaProximity = new AreaProximity( this, &structure );
    obtainSurfaceParameters( structure );
    getAreaOfConformer();
  }
  
  //SAW Oct 06 removing menu option- command lines only please!

  if (moveProb < 0.0) moveProb = 1.0;
  if (moveProb > 1.0) moveProb = 1.0;
  //catch nonphysical probabilities
  //SAW 2 DAN: may be relevant to your perturb system

  
  if (isTargeted) setUpTargetList(structure, target);

  //if it's local-targetting, find local friends:
  if ( isTargeted && localTargeting) setUpLocalTargeting( structure, target);


  cachedPos.resize(nTotalSites+1);
  
  
  if( parameters.verbose ){
    cout << "Simulating " << nConfigs << " configurations." << endl;
    if( parameters.output_freq == 1 )
      cout << "Outputing every conformation." << endl;
    else if( parameters.output_freq == 2 )
      cout << "Outputing every second conformation." << endl;
    else if( parameters.output_freq == 3 )
      cout << "Outputing every third conformation." << endl;
    else
      cout << "Outputing every " << parameters.output_freq << "th conformation." << endl;
  }
  
  outputCounter = 0;

  //check global rmsd if needed for morph frames or report-by-rmsd
  if (runMorph){
    distanceToTarget(restart, target);
    thisRMSD = lastRMSD = targetRMSD;
    RMSDSpacing = thisRMSD / ( 1.0+morphFrames);
  } else if (reportByRSMD){
    if (isTargeted){
      distanceToTarget(restart, target);
      thisRMSD = lastRMSD = targetRMSD;
    } else { thisRMSD = lastRMSD = 0.0;}
  }

  parameters.mapFromTaskNameToStatus["Initialize FRODA dynamics"] = "Complete";
  parameters.mapFromTaskNameToStatus["Generate conformations"] = "Running";
  outputStatus();


  //SAW 2 DAN : made the ED and scoreMC sections clearer by moving some logic into subroutines
    // set up electron-density specific objects
  GenericMap *targetMap = NULL;  // instantiate these if they're needed
  EDStructure *myPdb = NULL;
  if (useED && useTheoMap) {
    setupTheoEDObjects( structure, myPdb, targetMap );
  }
  if (useED && useEZDMap) {
    setupRealEDObjects( structure, myPdb, targetMap );
  }
  if (useED && trimMap) {
    setupEDTrimMap( structure, myPdb, targetMap );
  }
  // SAXS setup code begins
  SAXS *mySAXS = NULL; // instantiate if needed
  if (useSAXS) {
    mySAXS = new SAXS(*this,structure,edResFac,saxsFile);
    if (saxsFile != "") { // simulation is targeted to SAXS profile
      mySAXS->cachedCorrelation = mySAXS->correlate();
      cerr << "Initial correlation to SAXS target: " << mySAXS->cachedCorrelation << endl;
    }
  }
  // end SAXS setup code
  if (doAnneal) {
    mcRejected = 0;
  }
  //SCORE the initial conformer here
  if (isTargeted && doAnneal){
    distanceToTarget( restart, target);
    cachedRMSD = targetRMSD;   
  } 
   
  if (doScoringMC){
    clearMCFoldingVariables();

    getMCScore( structure );

    //keep these values for later
    storeMCFoldingVariables();
  }

  if ( usePairs ) {
    readPairFile( structure );
  }

  //SAW 2 DAN: this routine checks if the input structure has steric problems
  //May need checking if we're using your sterics system?
  bool initialBadSterics = false;
  if ( usesterics && stericMismatch( structure ) > 0 ) {
    initialBadSterics = true;
  }
  clearMismatchArrays();
  if ( initialBadSterics ) {
    cerr << "WARNING: some bad sterics in input structure." << endl;
  }
  
  bool isLethal = false;  
  maxFailsInRow = 9;
  failsInRow = 0;
  failLimit = 2*maxFailsInRow;

  ////////Froda2 initial minimization
  if ( froda2Hybrid ) {
    cout << "\n************" << endl;
    cout << "Froda2:" << endl;
    cout << "Listing of initial structure constraint violations (if any),\n";
    cout << "  listed by internal froda2 array indices 0..N-1: " << endl;
    constraintEnforcingPotential->mismatch()->setVerbose( true );
    constraintEnforcingPotential->mismatch()->calc();
    constraintEnforcingPotential->mismatch()->setVerbose( false );
  
    cout << "\nPerforming initial minimization...";
    minim->minimize();
    cout << "Done.\n" << endl;
  
    cout << "Post-relaxation constraint violations (if any): " << endl;
    constraintEnforcingPotential->mismatch()->setVerbose( true );
    constraintEnforcingPotential->mismatch()->calc();
    constraintEnforcingPotential->mismatch()->setVerbose( false );
    cout << "************\n" << endl;
    cout << "ConformationNumber|NumberOfConjGradSteps|MaxMismatch|RMSDfromInitial" << endl;
  }
  ////////End Froda2 initial minimization


  persistCounter = 0;
  
  timer.start();

  for(i=0;i<nConfigs;i++) { // MAIN PRODUCTION LOOP BEGINS HERE
    if ( stopWhenAllPairConstraintsSatisfied && allPairConstraintsSatisfied() ) break;
    if ( parameters.stopAtEffectiveTimeLimit && 
        pathLengthIntegrator->getIntegratedPathLength() / s_over_t_Carbon_Angstroms_per_ps >= parameters.effectiveTimeLimit )
      break;
    // cerr << "Making config " << i+1 << endl;

    checkingWeights = false;
    if (parameters.changeWeights) {
      resetWeights();
    }
    // output the current status 
    updateStatus();
    
    outputCounter++;
    keepTrying = true;

    //save previous (good) position for recovery   
    cachedPos = currentPos;
    //save variables
    if (isTargeted && doAnneal){
      cachedRMSD = targetRMSD;
    }
    if (doScoringMC){
      storeMCFoldingVariables(); 
    }

    if ( useElasticVector ) {
      makeResidueBasis();
    }


    //SAW 2 DAN: this section is my random perturbation routine
    //In your system, presumable we instead call your random perturb routine
    //and come back with new atom positions

    //Dan's comments: If running froda2 hybrid, then we will not
    //use Froda's random perturbation.  Instead, the random
    //perturbation will be performed directly by Froda2
    if ( !froda2Hybrid ) {
      if (persist ) {
        if (persistCounter >= nPersist) {
          persistCounter = 0;
        }
        if ( persistCounter != 0 ) {
          oldRandomMove();
        }
        else {
          newPerturbation(structure);
        }
        persistCounter++;
      }
      else {
        newPerturbation(structure);
      }
    }
    
    //apply targeted moves
    //SAW 2 DAN: this section is all the atom-based biases
    //my routines act only on atom positions (currentPos) 
    if (isTargeted && !doAnneal) moveToTarget(target);
    if (isTargeted && localTargeting) moveToLocalTarget();
    if (centrifugal) applyCentrifuge();
    if ( usePairs ) applyPairs();
    if ( useElasticVector ) addModeBias();
    // apply CM moves based on last turn's cyclesToFit
    if ( useCM && i > 1){
      applyCMbias();
    }

    //SAW 2 DAN: at this point the perturbations are all done
    //if we're using your system, this is where we refit ghosts to atom positions

    //Dan's comments: 
    //Here we do the following:
    //  Take the Froda perturbed atom positions
    //     and send them into Froda2.
    //  Perform a full Froda2 cycle (perturbation/relaxation)
    //  Retrieve the new atom positions from Froda2 for use
    //     by Froda.
    if ( froda2Hybrid ) {
      // load in the perturbed atom positions
      for ( SiteID s = 1; s <= nTotalSites; s++ ) {
        froda2_MeanPositions[s-1] = currentPos[s];
      }
      
      // generate conformer
      cycle->doCycle();
      
      // retrieve the new positions
      froda2_MeanPositions = rigidUnitSystem->meanPositions();
      for ( SiteID s = 1; s <= nTotalSites; s++ ) {
        currentPos[s] = froda2_MeanPositions[s-1];
      }
      
      //update other things that Froda needs
      fitCycles = minim->getFinalStepNum();
      getRMSD();
    }
    else
    { //begin regular Froda fitting

      // reset the cyclesToFit variable 
      cyclesToFit = 0;
  
      fitCycles=0;
      int do_more = atLeastCycles;
  
      while(keepTrying){ //THIS IS THE FITTING LOOP IN WHICH WE REESTABLISH CONSTRAINTS
  
        //SAW 2 DAN: my minimisation routine happens inside this loop
        //Your routine either takes place in here
        //or it replaces this loop entirely
        //if you replace this loop, check that we take care of
        //hatbox constraint and Rama checks
  
        fitCycles++;
  
        //cerr << "Fit cycle " << fitCycles << endl;
  
        problemReport.clear();
        stringstream probLine;
        probLine << "CONF " << setw(6) << i+1;
        problemReport.push_back(probLine.str() ); 
  
  
        //SAW 2 DAN
        //bodyResponse is my previous fitting approach
        //should be obsolete with your system
        //so this part can be skipped if we're using your system
        if ((use_bodyResponse) && (fitCycles%bodyResponseEvery ==0)) {
          do_bodyResponse=true;
        }
        else{
          do_bodyResponse = false;
        }
  
        regridCounter++;
  
        if (verbose && chatty_fitting) cerr << "Cycle " << setw(2) << fitCycles << ", conformer " << setw(6) << i+1 << ": ";
        notGood = 0;
        
        // zero the geometric mismatch arrays
        clearMismatchArrays();
        //SAW 2 DAN: your system doesn't use mismatches?
  
        //SAW 2 DAN: changeweights also not applicable to your system
        for ( unsigned int siteNumber = 1; siteNumber <= structure.total_sites; siteNumber++ ) {
          frodaAtom.at(siteNumber).weightIncreasing = false;
        }
        if ( parameters.changeWeights && fitCycles > 10 ) {
          checkingWeights = true;
        }
  
        //SAW 2 DAN: these next few calls are my minimiser
        //your system is an alternative to these calls
        // Update the ghost bodies to the new atomic positions
        fitGhosts();
       
        // Give each atom its mismatches to the ghost bodies
        notGood += findGhostMismatch();
        
        // Apply hydrophobic tethers
        notGood += phobicMismatch();
        
        sterics->update();
        // Apply vdw mismatches
        if (usesterics) notGood += stericMismatch( structure );
   
        //apply the hatbox constraint, if needed
        if ( parameters.useHatbox ) {
          applyHatbox();
        } 
        //SAW 2 DAN: do you have the hatbox constraint implemented?
        //it keeps a protein restricted to a cylindrical volume of space
        
  
        //SAW 2 DAN: by now my routine has finished one interation of fitting
   
        // check ramachandran constraints
        if (checkRama) checkRamaPlot( structure);
   
        getRMSD();
  
        if (verbose && chatty_fitting){
          cerr  << "  ALLATOM RMS-D = " << globalRMSD << endl;
        }
        
        if (RMSD != RMSD){ // not a number!
          // PUT IN FATAL ERROR WARNING
          cerr << "Numerical error (RMSD NAN) causing death immediately for conformer " << i+1 << endl << endl;
          outputConformer(structure, i);
          exit(1);
        }        
      
        //check weighting behaviour
        if ( parameters.changeWeights ) {
          checkWeights();
        }
  
        //cerr << "Updating atoms" << endl;
        updateAtoms(structure);
        //cerr << "Done update atoms" << endl;
  
        if(!notGood){
          do_more -= 1; //count down; structure must be good for (atLeastCycles)
        } else{
          if (do_more < atLeastCycles) do_more +=1;
          if(isBad) keepTrying=false;
        }
  
  
        if (do_more ==0) {
          keepTrying =false;
        }
   
        if (fitCycles > maxFitCycles){ //FAILURE TO CONVERGE- RECOVER
  
          //SAW 2 DAN: I assume your system has its own alternative to this?
          //and ideally your system never gets to this situation :)
  
          if (verbose){
            cerr << "FAILED AFTER MAX CYCLES: " << maxFitCycles << ", RECOVERING." << endl;
          }
  
          for ( unsigned int line = 0; line < problemReport.size(); line++) {
            problemLog << problemReport.at(line) << endl;
          } 
  
          currentPos = cachedPos;
  
          clearMismatchArrays();
          fitGhosts();
      
          sterics->update();
  
          if ( usesterics && stericMismatch( structure ) > 0 ) {
            isLethal = true;
          }
  
          clearMismatchArrays();
          updateAtoms( structure );
  
          nFailedCycles++;
          failsInRow++;
          //cerr << "Fails in row: " << failsInRow << endl;
          keepTrying = false; //leave with this recovery position as is
  
          getRMSD();
          didIFail = true; // was a problem
        } //end of recovery from jam
      } // found a good configuration

    
      //SOME POST-PROCESSING ON NEW CONFORMER
  
      //SAW 2 DAN: hopefully we don't need this lethality check if your system is working
      if ( !isLethal ) {
        failsInRow = 0;
      }
      isLethal = false;
      if ( failsInRow > maxFailsInRow ) {
        isLethal = true;
        //cerr << "Fails in row: " << failsInRow << endl;
      }
  
      if ( didIFail ) {
           failScore +=2;
           didIFail = false; // reset
      }
      else{
           if ( failScore > 0 ) failScore-- ; // decrement if we're running OK.
      }
  
      if (isBad){ // not a number!
        //FATAL ERROR
        cerr << "Numerical error causing death immediately for conformer " << i << endl << endl;
        cerr << "Printing wrong conformer and shutting down." << endl;
        outputConformer(structure, i);
        exit(1);
      } 
  
      if ( isLethal ) {
        cerr << "WARNING: Failed too many times on structure containing bad sterics." << endl;
        cerr << "Printing wrong conformer and shutting down." << endl;
        outputConformer(structure, i);
        exit(1);
      }
      
      if ( failScore > failLimit ){
           cerr <<"WARNING: too many successive failure in this run." << endl;
           cerr <<"Printing wrong conformer and shutting down." << endl;
           outputConformer(structure, i);
           exit(1);           
      }
    } //end of regular Froda Fitting

    if (isRestart || isTargeted ){ 
      distanceToTarget(restart, target);
    }
    
    getRMSDByResidue(); //get the current and running mobilities by residue;
    nConfsFound++; //count for resrmsd array
    running_allmain_RMSD = sqrt( running_allmain_MSD / nConfsFound );


    if ( parameters.stopAtEffectiveTimeLimit && !froda2Hybrid ) {
      pathLengthIntegrator->updatePath();
    }
    
    // BMH commented out. Need better memory management
    if ( doGetArea ) {
      getAreaOfConformer();
    }

    Rg = calculateRg();
    if ( doGetPhobicRg ) {
      phobicRg = calculatePhobicRg();
    }
    if (usePairs) checkPairs();

    double correlation = 0;
    if (useED && doAnneal) {
      correlation = checkMonteCarloED(*myPdb,*targetMap);
    }
    if (useSAXS && doAnneal) {
      correlation = checkMonteCarloSAXS(*mySAXS);
    }

    //check if a targeted run has finished
    if ( foundTarget ) {
      if (verbose) cout << " STRUCTURE CONVERGED TO TARGET: OUTPUTTING AND SHUTTING DOWN." << endl;

      if (reportByRSMD){
        frameNumber++;
        outputConformer (structure, frameNumber -1);
      } else{
        outputConformer(structure,i);
      }
      whenTargetFound = i;
      break; 
    }

    if (verbose && chatty_fitting){
      if (isTargeted and !foundTarget){
        cerr << "TARGET NOT YET FOUND." << endl;
      }
    }

    
    if (isTargeted && doAnneal){ 
      checkMonteCarloAnneal( structure);
    }
    
    double currentAccept = double ((i+1 - mcRejected)*1.0/(i+1));
    if (doAnneal && constantRate) { // adjust annealing scale to get correct accept rate
      if ((currentAccept > acceptRate) && (annealScale > 1e-10)) { // accept rate too high
        annealScale = annealScale * 0.5;
      }
      else if ((currentAccept < acceptRate) && (annealScale < 1e10)){ // accept rate too low
        annealScale = annealScale * 2.0;
      }
    }  
    if ( doScoringMC ){
        clearMCFoldingVariables();
        getMCScore( structure );
        checkMonteCarloBonding( structure);
    }

    if (parameters.flexweb) {
      *runningRMSDfile << showpoint << fixed << setprecision(3)  << setw(11) << (i+1) << " \t " 
            << setw(7) << allmain_RMSD << " \t " 
            << setw(10) << running_allmain_RMSD << endl;
    }
      
    if (verbose){
      cerr << "CONF. "  << setw(6) << i+1 << " FOUND at " << setw(2) << fitCycles << " cycles,";
      cerr << " ALLATOM "  << setw(8) << globalRMSD << "; ";

      if( isRestart){ cerr << "RESTART " << setw(8) << restartRMSD << "; "; }
      if (isTargeted){ cerr << "TARGET " << setw(8) << targetRMSD << "; ";}
      if (useED || useSAXS) { cerr << "CORR " << setw(8) << correlation << "; "; }
      if (doAnneal) { 
        cerr << "MCACCEPT " << setw(4) << currentAccept*100.0 << "%; "; 
      }
      if (doAnneal && constantRate) {
        cerr << "ANNEAL " << setw(4) << annealScale << "; ";
      }
      if ( parameters.stopAtEffectiveTimeLimit ) {
        cerr << "EFF TIME PS " << pathLengthIntegrator->getIntegratedPathLength() / s_over_t_Carbon_Angstroms_per_ps << "; "; 
      }
      cerr << endl;
    }
    //decide whether to report the conformer to file
    //SAW 2 DAN: my output routine involves reporting some ghost mismatches
    //and steric overlaps: the same information should be available if your system ran
    if (isTargeted && reportByRSMD){
      thisRMSD = targetRMSD;
      if ( (lastRMSD - thisRMSD) > RMSDSpacing ){
        frameNumber++;
        outputConformer(structure, frameNumber -1); outputCounter=0;
        lastRMSD = thisRMSD;
      }
    } else if (reportByRSMD){
      thisRMSD = globalRMSD;
      if ( (thisRMSD - lastRMSD) > RMSDSpacing ){
        frameNumber++;
        outputConformer(structure, frameNumber -1); outputCounter=0;
        lastRMSD = thisRMSD;
      }
    } else if (!parameters.stopAtEffectiveTimeLimit && outputCounter >= outputEvery){
      outputConformer(structure, i); outputCounter=0;
      if (useSAXS) {
        stringstream myfilename;
        myfilename << structure.base_name <<"_froda_" << setfill('0') << setw(8) << i+1 << ".saxs";
        mySAXS->writeSAXS(myfilename.str());
      }
    } else if ( parameters.stopAtEffectiveTimeLimit && 
		pathLengthIntegrator->getIntegratedPathLength() / s_over_t_Carbon_Angstroms_per_ps >= parameters.effectiveTimeFreq*(frameNumber+1) ) {
      frameNumber++;
      outputConformer(structure, frameNumber-1); 
    }
    if ( runMorph && ( frameNumber == morphFrames) ){
      cerr << "Final frame generated, shutting down." << endl;
      exit(0);
    }


  } //END OF PRODUCTION LOOP

  // Report any final messages and shut down
  
  //if (verbose) cerr << "Goodbye and thank you." << endl;
  timer.stop();
  
  if( parameters.verbose ){
    cout << "Total FRODA run time: " << showpoint << setiosflags(ios::fixed) << setprecision(3) << static_cast<double>(timer.getutime())/1000.0 << " seconds." << endl;
    cout << "Total Failed cycles: " << nFailedCycles << endl;

    if (doScoringMC){
      cout << "Folding Monte-Carlo statistics:" << endl;
      cout << "   Conformers accepted: " << n_accepted;
      cout << "   Conformers rejected: " << n_rejected;
      cout << "   accept:reject ratio = " << ((double) n_accepted) / ( (double) n_rejected);
      cout << endl;
    }
  
    if( whenTargetFound ){
      cout << "Target position reached after " << whenTargetFound << " cycles." << endl;
    }
  } 

  delete sterics;
 // delete electron-density specific pointers
  delete myPdb;
  delete targetMap;
  delete areaProximity;

  parameters.mapFromTaskNameToStatus["Generate conformations"] = "Complete";
  outputStatus();

  if (runningRMSDfile != NULL) {
    runningRMSDfile->close();
  }

  problemLog.close();

  return;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description: The kill handler - called whenever an outside program attempts 
// to kill FRODA before it has exited 
////////////////////////////////////////////////////////////////////////////////
void Froda::killFRODA(int signal) {
  parameters.mapFromTaskNameToStatus["Generate conformations"] = "Interrupted";
  outputStatus();

  exit(EXIT_SUCCESS);
}



////////////////////////////////////////////////////////////////////////////////
// Description: update status report to parameters maps
////////////////////////////////////////////////////////////////////////////////
void Froda::updateStatus() {

  stringstream global_RMSD;
  stringstream main_RMSD;
  stringstream total_confs;
  
  global_RMSD << "" << showpoint << fixed << setprecision(3) << globalRMSD << " &Aring; " << endl;
  
  main_RMSD << "" << showpoint << fixed << setprecision(3) << allmain_RMSD << " &Aring; " << endl;

  total_confs << i+1 << " / "<< nConfigs << endl;
  
  parameters.mapFromTaskNameToStatus["Current All-Atom RMSD"] = global_RMSD.str();
  parameters.mapFromTaskNameToStatus["Current Main Chain RMSD"] = main_RMSD.str();
  parameters.mapFromTaskNameToStatus["Conformations generated"] = total_confs.str();

  outputStatus();
}

