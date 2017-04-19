#ifndef _CLASS_FRODA_
#define _CLASS_FRODA_

#include "global_defs.h"
#include "MolFramework.h"
#include "interact.h"
#include "mt19937ar.h"
#include "flexweb.h"
#include <signal.h>

#include "SiteSelector.h"
#include "AllSiteSelector.h"
#include "CompositeSiteSelector.h"
#include "BackboneSiteSelector.h"
#include "TargetedSiteSelector.h"
#include "SidechainSiteSelector.h"
#include "SiteID.h"
#include "EDStructure.h"
#include "SAXS.h"
#include "GenericMap.h"
#include "Timer.h"
#include "SymReplica.h"
#include "SymmetryMatrices.h"

class Sterics;
class AreaProximity;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
struct Ghost{
  public:
  vector< unsigned int > atoms;
  vector< Vector > bondVector, pos, posAtom;
  Vector posCentral;
  double bodyweight, radMS, radRMS;
  bool isSidechain;
  Vector runningRotor, baseRotor;
  unsigned int bodysize;

  Vector oldPosCentral, olderPosCentral;
  Vector oldRotor, olderRotor;
  //int cyclesToFit;

  bool changedWeight;

  Ghost(){
    atoms.clear();
    bondVector.clear();
    pos.clear();
    posAtom.clear();
    posCentral = Vector(0,0,0);
    bodyweight =0;
    radMS = 0;
    radRMS = 0;
    isSidechain = false;
    runningRotor = Vector(0,0,0);
    baseRotor = Vector(0,0,0);
    bodysize = 0;

    oldPosCentral = Vector(0,0,0);
    olderPosCentral = Vector(0,0,0);
    oldRotor = Vector(0,0,0);
    olderRotor = Vector(0,0,0);
    //cyclesToFit = 0;

    changedWeight = false;
  }

};

struct myPhobe{
 public:

 unsigned int atomA, atomB; double AB;
 myPhobe(unsigned int A_in=0, unsigned int B_in=0, double AB_in=0){
  atomA=A_in; atomB=B_in; AB=AB_in;
 }
}; 

struct myPair{
 public:

 unsigned int atomA, atomB; double AB;
 double ABcurrent;
 bool tagA;
 bool tagB;
 bool tagGood;
 myPair(unsigned int A_in=0, unsigned int B_in=0, double AB_in=0,
        double ABcurrent_in = 0, bool tagA_in = true, bool tagB_in = true, bool tagGood_in = false){
  atomA=A_in; atomB=B_in; AB=AB_in;
  ABcurrent = ABcurrent_in;
  tagA = tagA_in;
  tagB = tagB_in;
  tagGood = tagGood_in;
 }
}; 

struct FrodaResidue{
  public:

  unsigned int cAlpha;
  vector< unsigned int > sameResidue; // moved from frodaAtom construct
  unsigned int n1;
  unsigned int n2;
  Vector e1, e2, e3;
  bool isBiased;
  Vector originalBias, internalBias, externalBias;
  FrodaResidue(){
    cAlpha = 0;
    n1 = 0;
    n2 = 0;
    e1 = Vector(0,0,0);
    e2 = Vector(0,0,0);
    e3 = Vector(0,0,0);
    isBiased = false;
    originalBias = Vector(0,0,0);
    internalBias = Vector(0,0,0);
    externalBias = Vector(0,0,0);
  }

  void makeInternalBias() {
    internalBias.x = dotProduct( e1, originalBias );
    internalBias.y = dotProduct( e2, originalBias );
    internalBias.z = dotProduct( e3, originalBias );
  }
  void makeExternalBias() {
    externalBias.x = internalBias.x*e1.x +
                     internalBias.y*e2.x +
                     internalBias.z*e3.x;
    externalBias.y = internalBias.x*e1.y +
                     internalBias.y*e2.y +
                     internalBias.z*e3.y;
    externalBias.z = internalBias.x*e1.z +
                     internalBias.y*e2.z +
                     internalBias.z*e3.z;
  }

};

////////////////////////////////////////////////////////////////////////////////
struct myRama{
 public:

 int Calpha, Ni, Ci, Ci_1, Ni1;
 double phi0, psi0;
 bool freePhi, freePsi, badStart;

 int mapValue;

 bool isPro;
 bool isGly;

 myRama(int Calpha_in=0, int Ni_in=0, int Ci_in=0, int Ci_1_in=0, int Ni1_in=0,
  double phi0_in=0.0, double psi0_in=0.0, bool freePhi_in=1, bool freePsi_in=1, bool badStart_in=0,
  int mapValueIn = 0, bool isProIn = false, bool isGlyIn = false
  ){ Calpha=Calpha_in; Ni=Ni_in; Ci=Ci_in; Ci_1 = Ci_1_in; Ni1=Ni1_in;
  phi0=phi0_in;psi0=psi0_in; freePhi=freePhi_in; freePsi=freePsi_in; badStart=badStart_in;
  mapValue = mapValueIn; isPro = isProIn; isGly = isGlyIn;
 }
};

////////////////////////////////////////////////////////////////////////////////
struct FrodaAtom{
  public:

  bool isCore;
  bool isMobile;
  bool isMain;
  bool isPolar;
  bool isHydrogen;
  bool isTerminal; //has only one covalent neighbor
  bool isPhobic;
  bool isAcceptor;
  bool isAlpha; //labels c-alpha atoms

  bool hasWeight;
  bool weightIncreasing;
  bool weightDecreasing;

  bool moved;

  double weight;
  double oldWeight;
  double offset;

  int charge;

  unsigned int myTarget;

  Vector mainMismatch;
  int nMainMismatches;

  Vector randomMove;

  vector< unsigned int > friends;
  vector< unsigned int > ghostlist;
  vector< unsigned int > neighbors;

  vector< double > surfaceCk;
  double surfaceRho;
  double SASA;
  double polarPairSA;

  FrodaAtom() {
    isCore = false;
    isMobile = true;
    isMain = false;
    isPolar = false;
    isHydrogen = false;
    isTerminal = false;
    isPhobic = false;
    isAcceptor = false;
    isAlpha = false;
    hasWeight = false;
    weightIncreasing = false;
    weightDecreasing = false;
    weight = 1.0;
    oldWeight = 1.0;
    offset = 0;
    charge = 0;

    myTarget = 0;

    randomMove = Vector( 0,0,0 );
    
    friends.clear();
    ghostlist.clear();
    surfaceCk.clear();
    surfaceRho = 0;
    SASA = 0;
    polarPairSA = 0;
  }
 
};

////////////////////////////////////////////////////////////////////////////////
class Froda {
 public:

  Timer timer;

  vector< string > problemReport;

  int nConfigs;
  int notGood;
  int whenTargetFound;
  int i; // configuration counter
  bool keepTrying;
  bool debug_vdw, debug_rigid, debug_connect, debug_phobe, debug_hbond, debug_ghost;
  bool chatty_fitting; // report every cycle
  bool isBad;
  bool isRestart, isTargeted, foundTarget; // for targeted dynamics routines
  double tolTarget; // get this close to the target 
  double targetStep; // move atoms this far during targeted dynamics;
  double moveProb;
  
  bool doAnneal; // for simulated annealing routine
  double annealScale; //scaling on distance changes
  bool constantRate; // maintain constant MC accept rate?
  double acceptRate; // constant MC accept rate to maintain
  int mcRejected;   // number of rejected Monte Carlo moves
  vector< Vector > cachedPos; // store positions for recovery
  double cachedRMSD;
  
  bool cylindrical; // for cylindrical centrifuge for membrane proteins
  double lowZRange, highZRange; // range of z coordinates in which to apply centrifuge
  
  bool localTargeting; // target locally
  double lstep; // local targetting step size, differs from dstep
  
  bool propTotarget;//scales target steps according to distance from target
  
  // tolerances during fitting
  double ghostTol; // mismatch to ghost positions, Angstroms
  double phobeTol; // tolerance on hydrophobic variation
  double vdwTol; // acceptable vdw overlap
  double vdwNotice; // enough contact to matter
  double epsilon; // very small vector length, equivalent to zero
  double bounce; // random step size, Angstroms
  int cyclesToFit;

  double worstGhostMismatch;
  unsigned int worstGhostAtom;
  unsigned int worstGhost;

  double worstClash;
  unsigned int worstClashAtom1;
  unsigned int worstClashAtom2;

  double worstTetherStretch;
  unsigned int worstStretchAtom1;
  unsigned int worstStretchAtom2;
  
  double verysmall;
  
  int regridEvery, regridCounter; // update vdW coarse gridding.
  
  int failsInRow, maxFailsInRow; // stop when jammed
  bool didIFail;
  int failScore, failLimit;
  
  bool centrifugal; // if yes, bias the random steps to unfold the protein
  double centri_maxr; // only apply centrifuge if this close to center
  bool localCentrifuge; //use local centers not geometric center
  vector< unsigned int > centerID; //use to identify "local center" atoms
  vector< int > centerDirection; // towards or away?
  
  int outputEvery; // report every nth configuration
  int outputCounter;
  
  double thisRMSD, lastRMSD;
  int frameNumber;
  
  int cacheEvery;
  int fitCycles;
  int maxFitCycles;
  int atLeastCycles;
  unsigned int maxRotorCycles;
  int nMobileAtoms, nMainchainAtoms, nSidechainAtoms, nMobileMainchainAtoms, nMobileSidechainAtoms;
  int nTargetedMainchainAtoms;
  
  unsigned int nTotalSites; //same as structure.total_sites;
  unsigned int nTotalClusters; //same as structure.total_clusters
  
  double MSD, RMSD; // rmsd on mobile atoms
  double globalMSD, globalRMSD; //rmsd on all atoms
  double restartMSD, restartRMSD; // rmsd to restart point
  double targetMSD, targetRMSD; // rmsd to target point
  
  bool fancy_targetting; // check original pdb numbers for partial target files
  
  double allmain_MSD, allmain_RMSD; // rmsd on all mainchain atoms
  double allTargetedMainchainMSD, allTargetedMainchainRMSD; //targeted mainchain atoms  
  double running_allmain_MSD, running_allmain_RMSD; //cumulative RMSDs for flexweb graphs
  double targmain_RMSD; // mainchain target
  double allside_MSD, allside_RMSD; // rmsd on all sidechain atoms
  double mobmain_MSD, mobmain_RMSD; // rmsd on mobile mainchain atoms
  double mobside_MSD, mobside_RMSD; // rmsd on mobile mainchain atoms
  vector< bool > hasBadClash;
  
  vector< double > currentResidueRMSD, runningResidueRMSD;
  int nConfsFound;
  
  Vector pq, pqdash,b, bond1, bond2; // my Vector object, cartesian
  Vector rotor, tempvec, dummyvec, myvec, gradient;
  static const Vector NULL_VEC;
  double veclength, radius, temp, dot, myX, modgrad,smallgrad, large_rotor_2;

  vector< Ghost > ghost;
  vector< FrodaAtom > frodaAtom;
  
  double range;
  
  int nAtomsWithTargets;

  int bodyResponseEvery; // do the response thing this often 
  bool do_bodyResponse;
  bool use_bodyResponse;

  bool checkingWeights;
  
  int n_myHbonds, n_myPhobes;
  
  vector<myPhobe> Phobes;
  myPhobe dummyPhobe;
  vector<int> chargedAtoms, phobicAtoms;
  
  double myXLength, myYLength, myZLength;
  int myXCells, myYCells, myZCells;
  Vector myCentralPoint;
  
  //these variable are used to count VDW overlaps and sort them by severity;
  //e.g. an overlap goes into the bin (overlap/overlapStep)
  vector< int > overlapDegrees;
  double overlapStep;
  int nOverlaps;
  
  //per-residue list: get all the calphas and the phobic CBetas
  vector< unsigned int > alphaCarbon;
  vector< unsigned int > phobicBetaCarbon;

  bool moveByCluster; // bias collective motions of large rigid bodies
  bool makeSingleMove; // move only one atom or body at a time
  bool moveDihedral; //use angular perturbations

  bool useInternalBias;
  bool useElasticVector;
  vector< unsigned int > nodalAtoms; //.resize( nNodes );
  vector<  Vector  > eigenvector; //
  vector< FrodaResidue > nodeResidues; //.resize to suit alphas

  
  int nFailedCycles;
  
  myRama dummyRama;
  vector< myRama> ramaList;
  bool checkRama;
  static const double RAD_TO_DEG;
  vector<vector< int > > generalRamaPlot, prolineRamaPlot, glycineRamaPlot;
  int nBadRama, nGenerousRama, nAllowedRama, nCoreRama;
  
  unsigned int myNClusters; // structure.total_clusters; add 1 if it's mobile rc1
  
  static const double TWOPI;
  
  bool   useCM ;
  bool   useCMforward ;
  bool   useCMreverse ;
  bool   cmdebug;
  double cmforwardfactor ;
  double cmreversefactor ;
  int    cmlowerFitcycles ;
  int    cmupperFitcycles ;
  double cmmaxstep ;
  int CMWorstFitCycles;

  bool persist;
  int nPersist;
  int persistCounter;
  
  bool useSeed;
  long int seed;
  bool mobileRC1;
  bool verbose;
  bool interactive;
  bool runMorph;
  int morphFrames;
  double RMSDSpacing;
  bool reportByRSMD;
  bool polarGeometry;
  double polarHRadius;
  bool use_group_id;
  double energy_cutoff;
  
  bool doScoringMC; //should we run the folding monte carlo
  double frodaMCScale; //normaliser for delta-e in MC accept-reject

  double Rg;
  double oldRg;

  double phobicRg;
  double oldPhobicRg;
  bool doGetPhobicRg;

  double RgWeight;
  bool mcRg;
  double minimumRg;

  bool mcPhobicRg;
  double phobicRgWeight;
  double minimumPhobicRg;

  bool mcSASA;
  double SASAWeight;
  double minimumSASA;

  bool mcPhobicSASA;
  double phobicSASAWeight;
  double minimumPhobicSASA;
  
  bool mcPolarSASA;
  double polarSASAWeight;

  bool doHBScoring;
  double HBWellDepth;
  double SBWellDepth;

  bool doGetArea;
  double totalSASA;
  double phobicSASA;
  double polarSASA;
  double polarPairSA;
  double oldTotalSASA;
  double oldPhobicSASA;
  double oldPolarSASA;
  double oldPolarPairSA;

  double frodaMCScore, oldFrodaMCScore;
  double RgScore, oldRgScore;
  double phobicRgScore, oldPhobicRgScore;
  double SASAScore, oldSASAScore;
  double phobicSASAScore, oldPhobicSASAScore;
  double polarSASAScore, oldPolarSASAScore;
  double HBScore, oldHBScore;

  vector< int > donorHydrogenList; //hbonding_atoms_list; //scan only over these atoms for finding hbonds
  
  int n_hb_3_10, n_hb_alpha; //count the number of i -> i+3 and i -> i+4 backbone hbonds for helix tracking
  int old_n_hb_3_10, old_n_hb_alpha;
  int n_accepted, n_rejected;
  
  vector< Vector> initialPos, currentPos;

  bool nosterics;
  //bool useTalksto;

// for electron-density options
  bool useED;
  bool useTheoMap;
  bool useEZDMap;
  bool trimMap;  // use TrimMap class
  double edResFac;  // electron density resolution factor
  double edNoise; // noise level in theoretical ED map
  string theoMapFile; // PDB file from which to generate an ED map
  string ezdMapFile;
  double tolED;   // get this close to ED map when targeting
  int trimMapFactor;
// for SAXS profiles
  bool useSAXS;
  double saxsTol;
  string saxsFile;

  vector< unsigned int > possibleContacts;
  vector< unsigned int > possiblePhobic;
  vector< unsigned int > possiblePolar;

  double initialToTargetL;

  bool usePairs;
  bool stopWhenAllPairConstraintsSatisfied;
  vector< myPair > myPairs;
  int nMyPairs;

  Sterics *sterics;
  AreaProximity *areaProximity;
  
  bool froda2Hybrid;

  SymmetryMatrices *symMat;

  Froda();
  ~Froda();
  void runFRODA( MolFramework &structure, MolFramework &target, MolFramework &restart );
  static void killFRODA(int signal);
  void newPerturbation( MolFramework &structure);
  void newRandomMove();
  void oldRandomMove();
  void fitGhosts();
  void saveGhostMomentum();
  int findGhostMismatch();
  int phobicMismatch();
  int stericMismatch( MolFramework &structure );
  void updateAtoms( MolFramework &structure );
  void updateStatus();
  void outputConformer( MolFramework &structure, int i );
  void outputGhosts(MolFramework &structure );
  //void initialiseGhosts( MolFramework &structure );
  void stripNonCovalentNeighbors( MolFramework &structure );
  void getGhostMembershipFromFile( MolFramework &structure );
  void getRCDFromGhosts( MolFramework &structure );
  void getGhostMembershipFromRCD( MolFramework &structure );
  void buildGhostGeometry( MolFramework &structure );
  void initialisePhobic( MolFramework &structure );
  void findPhobicBetaCarbons( MolFramework &structure );
  void centrifuge(int an_atom);
  void cylinderCentrifuge(int an_atom);
  int checkRamaPlot(MolFramework &structure );
  void makeRamaList( MolFramework &structure);
  void makeBackboneList( MolFramework &structure );
  void getGeneralRamaPlot();
  void getGlycineRamaPlot();
  void getProlineRamaPlot();
  void moveToTarget( MolFramework &aim);
  void distanceToTarget( MolFramework &restart, MolFramework &target);
  void setUpLocalTargeting ( MolFramework &structure, MolFramework &target );
  void moveToLocalTarget ();
  void getRMSDByResidue();
  
  void initialiseMainArrays( MolFramework &structure, MolFramework &restart );
  void defineAtomsAsMainOrSide( MolFramework &structure);
  void defineClustersAsMainOrSide();
  //void setUpTalkstoList( MolFramework &structure );
  void setUpRMSDByResidueArrays();
  void identifyPolarAtoms(MolFramework &structure);
  void obtainSurfaceParameters(MolFramework &structure);
  void identifyPhobicAtoms(MolFramework &structure);
  void setUpTargetList(MolFramework &structure, MolFramework &target);
  void applyCentrifuge();
  void clearMismatchArrays();
  void getRMSD();
  void checkMonteCarloAnneal( MolFramework &structure);
  double checkMonteCarloED(EDStructure &currentStructure, GenericMap &edMap);
  double checkMonteCarloSAXS(SAXS &saxs);
 
  void getAreaOfAtom( unsigned int atomID );
  void getAreaOfConformer();
  double returnAreaCalculation( double c0, double c1, double c2, double c3,
                                double c4, double Ai );
 
  void findDonorHydrogens( MolFramework &structure );
  void findHBondingAtoms( MolFramework &structure );
  void findChargedInteractions( MolFramework &structure );
  void checkMonteCarloBonding( MolFramework &structure);
  void storeMCFoldingVariables();
  void restoreMCFoldingVariables();
  void clearMCFoldingVariables();
  void getMCScore( MolFramework &structure);


  bool allPairConstraintsSatisfied();
  bool isThirdNeighbor( unsigned int atom1, unsigned int atom2 );
  bool hasMoved( unsigned int atom, MolFramework &structure );
  void findNeighbors( MolFramework &structure);
  bool bonded( vector<unsigned int> &vec1, vector<unsigned int> &vec2 );
  bool hasHBondGeometry( unsigned int atom1, unsigned int atom2, vector<unsigned int > &nlist1, vector<unsigned int> &nlist2 );

  Vector getNormalVector( unsigned int atom1, bool &isFlat, MolFramework &structure );
  double computeAngle( SiteID atom1, SiteID atom2, SiteID atom3 );
  double getDist2( SiteID atom1, SiteID atom2 );

  bool isSaltBridge( SiteID hydrogen, SiteID acceptor, MolFramework &structure );
  bool maybeHBond( SiteID hydrogen, SiteID acceptor, MolFramework &structure );
  int HBType( SiteID hydrogen, SiteID acceptor, MolFramework &structure );
  double SBEnergy( SiteID donor, SiteID acceptor );
  void scanPhi( SiteID hydrogen, SiteID acceptor, MolFramework &structure, 
                     int hbType, double &phi, bool &isOK );
  double scoreHBond( SiteID hydrogen, SiteID acceptor, MolFramework &structure );

  void applyMomentum( int ghostid );
  void reverseMomentum( int ghostid );

  void getStericParms( const MolFramework& structure, SiteID atom1, SiteID atom2, double& rInner, double& rOptimum, double& rOuter );

  bool isBlocking( unsigned int atom1, unsigned int neighbor, unsigned int atom2, double margin );

  void applyHatbox();

  double calculateRg(); //gets the rg of the protein
  double calculatePhobicRg(); //gets the phobic rg of the protein

  void applyPairs();
  void readPairFile( MolFramework &structure);
  void checkPairs();

  void readElasticNetworkMode( MolFramework &structure );
  void makeResidueList( MolFramework &structure);
  void makeResidueBasis();
  void formBiasVectors();
  void addModeBias();
  void applyCMbias();

  void pickPeptideDihedral( SiteID &CA, SiteID &base, MolFramework &structure );
  void perturbBackboneDihedral( SiteID CA, SiteID base, MolFramework &structure );
  void pickPeptidePlane( SiteID &CA1, SiteID &b1, SiteID &b2, SiteID &CA2, MolFramework &structure );
  void perturbPeptidePlane( SiteID CA1, SiteID b1, SiteID b2, SiteID CA2, MolFramework &structure );
  void twitchSidechain( MolFramework &structure );

  void setupTheoEDObjects( MolFramework &structure, EDStructure* &myPdb, GenericMap* &targetMap );
  void setupRealEDObjects( MolFramework &structure, EDStructure* &myPdb, GenericMap* &targetMap );
  void setupEDTrimMap( MolFramework &structure, EDStructure* &myPdb, GenericMap* &targetMap );

  void resetWeights( );
  void checkWeights( );

  Vector estimateRotor(int ghostid);
  
  void obtainParameters( MolFramework &structure);
  double Xval(Vector &b);
  double epsx(Vector &pq, double pqdash_x, Vector &b, double myX);
  double epsy(Vector &pq, double pqdash_y, Vector &b, double myX);
  double epsz(Vector &pq, double pqdash_z, Vector &b, double myX);
  double depsxdbx(Vector &pq, Vector &b, double myX);
  double depsxdby(Vector &pq, Vector &b, double myX);
  double depsxdbz(Vector &pq, Vector &b, double myX);
  double depsydbx(Vector &pq, Vector &b, double myX);
  double depsydby(Vector &pq, Vector &b, double myX);
  double depsydbz(Vector &pq, Vector &b, double myX);
  double depszdbx(Vector &pq, Vector &b, double myX);
  double depszdby(Vector &pq, Vector &b, double myX);
  double depszdbz(Vector &pq, Vector &b, double myX);
  double deps2dbx(Vector &pq, Vector &pqdash, Vector &b);
  double deps2dby(Vector &pq, Vector &pqdash, Vector &b);
  double deps2dbz(Vector &pq, Vector &pqdash, Vector &b);
  Vector crossProduct(Vector &bond1, Vector &bond2);
  Vector randomUnitVector( );
  double dotProduct (Vector &bond1, Vector &bond2);
 
  Vector getGradient( double &weight, Vector &pq, Vector &pqdash, Vector &b );
  
  void setMatrices(vector <Matrix> *matrices_);
 
  SiteSelector *allSiteSelector;
  SiteSelector *backboneSiteSelector;
  SiteSelector *sidechainSiteSelector;
  SiteSelector *targetedSiteSelector;
  SiteSelector *targetedBackboneSiteSelector;
  
  SiteSelector *preferredSiteSelector;
  
  float currentRMSDfromInitial;
  float currentRMSDfromTarget;
};

// Here are all the geometric algebra equations for rotating bonds
// and the vector dot and cross product operations
////////////////////////////////////////////////////////////////////////////////
inline double Froda::returnAreaCalculation( double c0, double c1, double c2, double c3,
                                double c4, double Ai ) {
  double SASA;
  double runningA;
  runningA = 1;
  SASA = c0;
  runningA *= Ai;
  SASA += c1*runningA;
  runningA *= Ai;
  SASA += c2*runningA;
  runningA *= Ai;
  SASA += c3*runningA;
  runningA *= Ai;
  SASA += c4*runningA;

  return SASA;
}

//////////////////////////////////////////////////////////////////////////////////
// Decription: Scalar X value 1 - 0.25*|b|^2 for a rotor b
////////////////////////////////////////////////////////////////////////////////
inline double Froda::Xval(Vector &b){
  return sqrt( 1.0 - 0.25 * (b.x*b.x + b.z*b.z + b.y*b.y)); 
}

////////////////////////////////////////////////////////////////////////////////
// Description: x component of mismatch between a vector pqdash
//              and a vector pq rotated by a rotor b
//              only the x component of pqdash is passed
//              myX is the X scalar belonging to b
////////////////////////////////////////////////////////////////////////////////
inline double Froda::epsx(Vector &pq, double pqdash_x, Vector &b, double myX){
  return pq.x*(1.0 - 0.5*(b.y*b.y+ b.z*b.z))
    + pq.y*(-myX*b.z + 0.5*b.x*b.y)
    + pq.z*(myX*b.y + 0.5*b.x*b.z)
    - pqdash_x;
  
}

////////////////////////////////////////////////////////////////////////////////
// Description: y component of mismatch between a vector pqdash
//              and a vector pq rotated by a rotor b
//              only the y component of pqdash is passed
//              myX is the X scalar belonging to b
////////////////////////////////////////////////////////////////////////////////
inline double Froda::epsy(Vector &pq, double pqdash_y, Vector &b, double myX){
  return pq.y*(1.0 - 0.5*(b.x*b.x + b.z*b.z))
    + pq.z*(-myX*b.x + 0.5*b.y*b.z)
    + pq.x*(myX*b.z + 0.5*b.x*b.y)
    - pqdash_y;
  
}

////////////////////////////////////////////////////////////////////////////////
// Description: z component of mismatch between a vector pqdash
//              and a vector pq rotated by a rotor b
//              only the z component of pqdash is passed
//              myX is the X scalar belonging to b
////////////////////////////////////////////////////////////////////////////////
inline double Froda::epsz(Vector &pq, double pqdash_z, Vector &b, double myX){
  return pq.z*(1.0 - 0.5*(b.x*b.x + b.y*b.y))
    + pq.x*(-myX*b.y + 0.5*b.x*b.z)
    + pq.y*(myX*b.x + 0.5*b.y*b.z)
    - pqdash_z;
  
}

// Now here are a series of terms for differentials
// deps{m}db{n}(pq,b) is the differential of eps{m} by b{n}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depsxdbx(Vector &pq, Vector &b, double myX){
  return pq.y*( (0.25*b.x*b.z/myX) + 0.5*b.y) 
    + pq.z*( (-0.25*b.x*b.y/myX) + 0.5*b.z);
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depsxdby(Vector &pq, Vector &b, double myX){
  return pq.x*(-b.y)
    + pq.y*( (0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.z*(myX - (0.25*b.y*b.y/myX));
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depsxdbz(Vector &pq, Vector &b, double myX){
  return pq.x*(-b.z)
    + pq.z*( (-0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.y*(-myX + (0.25*b.z*b.z/myX));
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depsydbx(Vector &pq, Vector &b, double myX){
  return pq.x*((-0.25*b.x*b.z/myX) + 0.5*b.y)
    + pq.y *(-b.x)
    + pq.z*(-myX + (0.25*b.x*b.x/myX)); 
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depsydby(Vector &pq, Vector &b, double myX){
  return pq.x*((-0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.z*((0.25*b.x*b.y/myX) + 0.5*b.z);
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depsydbz(Vector &pq, Vector &b, double myX){
  return pq.x*(myX - (0.25*b.z*b.z/myX))
    + pq.y*(-b.z)
    + pq.z*((0.25*b.x*b.z/myX) + 0.5*b.y);
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depszdbx(Vector &pq, Vector &b, double myX){
  return pq.x*((0.25*b.x*b.y/myX) + 0.5*b.z)
    + pq.y*(myX - (0.25*b.x*b.x/myX))
    + pq.z*(-b.x);
} 

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depszdby(Vector &pq, Vector &b, double myX){
  return pq.x*(-myX + (0.25*b.y*b.y/myX))
    + pq.y*((-0.25*b.x*b.y/myX) + 0.5*b.z)
    + pq.z*(-b.y);
} 

////////////////////////////////////////////////////////////////////////////////
inline double Froda::depszdbz(Vector &pq, Vector &b, double myX){
  return pq.x*((0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.y*((-0.25*b.x*b.z/myX) + 0.5*b.y); 
} 

// The functions deps2dbn are the mismatch 
// eps squared differentiated by rotor component bn 
// These are used to minimise eps squared wrt the rotor b

////////////////////////////////////////////////////////////////////////////////
inline double Froda::deps2dbx(Vector &pq, Vector &pqdash, Vector &b){
  double myX = Xval(b);
  
  return 2*epsx(pq, pqdash.x, b, myX)*depsxdbx(pq, b, myX)
    + 2*epsy(pq, pqdash.y, b, myX)*depsydbx(pq, b, myX)
    + 2*epsz(pq, pqdash.z, b, myX)*depszdbx(pq, b, myX);  
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::deps2dby(Vector &pq, Vector &pqdash, Vector &b){
  double myX = Xval(b);
  
  return 2*epsx(pq, pqdash.x, b, myX)*depsxdby(pq, b, myX)
    + 2*epsy(pq, pqdash.y, b, myX)*depsydby(pq, b, myX)
    + 2*epsz(pq, pqdash.z, b, myX)*depszdby(pq, b, myX);
}

////////////////////////////////////////////////////////////////////////////////
inline double Froda::deps2dbz(Vector &pq, Vector &pqdash, Vector &b){
  double myX = Xval(b);
  
  return 2*epsx(pq, pqdash.x, b, myX)*depsxdbz(pq, b, myX)
      + 2*epsy(pq, pqdash.y, b, myX)*depsydbz(pq, b, myX)
      + 2*epsz(pq, pqdash.z, b, myX)*depszdbz(pq, b, myX);
}


////////////////////////////////////////////////////////////////////////////////
// Description: scalar (dot) product of two vectors bond1 and bond2
////////////////////////////////////////////////////////////////////////////////
inline double Froda::dotProduct (Vector &bond1, Vector &bond2) {
  return bond1.x*bond2.x + bond1.y*bond2.y+bond1.z*bond2.z;
}

////////////////////////////////////////////////////////////////////////////////
// Description: vector (cross) product of two vectors bond1 and bond2
////////////////////////////////////////////////////////////////////////////////
inline Vector Froda::crossProduct(Vector &bond1, Vector &bond2){
  Vector res(0,0,0);
  res.x = bond1.y*bond2.z - bond1.z*bond2.y;
  res.y = bond1.z*bond2.x - bond1.x*bond2.z;
  res.z = bond1.x*bond2.y - bond1.y*bond2.x;
  
  return res;
} 

////////////////////////////////////////////////////////////////////////////////
// Description: tell whether two atoms are 3rd neighbors or not
////////////////////////////////////////////////////////////////////////////////
inline bool Froda::isThirdNeighbor(unsigned int atom1, unsigned int atom2){

  if ( atom1 == atom2 ) return false; //identical!

  //use loop over neighbors; the thirdneighbor array is too memory costly
  unsigned int atomA, atomB;
  
  if ( frodaAtom.at(atom1).neighbors.size() < frodaAtom.at(atom2).neighbors.size() ) {
    atomA = atom1; atomB = atom2;         
  }
  else {
    atomA = atom2; atomB = atom1;         
  }

  bool maybe = false;
  for ( unsigned int whichN = 0; whichN < frodaAtom.at(atomA).neighbors.size();
        whichN++ ) {
    if ( frodaAtom.at(atomA).neighbors.at(whichN) == atomB ) return false;
    //is a first neighbor               
    unsigned int theN = frodaAtom.at(atomA).neighbors.at(whichN);
    for ( unsigned int whatN = 0; whatN < frodaAtom.at(theN).neighbors.size();
          whatN++ ) {
      if ( frodaAtom.at(theN).neighbors.at(whatN) == atomB ) return false;
      // is a second neighbor
                 
      for ( unsigned int yetN = 0; yetN < frodaAtom.at(atomB).neighbors.size();
            yetN++ ) {
        if ( frodaAtom.at(theN).neighbors.at(whatN) ==
             frodaAtom.at(atomB).neighbors.at(yetN) ) {
          maybe = true; //could be a third neighbor!                                            
        }            
      }                 
                 
    }
  }

  return maybe;

}

#endif
