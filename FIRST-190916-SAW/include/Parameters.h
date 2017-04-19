#ifndef _PARAMETERS_
#define _PARAMETERS_

#include "global_defs.h"
using namespace std;

enum XAxisType {perSite, perResidue, perRigidLabel};

////////////////////////////////////////////////////////////////////////////////
// Description
//   This data structure holds all the free parameters for running FIRST. Most of 
//   of them dictate various run-type options, but there are also parameters that
//   affect the analysis, such as the geometric criteria for the hydrogen bond
//   energy function, and the default number of bars to use for modeling different
//   microscopic forces. 
////////////////////////////////////////////////////////////////////////////////
struct Parameters {

  bool interactive; 
  bool exclude_water;
  bool run_timme;
  bool useSpectralColoring;
  bool saveSigma;
  bool diluteWithRigidClusterFactory;
  bool make_raw_data_file;
  
  bool keep_missing_bbg_vertices;
  bool do_not_validate;
  bool run_query;
  bool lock_disulfide;

  bool hbin;
  bool phin;
  bool srin;
  bool read_user_defined_constraint_file;
  bool runFRODA;

  bool using_target;
  bool using_restart;
  bool timme_backbone;
	bool timme_complement;

  bool outputNeighborList;
  bool useLinearTimme;
  bool useLinearCutoffScale;

  float xPeriod;
  float yPeriod;
  float zPeriod;
  
  bool publicationStyle;
  bool covalentTIMME;
  unsigned short timme_nThNeighborForFlexibility;
  unsigned short timme_nThNeighborForPlasticity;
  unsigned short timme_bodyClusterSize;
  float timme_neighborCutoffDistance;
  bool timme_useAssembly;

  float timmeScaleFactorA;

  bool includeHydrogenbonds;
  bool includeHydrophobicTethers;
  bool run_first;

  bool body;
  bool dihedral;
  bool found_alt_side_chain_loc;
  bool output_bbg;
  bool skip_identify_hydrogen_bonds;
  bool skip_identify_hydrophobic_tethers;
  bool skip_identify_aromatics;

  bool fancy_target;
  bool ghost_debug;
  bool every_cycle;

  bool mobileRC1;
  bool centri;
  double centri_maxr;
  bool localCenter;
  vector< int > localCenters;

  bool phout;
  bool hbout;
  bool srout;
  bool outputGhosts;
  bool inputGhosts;

  bool makeSingleMove;

	bool use_first_numbering;
	
  bool covin;
  bool covout;
  bool all_5_bars;
  bool energy_set_on_command_line;

  bool polarGeometry;
  double polarHRadius;

  double phobic_scaling;

  bool flexweb;
  bool runMorph;
  int morphFrames;
  bool use_unpruned_mean_coord_in_stripy_plot;

  bool correction_for_2_4;
  bool energyFxnA;

  short int hphobe_fxn;
  short int run_number;

  short int output_level;
  short unsigned int min_output_cluster_size;

  short int verbose;
  short int hp_bars;

  short int hb_bars;
  short int bond_dilution;

  int cov_bars;
  int ud_bars;
  int rama;
  int output_freq;
  int total_conformations;
  int use_model_number;
  int maxFitCycles;
  int bodyResponseEvery;
  int atleast;
  
  int numberOfFRODAwebOverlaySnapshots;

  float energy_cutoff;
  float PH_cutoff;
  float PH_nucleic_cutoff;
  float resolution_factor;
  float grid_length;
  float maxdist;
  float cutoff_SB_hyd_accpt_dist;
  float cutoff_SB_donor_accpt_dist;
  float cutoff_SB_donor_hyd_accpt_angle;
  float cutoff_SB_hyd_accpt_base_angle;
  float cutoff_HB_hyd_accpt_dist;
  float cutoff_HB_donor_accpt_dist;
  float cutoff_HB_donor_hyd_accpt_angle;

  float step_size;
  //new parameters for stacked aromatic rings:  AJR
  float cutoff_SR_normal_plane_angle;
  float cutoff_SR_center_normal_angle;
  float cutoff_SR_distance;
  float cutoff_SRnucleic_distance;
  short int rs_bars;

  double dtol;
  double dstep;
  double prob;

  double suppress_mc;
  double ghostTol;
  double ph_radius;
  double vdw_overlap;

  string start_time;
  string command_line;
  string alt_location_label;
  string target_name;
  string restart_name;
  string infile_name;
  string path;
  string timme_outputPrefix;
  string path_to_working_dir;
  
  bool cylindrical; 
  double cyl_low_z;
  double cyl_high_z;

  bool local;
  double lstep;

  double vdwNotice;

  // simulated annealing options
  bool doAnneal;
  double annealScale;
  bool constantRate;
  double acceptRate;
  bool doChill;
  int chillSteps;

  float stripeThickness;
  float separationBetweenStripes;
  size_t minimumRigidClusterSize;
  
  bool useSeed;
  unsigned long seed;

  bool propto; // proportional targetting

  bool reportByRSMD;
  double RMSDSpacing;
  
  bool use_group_id; // extra identification field

  vector< pair<string,int> > taskNames;
  map< string, string > mapFromTaskNameToStatus;

  bool scoreMC;
  double frodaMCScale;

  bool doGetArea;
  bool doGetPhobicRg;
  
  bool mcRg;
  double RgWeight;
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

  bool nosterics;

  bool report_resrmsd;

  bool chatty;
  bool bondnew;

  bool useCM;
  bool useCMforward;
  bool useCMreverse;
  double CMforward;
  double CMreverse;
  int CMlower;
  int CMupper;
  double CMmax;

  bool persist;
  int nPersist;

  bool findExposed; //check for surface exposure?
  int exposedBonds; //include exposed bonds?
  
  XAxisType dilutionPlotXAxisType;

  bool useHatbox;
  double hatboxZa;
  double hatboxZb;
  double hatboxXY;

  bool usePairs;
  bool stopWhenAllPairConstraintsSatisfied;

  bool useElasticVector;
  bool useInternalBias;

  bool useTalksto;
  bool changeWeights;

  // flags for electron density options
  bool useED;
  bool useTheoMap;
  bool useEZDMap;
  bool trimMap;
  double edResFac;
  string theoMapFile;
  string ezdMapFile;
  double edTol;
  double edNoise;
  int trimMapFactor;
 
  // SAXS options
  bool useSAXS;
  double saxsTol;
  string saxsFile;
  
  //Froda2 options
  bool froda2;
  bool froda2Hybrid;
  bool froda2UseGhostTol;
  bool froda2Momentum;
  double froda2Tol;
  string froda2RepulsionType;
  string froda2AmberPrmtopFile;
  bool froda2AmberTrajOutput;
  
  bool stopAtEffectiveTimeLimit;
  double effectiveTimeLimit;
  double effectiveTimeFreq;

  //Symmetry-enforcing options
  bool useSymmetry;
  
  // Constructor function for type arg_list. 
  Parameters();
  ~Parameters(){};
  bool isParameter( char *arg );
  bool checkNonBoolParameter( int argNumber, int argc, char **argv, string myArg );
  bool checkBoolParameter( int argNumber, int argc, char **argv, string myArg );
  void readCommandLineOptions( int argc, char **argv );
  void printUsage();
};

#endif
