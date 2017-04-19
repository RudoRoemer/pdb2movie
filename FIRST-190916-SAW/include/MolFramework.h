#ifndef _MOL_FRAMEWORK_
#define _MOL_FRAMEWORK_

#include "global_defs.h"
#include "AllSiteSelector.h"
#include "Structure_Info.h"

#include "SiteSelector.h"

#include "VectorAlgebra.h"
#include "Chain.h"
#include "Residue.h"

////////////////////////////////////////////////////////////////////////////
class MolFramework: public Structure_Info {

public: 

  MolFramework() : Structure_Info() {
    chain_counter = 1;
    prev_topology = NULL;
    total_hbonds = 0;
    total_hphobes = 0;
    total_residues = 0;
    total_stacked_rings = 0;
  };
  ~MolFramework();
  
  inline void  setListAsPruned( vector<unsigned int> atom_list );

  void  addBond( int this_bond );
  void  addConstraint( unsigned int atom_1 = 0 );
  void  addUserDefinedConstraints();
  void  assignAtomsTo3DGrid();
  void  countNeighbors();
  void  assignHydrogenBondStatus();
  void  checkForPreviousBond( int, int, int );
  void  compareTopologies();
  void  computeMeanCoordination( bool pruned = false );
  void  findAllAtomTriples( int, vector<angle_triple> * );
  void  getUnitNormalToPlane( float *, int, int, int );
  void  identifyCovalentBonds();
  void  identifyHydrogenBonds();
  void  identifyHydrophobicTethers();
  void  initialize( int total_sites );
  void  outputCovalentBondList();
  void  outputHydrogenBondList( bool full_output = false );
  void  outputHydrophobicTetherList();
  void  outputStackedRingsList();
  void  listSideChainAtoms( int current_atom,
                               int current_residue,
                               int current_chain_ID );
  
  void  MakeTmpFileForHBDilute();
  void  outputBondData();
  void  outputCONECTRecords( ostream &f );
  void  outputFlexibilityIndexPDBFormat();
  void  outputHARLEMFormat();
  void  outputPyMolScriptCMD_PDBFormat();
  void  outputPyMolScriptRCD_MOL2Format();
  void  outputPyMolScriptRCD_PDBFormat();
  void  output_rasmol_rcd();
  void  outputRawData();
  void  outputRCD_MOL2FORMAT();
  void  outputRCD_PDBFormat( string file_extension, 
                               bool RCD_in_occu = false );
  void  outputText();
  void  prepareFrameworkForAnalysis();
  void  prepareFrameworkForFreeEnergyDilution();
  void  printInfo();
  void  printInfoLine1( int, int &, int position = 0,
                           string bond_type = "CB", 
                           string demark = "" );
  void  pruneDanglingEnds();
  void  pruneSideChains();
  void  queryNetwork();
  void  queryOptionsScreen();
  void  querySiteNumber( vector<unsigned int> atom_number, 
                           bool interresidue_bonds_only = false );
  void  querySequenceNumber( int residue_number = 0 );
  void  readHydrogenBondList();
  void  readHydrophobicTetherList();
  void  read_ringstack_file();
  void  readCovalentBondList();
  void  removeConstraint( int, int, int );
  void  removeBond( int this_bond );
  void  resetSearchList();
  void  resetCheckedLabel();
  void  saveTopology( unsigned int hbond_number );
  void  userInteraction2( vector<int> mult_side_chain_conf_list, 
                            vector<string> alt_loc_label_list, 
                            Site_Info *site_info );
  void  userInteractionEnergyCutoff();
  void  user_interaction_6( int hydrogen );
  void  validateStructure();
  void  viewBondDetails( int, int, int, float );

  bool  isHydrogenBondAcceptor( unsigned int );
  bool  isBonded( SiteID atom1, SiteID atom2 );
  bool  isConnectedToOnlyHydrophobicAtoms( unsigned int );
  bool  isDifferentResidue(SiteID current_atom,
                             SiteID neighbor_atom );
  bool  isDonorHydrogen( int, int * );
  bool  isDoubleBond( SiteID atom_1, 
                        SiteID atom_2 );
  bool  isHydrogenBond( unsigned int, unsigned int );
  bool  isHydrophobicAtom( unsigned int );
  bool  isIsolated( unsigned int );
  bool  isOnlyBondedToHydrogen( unsigned int atom_1 );
  bool  isSaltBridge( unsigned int, unsigned int, unsigned int );
  bool  isTerminalNitrogen( SiteID current_atom );
  bool  isWithinHydrophobicCutoff( unsigned int, unsigned int );
  bool  noSideChainBonds();
  bool  oneInterresidueHydrophobicTetherPerAtom( unsigned int current_atom,
                                          unsigned int neighbor_atom );

  // AJR 10.25.05  new nucleic acid functions
  void  list_aromatics(); // AJR 10.25.05 new output listing for stacked bases
  void  identify_stacked_rings(); // AJR 10.25.05
  void  get_ring_center( float *,int, int, int, int);// AJR 10.31.05
  void  define_ring_center( aromatic_ring &ring_1); // AJR 10.31.05
  void  list_ring_atoms( int, int, int);
  void  list_ring_atoms2();
  bool  is_sugar_phosphate_bond( unsigned int, unsigned int );
  bool  is_ribose_atom( unsigned int atom_1 );
  bool  is_nucleobase_atom( unsigned int atom_1 );
  bool  is_modifiedNucleobase_atom( unsigned int atom_1 );
  bool  is_ring_atom( unsigned int );
  bool  is_stacked_ring_neighbor(int, aromatic_ring );
  bool  is_stacked_ring_pair(aromatic_ring, aromatic_ring);
  bool  is_nucleic_backbone(unsigned int );
  bool  in_ring_stacked_pair(int, int );
  int   get_ringid( unsigned int, int, aromatic_ring * );
  int   is_ring_anchor( unsigned int );
  int   get_ring_atom( unsigned int, unsigned int );
  //--------------------------------------- end of aromatic_ring functions

  short int isBackbone( unsigned int atom_1 );
  short int isMainchain( unsigned int atom_1 );
  bool isNucleicAcidAtom(SiteID atom_1);

  int  alreadyBonded( unsigned int, unsigned int );
  int  determineNitrogenHydrogenBondStatus( unsigned int );
  int  determineOxygenHydrogenBondStatus( unsigned int );
  int  determineSulfurHydrogenBondStatus( unsigned int );
  int  getFIRSTNumber( int );
  int  gridPosition( int current_atom );
  int  gridPosition( int grid_X, int grid_Y, int grid_Z );
  int  getHistidineProtonationState( SiteID );
  int  isChargedResidue( SiteID atom_1 );
  unsigned int  NthNearestNeighbor( unsigned int, 
				    unsigned int, 
				    int, 
				    unsigned int root_node = 0 );
  
  unsigned int  NthNearestNeighbor( unsigned int,
				    string, 
				    int, 
				    unsigned int root_node = 0 );
  
  int  removeFromBondList( unsigned int,
			   unsigned int, 
			   vector<new_bonds> &bond_list );
  
  int  userInteraction3( SiteID atom_1, 
			 SiteID atom_2,
			 Site_Info *site_info,
			 float distance );
  
  float  computeAngle( SiteID atom_1, 
		       SiteID atom_2,
		       SiteID atom_3 );
  
  float  computeOutOfPlaneAngle( angle_triple *donor_list, 
				 angle_triple *accpt_list );

  float  getDistance( SiteID atom_1, 
		      SiteID atom_2 );
  
  float  saltBridgeEnergy( int, int );
  float  hydrogenBondEnergy( int, int );
  
  string atomNamePDBFormat( const string &atom_name, 
			    const string &element );
  string getElementNamePDBFormat( string &atom_name, 
				  bool ambiguous_atom_type = false );
  string getElementName( string &atom_name,
			 bool ambiguous_atom_type = false );
  string getElementName( unsigned int );
  
  static float computePairMSD(MolFramework &, 
                              MolFramework &, 
                              SiteSelector &);
  
  static float computePairRMSD(MolFramework &, 
                               MolFramework &,
                               SiteSelector &);

  // FIXME - refactor to eliminate the need for polymorphism 
  static float computePairMSD(MolFramework &, 
                              vector<Vector> &, 
                              SiteSelector &);
  
  static float computePairRMSD(MolFramework &, 
                               vector<Vector> &,
                               SiteSelector &);
  
  static float computePairMSD(vector<Vector> &, 
                              vector<Vector> &, 
                              SiteSelector &);
  
  static float computePairRMSD(vector<Vector> &,
                               vector<Vector> &, 
                               SiteSelector &);
  
  void stripNonCovalentNeighbors();
  
  //  inline vector<int> getGrids( int current_atom );
  vector<int> getGrids( int current_atom );
  void getGrids( int current_atom, int *grids );
  vector<int>::iterator current_grid;
  vector<new_bonds> hydrogen_bonds;
  vector<new_bonds> hydrophobic_tethers;
  vector<new_bonds> user_defined_constraints;
  vector<new_bonds> stacked_rings; //AJR 10.25.05
  
  map<int,int> orig_2_FIRST_atom_number;
  map<int,int> orig_2_FIRST_chain_ID;
  
  set<int> seq_numbers;
  map<string,int> unique_res_id;
  
  void populateStructureHierarchy();
  map <unsigned int, Chain*> mapFromChainIDtoChain;
  
  map<SiteID, vector<SiteID> > conect_records;

  map <Chain*, unsigned int> mapFromChainToChainID;
  map <Residue*, Chain*> mapFromResidueToChain;
  map <Residue*, unsigned int> mapFromResidueToResidueID;
  map <SiteID, Residue*> mapFromSiteIDToResidue;
      
 private:
  
  int chain_counter;
  int total_hbonds;
  int total_hphobes;
  int hbond_number;
  int total_stacked_rings; // AJR 10.25.05

  stringstream hbond_data;
  stringstream included_hbonds;
  stringstream included_hphobes;
  stringstream included_aromats; // AJR 11.11.05

  vector<int> search_list;
  vector<unsigned int> side_chain_list;
  vector<unsigned int> ring_atom_list; //AJR 10.25.05

  vector<unsigned int>::iterator neighbor_atom;


  Site_Info *prev_topology;
  
};
#endif
