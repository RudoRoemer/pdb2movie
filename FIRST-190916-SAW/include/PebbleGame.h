#ifndef _PEBBLE_GAME_
#define _PEBBLE_GAME_

#include "global_defs.h"
#include "MolFramework.h"

////////////////////////////////////////////////////////////////////////////////
class PebbleGame {

public:

  PebbleGame( MolFramework *molFramework );
  ~PebbleGame();
  
  // These functions are in PebbleGame.cpp
  //////////////////////////////////////////////////////////////////////
  void runPebbleGame( bool run_decomps );
  void initializePebbleGame();
  void resetPebbleGame();
  void printBBG_Line( SiteID site_1, SiteID site_2, int bars );
  void buildNetwork( unsigned int, unsigned int, int );
  void coverBonds( unsigned int, unsigned int, int );
  void collectGroupMaxPebbles( unsigned int );
  void findPebble( int );
  int  getStressedRegionLabel( unsigned int );
  void expandStressedRegion( int );
  void printPebbleCovering();
  void printInfo();

  // These functions are in PebbleGameDecompRigid.cpp
  //////////////////////////////////////////////////////////////////////
  void rigidClusterDecomposition();
  void saveThisSiteToSearchList( int );
  void initializeTopologySearchAtSite( int );
  void initializePebbleSearchAtSite( int );
  void startPebbleSearchAtSite( int );
  void resetBlockedLabels();

  // These functions are in PebbleGameDecompStressed.cpp
  //////////////////////////////////////////////////////////////////////
  void stressedRegionDecomposition();

  // These functions are in PebbleGameDecompHinges.cpp
  //////////////////////////////////////////////////////////////////////
  void dependentHingesDecomposition();
  int  lockHinge( unsigned int, unsigned int, int );

private:

  int total_included_sites;
  vector<SiteID>::iterator  neighbor_site;
  
  // Pebble game variables
  //////////////////////////////////////////////////////////////////////
  int **pebbles;
  int *mult;
  int *back_track;
  int *blocked_label;
  int *sites_to_check;
  int blocked_cutoff;
  int total_sites_found;
  bool failed_pebble_search;

  // Rigid cluster decomposition variables
  //////////////////////////////////////////////////////////////////////
  int is_floppy;
  int sites_checked;
  int pebbles_followed;
  int pebble_sites_searched;

  // Collective motion variables
  //////////////////////////////////////////////////////////////////////
  unsigned int total_hinges;
  int *hinge_label;
  int **hinges;
  int **dihedral_angle;

  MolFramework *molFramework;
  Site_Info *site_info;

  // Copy local values in the current MolFramework structure
  //////////////////////////////////////////////////////////////////////
  int floppy_modes;
  unsigned int total_sites;
  int total_bonds;
  int isolated_sites;
  int isolated_dimers;
  unsigned int total_clusters;
  unsigned int total_stressed_regions;
  unsigned int total_clusters_larger_than_one;
  unsigned int total_clusters_larger_than_min_cluster_size;

};
#endif

