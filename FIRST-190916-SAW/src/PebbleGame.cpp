#include "global_defs.h"
#include "Parameters.h"
#include "PebbleGame.h" 

extern const Parameters parameters;


////////////////////////////////////////////////////////////////////////////////
// Description:
//   PebbleGame class constructor
////////////////////////////////////////////////////////////////////////////////
PebbleGame::PebbleGame( MolFramework *molFramework ){
  this->molFramework = molFramework;
  this->site_info = molFramework->site_info;
  this->total_sites = molFramework->total_sites;
  floppy_modes = 0;
  total_bonds = 0;
  isolated_sites = 0;
  isolated_dimers = 0;
  total_clusters = 0;
  total_stressed_regions = 0;
  total_clusters_larger_than_one = 0;
  total_clusters_larger_than_min_cluster_size = 0;
  
  initializePebbleGame();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   PebbleGame class descructor
////////////////////////////////////////////////////////////////////////////////
PebbleGame::~PebbleGame(){
  for( SiteID a = 0; a <= total_sites; a++ )
    delete [] pebbles[a];
  
  delete [] pebbles;
  delete [] mult;
  delete [] back_track;
  delete [] blocked_label;
  delete [] sites_to_check;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Build up the network, and run the pebble game. Constraints are only added to 
//   the network if both sites have the prune label is set to "false". 
// Parameters:
//   run_decomps - optional parameter to not run the decompositions. Set to "true"
//     by default.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::runPebbleGame( bool run_decomps ){

  if( parameters.verbose &&
      !parameters.bond_dilution )
    cout << "Running pebble game..." << endl;

  unsigned int neighbor_site = 0;
  int bars_to_use = 0;

  // Assign DOF for each site before building the network.
  //////////////////////////////////////////////////////////////////////
  for( SiteID current_site = 1; current_site <= total_sites; current_site++ ){
    if( site_info[current_site].neighbor_list.size() > 1 ){
      pebbles[current_site][PEBBLE_1] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_2] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_3] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_4] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_5] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_6] = FREE_PEBBLE;
    }
    else if( site_info[current_site].neighbor_list.size() == 1 ){
      pebbles[current_site][PEBBLE_1] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_2] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_3] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_4] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_5] = FREE_PEBBLE;
    }
    else{
      pebbles[current_site][PEBBLE_1] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_2] = FREE_PEBBLE;
      pebbles[current_site][PEBBLE_3] = FREE_PEBBLE;
    }
  }

  for( SiteID current_site = 1; current_site <= total_sites; current_site++ ){

    if( !(site_info[current_site].pruned) ){
      
      // 1. Each site starts off with 3 degrees of freedom.
      //////////////////////////////////////////////////////////////////////
      floppy_modes += 3;
      
      // 2. For each bond incident with the current_site, add it to the network if
      //    the site number of the neighbor atom is larger than the current_site. 
      //    This prevents double counting bonds. 
      //////////////////////////////////////////////////////////////////////
      for( unsigned int a = 0; a < site_info[current_site].neighbor_list.size(); a++ ){
	neighbor_site = site_info[current_site].neighbor_list[a];
	
	if( neighbor_site > current_site && 
	    !(site_info[neighbor_site].pruned) ){
	  
	  bars_to_use = site_info[current_site].number_of_bars[neighbor_site];
	  if( bars_to_use >= 0 && bars_to_use <= 6 ){
	    buildNetwork( current_site, neighbor_site, bars_to_use );
	    
	  }
	  else{
	    cout << "Error: The number of bars between sites should be 0-6." << endl;
	    exit(1);
	  }

	  if( parameters.output_bbg )
	    printBBG_Line( current_site, neighbor_site, bars_to_use );

	}
	else if( neighbor_site == current_site) {
	  cout << " Error: current site = neighbor site " << current_site << " " << neighbor_site << endl;
	}
	neighbor_site++;
      }

    }

  }

  // Only run these decompositions when analyzing the complete network, not
  // pruned network or during dilution. 
  //////////////////////////////////////////////////////////////////////
  if( run_decomps ){
    
    // Rigid region decomposition
    //////////////////////////////////////////////////////////////////////
    if( parameters.verbose >= 2 && 
	!parameters.bond_dilution )
      cout << "Rigid cluster decomposition started." << endl;
    rigidClusterDecomposition();
    if( parameters.verbose >= 2 &&
	!parameters.bond_dilution )
      cout << "Rigid cluster decomposition finished." << endl;

    // TODO specify which combination of decompositions to do as parameter
    //      to run_pebble_game, or take these out of this function, and call
    //      them explicitly in general_utils..
    //if( parameters.bond_dilution ) 
    //return;

    // Stressed region decomposition
    //////////////////////////////////////////////////////////////////////
    if( parameters.verbose >= 2)
      cout << "  Stressed region decomposition started." << endl;
    stressedRegionDecomposition();
    if( parameters.verbose >= 2)
      cout << "  Stressed region decomposition finished." << endl;

    // Collective motion decomposition.
    //////////////////////////////////////////////////////////////////////
    if( parameters.verbose >= 2 )
      cout << "  Collective mode decomposition started." << endl;
    dependentHingesDecomposition();
    if( parameters.verbose >= 2 )
      cout << "  Collective mode decomposition finished." << endl;
  }

  // Assign the results of the pebble game decompositions to the input
  // structure. 
  //////////////////////////////////////////////////////////////////////
  molFramework->floppy_modes = floppy_modes;
  molFramework->total_bonds = total_bonds;
  molFramework->total_clusters = total_clusters;
  molFramework->total_stressed_regions = total_stressed_regions;
  molFramework->isolated_sites = isolated_sites;
  molFramework->isolated_dimers = isolated_dimers;
  molFramework->total_clusters_larger_than_one = total_clusters_larger_than_one;
  molFramework->total_clusters_larger_than_min_cluster_size = total_clusters_larger_than_min_cluster_size;

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   For each bond added to the network, this routine will augment the number of 
//   pebbles at each site, if necessary, and then rearrange the pebbles to 
//   determine if the bond is independent. If the bond is between two sites 
//   that are in the same stressed region, then the bond is redundant, and no 
//   change in the pebble covering or the number of floppy modes occurs.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::buildNetwork( SiteID site1, SiteID site2, int bars ){

  //cout << "adding " << site1 << " " << site2 << " " << bars << endl;

  // If the bonds has no bars to place between the sites, do nothing.
  //////////////////////////////////////////////////////////////////////
  if( !bars ){
    cout << "Warning: Constraint between " << site1 << " and " << site2 << " has 0 bars." << endl;
    return;
  }

  // If site1 is not yet connected to anything
  //////////////////////////////////////////////////////////////////////
  if( mult[site1] == 0 ) {

    for( int a = PEBBLE_1; a < bars; a++ )
      pebbles[site1][a] = site2;
    
    if( mult[site2] > 1 )
      floppy_modes += (2 - bars); // 2 new pebbles - # of bars
    
    else if( mult[site2] == 1 )
      floppy_modes += (3 - bars); // 3 new pebbles - # of bars

    else
      floppy_modes += (4 - bars); // 4 new pebbles - # of bars 
  }   
   
  // If site2 is not yet connected to anything (but site1 is)
  //////////////////////////////////////////////////////////////////////
  else if( mult[site2] == 0 ){

    for( int a = PEBBLE_1; a < bars; a++ )
      pebbles[site2][a] = site1;

    if( mult[site1] == 1 )
      floppy_modes += (3 - bars); // 3 new pebbles - # of bars 

    else
      floppy_modes += (2 - bars); // 2 new pebbles - # of bars
  }

  // If both site1 and site2 are already in the network, proceed if they
  // do not belong to the same stressed region
  //////////////////////////////////////////////////////////////////////
  else{
    if( getStressedRegionLabel(site1) != getStressedRegionLabel(site2) ){
      
      if( mult[site1] == 1 )
	floppy_modes++;
      
      if( mult[site2] == 1 )
	floppy_modes++;
      
      coverBonds( site1, site2, bars );
    }
  }

  mult[site1]++;
  mult[site2]++;
  total_bonds++;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the label of the stressed region to which a site belongs. The stress
//   label for a region is set to the smallest site number in the region. The
//   loop below keeps backtracking through stress labels until it finds one
//   where stress_label == site_number. This should (but not necessarily) be
//   the smallest site number. Once the complete pebble game has been played, 
//   a single unique number will be assi 
//   Initially, the stressed region label for every site is set equal to the 
//   site number.
////////////////////////////////////////////////////////////////////////////////
int PebbleGame::getStressedRegionLabel( SiteID site ){

  unsigned int label = site_info[site].stressed_label;

  while( site_info[label].stressed_label != label ){
    label = site_info[label].stressed_label;
  }

  return label;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Rearrange the pebbles in the current network to see if the bond being 
//   added to the network is independent or not. It begins by moving 6 pebbles 
//   to site1, which is always possible. Then we look for N pebbles to place 
//   at site2. If we find N, the bond is independent. If we can not find N 
//   pebbles, this bond creates a stressed region, and all the sites in this 
//   region (laman subgraph) labeled accordingly.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::coverBonds( SiteID site1, SiteID site2, int total_bars ){

  int 
    pebbles_needed = total_bars,
    sites_checked = 0,
    current_pebble = 0,
    uncovered_bars = 0;

  unsigned int
	  smallest_site_number = 0,
	  old_label = 0,
	  new_label = 0,
	  current_site = 0;

  // 1. Collect 6 pebbles on site1. This is always possible
  //////////////////////////////////////////////////////////////////////
  collectGroupMaxPebbles( site1 );

  // 2. Determine how many pebbles we need on site2. 
  //////////////////////////////////////////////////////////////////////
  for( int a = PEBBLE_1; a <= PEBBLE_6; a++ ){
    if( pebbles[site2][a] == FREE_PEBBLE )
      pebbles_needed--;
  }

  // 3. Start looking for free pebbles in the network
  //////////////////////////////////////////////////////////////////////
  back_track[site2] = -1;
  for( int a = 1; a <= pebbles_needed; a++ ){
    
    blocked_cutoff++;
    blocked_label[site1] = blocked_cutoff;
    blocked_label[site2] = blocked_cutoff;

    findPebble( site2 );

    // A. If we fail a pebble search
    if( failed_pebble_search ){

      //cout <<"Failed a pebble search" << endl;
      failed_pebble_search = false;

      sites_to_check[0] = site1;
      sites_to_check[1] = site2;
      smallest_site_number = site1;

      floppy_modes = floppy_modes -total_bars +pebbles_needed -a +1;

      if( total_sites_found < MIN_SEARCH_RANGE ){

	for( int b = 0; b < total_sites_found; b++ ){

	  current_site = sites_to_check[b];
	  while( site_info[current_site].stressed_label < current_site ) {
	    current_site = site_info[current_site].stressed_label;
	  }
                       
	  if( current_site < smallest_site_number )
	    smallest_site_number = current_site;
	}
      }
      else{

	sites_checked = 0;
	while( sites_checked < total_sites_found ){
	  current_site = sites_to_check[sites_checked];

	  while( site_info[current_site].stressed_label < current_site ){
	    current_site = site_info[current_site].stressed_label;
	    if( blocked_label[current_site] < blocked_cutoff )
	      expandStressedRegion( current_site );
	  }
	  
	  if( current_site < smallest_site_number )
	    smallest_site_number = current_site;

	  sites_checked++;
	}
      }

      // assign stressed region labels to sites in the search
      //////////////////////////////////////////////////////////////////////
      new_label = site_info[site1].stressed_label;
      site_info[site1].stressed_label = smallest_site_number;
      if( new_label < site1 ){
	do{
	  old_label = new_label;
	  new_label = site_info[old_label].stressed_label;
	  site_info[old_label].stressed_label = smallest_site_number;
	} while( new_label < old_label );
      }

      for( int b = 1; b < total_sites_found; b++ ){
	current_site = sites_to_check[b];

	new_label = site_info[current_site].stressed_label;
	site_info[current_site].stressed_label = smallest_site_number;
	if( new_label < site1 ){
	  do{
	    old_label = new_label;
	    new_label = site_info[old_label].stressed_label;
	    site_info[old_label].stressed_label = smallest_site_number;
	  } while( new_label < old_label );
	}

	// point all pebbles in the stressed region to site1, this will make
	// subsequent pebble searches in this region very fast. 
	//////////////////////////////////////////////////////////////////////
	for( int c = 0; c < 6; c++ )
	  pebbles[current_site][c] = site1;
      }
     
      return;

    }
  }

  // This next section of code is only executed if the pebble search did 
  // not fail for site2
  //////////////////////////////////////////////////////////////////////
  floppy_modes -= total_bars;

  uncovered_bars = total_bars;
  current_pebble = 0;

  while( uncovered_bars ){

    if( pebbles[site2][current_pebble] == FREE_PEBBLE ){
      pebbles[site2][current_pebble] = site1;
      uncovered_bars--;
    }
    current_pebble++;

    if( current_pebble > 6 ){
      cout << "Error in cover_bonds routine. Insufficient pebbles found" << endl;
      exit(80);
    }
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   If any of 
//   the pebbles associated with site1 are covering bonds, this routine will 
//   follow the bond, and any necessary subsequent bonds (in a breadth-first 
//   search) to find a site with a free pebble. This free pebble is then places 
//   on the last bond searched, which displaces the pebble currently on the bond, 
//   which is used to displace the previous bond, and so on, until the pebble 
//   on the initial bond searched is displaced onto site1, which will now have 
//   one additional free pebble. This process is repeated until site1 has 6 
//   free pebbles.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::collectGroupMaxPebbles( SiteID site1 ){

  int pebbles_needed = 0;

  for( int a = 0; a < 6; a++ ){
    if( pebbles[site1][a] > 0 ) 
      pebbles_needed++;
  }

  back_track[site1] = -1;
  
  for( int a = 1; a <= pebbles_needed; a++ ){
    blocked_cutoff++;
    blocked_label[site1] = blocked_cutoff;

    findPebble( site1 );
  }

  if( failed_pebble_search ){
    cout << "ERROR: Failed to collect " << pebbles_needed 
	 << "pebbles on site " << site1 << ". (" << site_info[site1].orig_atom_number << ")" << endl;
    exit(1);
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Begining with bonds covered by the site "search_root", traverse the current 
//   network (all bonds have not necessarily been placed yet) until a free pebble 
//   is found. Keep track of each site that has been searched. The routine uses a 
//   depth-first search of the network. Since a record of the search is kept, when
//   a free pebble is found, place it on the last bond searched, and successively 
//   displace-replace pebbles from each site in the search until we reach the 
//   "search_root", and let the free pebble stay there.
// Parameters:
//   root_node - site where the pebble search will begin. 
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::findPebble( int root_node ){

  int current_pebble = 0;
  int sites_checked = 1;
  int node = 0;
  int previous_node = 0;

  total_sites_found = 2; // need to count the two sites of the current bond
 
  // store each unique site that is connected to root_node AND has
  // a pebble from root_node covering the intervening bond. 
  //////////////////////////////////////////////////////////////////////
  for( int a = 0; a < 6; a++ ){
    current_pebble = pebbles[root_node][a];
    if( current_pebble != FREE_PEBBLE &&
	blocked_label[current_pebble] < blocked_cutoff ){

      sites_to_check[total_sites_found] = current_pebble;
      total_sites_found++;
      blocked_label[current_pebble] = blocked_cutoff;
      back_track[current_pebble]  = root_node;

    }
  }

  // the array "sites_to_check[]" now contains all the sites that root_node 
  // points to with one of its pebbles (the pebble graph is a directed graph). 
  //////////////////////////////////////////////////////////////////////
  sites_checked = 2; // skip the first two positions.
  while( sites_checked < total_sites_found ){
    
    node = sites_to_check[sites_checked];
    for( int a = 0; a < 6; a++ ){
      current_pebble = pebbles[node][a]; 
 
      if( current_pebble != FREE_PEBBLE && // this is an unchecked site, store it to check later
	  blocked_label[current_pebble] < blocked_cutoff ){ 
	sites_to_check[total_sites_found] = current_pebble;
	total_sites_found++;
	blocked_label[current_pebble] = blocked_cutoff; // now this site is blocked
	back_track[current_pebble]  = node; 
      }
      else{ // we found a free pebble, or this site (current_pebble) is blocked.
	if( current_pebble == FREE_PEBBLE ){
	  total_sites_found = -1;

	  // When we find the free pebble, place it on the bond, pointing at the previous node,
	  // then begin backtracking (and displacing/replacing pebbles) until we reach root_node.
	  ////////////////////////////////////////////////////////////////////////////////
	  previous_node = back_track[node];
	  pebbles[node][a] = previous_node;
	 
	  do{
	    for( int a = 0; a < 6; a++ ){
	      if( pebbles[previous_node][a] == node ){
		node = previous_node;
		previous_node = back_track[node];
		pebbles[node][a] = previous_node;
		break;
	      }
	    }
	  } while( previous_node != -1 ); // back_track = -1 for the first site in the list, ie. root_node.

	  return;
	}

      }
    }

    sites_checked++;
  }

  failed_pebble_search = true;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   If a failed pebble search covers more then MIN_SEARCH_RANGE sites, then 
//   check each site to see if it already belongs to a stressed region, but 
//   wasn't labeled yet.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::expandStressedRegion( int root_node ){

  int 
    current_site = 0,
    current_pebble = 0,
    sites_checked = 0;
  
  sites_checked = total_sites_found;
  
  sites_to_check[total_sites_found] = root_node;
  blocked_label[root_node] = blocked_cutoff;
  total_sites_found++;
  
  while( sites_checked < total_sites_found ){
    
    current_site = sites_to_check[sites_checked];
    current_pebble = pebbles[current_site][PEBBLE_1];
    
    if( current_pebble != FREE_PEBBLE &&
	blocked_label[current_pebble] < blocked_cutoff ){
      sites_to_check[total_sites_found] = current_pebble;
      blocked_label[current_pebble] = blocked_cutoff;
      total_sites_found++;
    }
    sites_checked++;
  }
  
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void PebbleGame::initializePebbleGame(){

  // Create the pebble array and initialize each site to have 3 free pebbles.
  // Assign the rigid and stressed region label for each site equal to the site number.
  ////////////////////////////////////////////////////////////////////////////////
  pebbles = new int *[total_sites +1];   

  for( SiteID siteNumber = 0; siteNumber <= total_sites; siteNumber++ ){
    pebbles[siteNumber] = new int[6];
    
    for( int b = 0; b < 3; b++ )
      pebbles[siteNumber][b] = FREE_PEBBLE;
    for( int b = 3; b < 6; b++ )
      pebbles[siteNumber][b] = 0;
  }

  mult = new int[total_sites +1];
  back_track = new int[total_sites +1];
  blocked_label = new int[total_sites +1];
  sites_to_check = new int[total_sites +1];
  
  for( SiteID siteNumber = 0; siteNumber <= total_sites; siteNumber++ ){
    mult[siteNumber] = 0;
    back_track[siteNumber] = 0;
    blocked_label[siteNumber] = 0;
    sites_to_check[siteNumber] = 0;
  }
  
  blocked_cutoff = 3;
  total_sites_found = 0;
  failed_pebble_search = 0;
  
  is_floppy = 0;
  isolated_sites = 0;
  isolated_dimers = 0;
  sites_checked = 0;

  total_bonds = 0;
  total_clusters = 0;
  total_clusters_larger_than_one = 0;
  total_stressed_regions = 0;
  total_hinges = 0;
};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::resetPebbleGame(){

  total_bonds  = 0;
  floppy_modes = 0;

  int b = 0;
  for( SiteID siteNumber = 0; siteNumber <= total_sites; siteNumber++ ){

    for( b = 0; b < 3; b++ )
      pebbles[siteNumber][b] = FREE_PEBBLE;

    for( b = 3; b < 6; b++ )
      pebbles[siteNumber][b] = 0;

    mult[siteNumber] = 0;
    back_track[siteNumber] = 0;
    blocked_label[siteNumber] = 0;
    sites_to_check[siteNumber] = 0;
  }
  
  blocked_cutoff = 3;
  total_sites_found = 0;
  failed_pebble_search = 0;
  
  is_floppy = 0;
  isolated_sites = 0;
  isolated_dimers = 0;
  sites_checked = 0;
  
  total_clusters = 0;
  total_clusters_larger_than_one = 0;
  total_stressed_regions = 0;
  total_hinges = 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::printInfo(){

  neighbor_site = (site_info[3].neighbor_list).begin();

  while( neighbor_site != (site_info[3].neighbor_list).end() ){

    cout << "   " << *neighbor_site;
    neighbor_site++;
  }
  cout << endl;

  cout <<  pebbles[38][0] << " "  <<  pebbles[38][1] << " " <<  pebbles[38][2] << endl;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Debugging function for the "PebbleGame" class. Prints the current state 
//   of the pebble network. One line for each site, one column for each pebble.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::printPebbleCovering(){

  ofstream out;

  out.open("pebble_covering.txt", ios::app );
  out << "total bonds " << total_bonds << endl;
  out << "total sites " << total_sites << endl << endl;
  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ )
    out << setw(8) << siteNumber 
	<< setw(8) << pebbles[siteNumber][0] << setw(8) << pebbles[siteNumber][1] 
	<< setw(8) << pebbles[siteNumber][2] << setw(8) << pebbles[siteNumber][3] 
	<< setw(8) << pebbles[siteNumber][4] << setw(8) << pebbles[siteNumber][5] << endl;
  
  out.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Write the current constraint to a file for analysis by the pebble3D program.
//   Used primarily for verification and debugging, the pebble3D program only
//   allows constraints of 5 bars. Therefore, it is necessary to construct 
//   additional network elements in order to reproduce the proper constraint 
//   network.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::printBBG_Line( SiteID site_1, SiteID site_2, int bars){

  static unsigned int extra_site = total_sites +1;
  ofstream test("bonds.dat", ios::app );

  if( !parameters.all_5_bars ){
    test << setw(8) << site_1 
	 << setw(8) << site_2 
	 << setw(8) << bars << endl;
  }
  else{
    if( bars == 6 ){
      test << setw(8) << site_1 << setw(8) << extra_site << endl;
      test << setw(8) << extra_site << setw(8) << extra_site+1 << endl;
      extra_site++;
      test << setw(8) << extra_site << setw(8) << extra_site+1 << endl;
      extra_site++;
      test << setw(8) << extra_site << setw(8) << extra_site+1 << endl;
      extra_site++;
      test << setw(8) << extra_site << setw(8) << site_2 << endl;
      extra_site++;
      test << setw(8) << site_1 << setw(8) << site_2 << endl;
    }
    else if( bars == 5 )
      test << setw(8) << site_1 << setw(8) << site_2 << endl;
    else if( bars == 4 ){
      test << setw(8) << site_1 << setw(8) << extra_site << endl;
      test << setw(8) << extra_site++ << setw(8) << site_2 << endl;
    }
    else if( bars == 3 ){
      test << setw(8) << site_1 << setw(8) << extra_site << endl;
      test << setw(8) << extra_site++ << setw(8) << extra_site << endl;    
      test << setw(8) << extra_site++ << setw(8) << site_2 << endl;
    }
    else if( bars == 2 ){
      test << setw(8) << site_1 << setw(8) << extra_site << endl;
      test << setw(8) << extra_site << setw(8) << extra_site+1 << endl;    
      extra_site++;
      test << setw(8) << extra_site << setw(8) << extra_site+1 << endl;
      extra_site++;
      test << setw(8) << extra_site << setw(8) << site_2 << endl;
      extra_site++;
    }
    else if( bars == 1 ){
      test << setw(8) << site_1 << setw(8) << extra_site << endl;
      test << setw(8) << extra_site++ << setw(8) << extra_site << endl;    
      test << setw(8) << extra_site++ << setw(8) << extra_site << endl;
      test << setw(8) << extra_site++ << setw(8) << extra_site << endl;
      test << setw(8) << extra_site++ << setw(8) << site_2 << endl;
    }
    else{
      cout << "ERROR: No bars for this bond!" << endl;
      exit(1);
    }
  }

}
////////////////////////////////////////////////////////////////////////////////

