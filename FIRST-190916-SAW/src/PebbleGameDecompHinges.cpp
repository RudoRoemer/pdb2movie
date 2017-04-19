#include "global_defs.h"
#include "PebbleGame.h" 
#include "generalUtils.h"

const short int UNLOCKED = -1;
const short int LOCKED = 0;

unsigned int *neighbors;
int *pointer;
int *labels_array;

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::dependentHingesDecomposition(){

  int site1 = 0;
  int site2 = 0;
  int bars = 0;
  int new_label = 0;

  float numerator = 0.0;
  float denominator = 0.0;    
  
  dihedral_angle = new int *[total_bonds+1];
  hinge_label = new int[total_bonds+1];

  for( int a = 0; a <= total_bonds; a++ ){
    dihedral_angle[a] = new int[2];
    hinge_label[a] = 0;
  }

  // Map the neighbor_list data into a linear array for fast access.
  ////////////////////////////////////////////////////////////////////////////////
  try{
    neighbors = new unsigned int[2*total_bonds+1];
  } catch (bad_alloc xa){
    cout << "can't allocate neighbors array" << endl;
    exit(1);
  }

  labels_array = new int[2*total_bonds+1];
  pointer = new int[total_sites+1];
  int counter = 0;

  for( unsigned int a = 1; a <= total_sites; a++ ){
    pointer[a] = counter;
    neighbor_site = (site_info[a].neighbor_list).begin();
    while( neighbor_site != (site_info[a].neighbor_list).end() ){
      labels_array[counter] = 0;
      neighbors[counter] = *neighbor_site;
      counter++;
      
      neighbor_site++;
    }
  }

  // Error check for total bonds
  //////////////////////////////////////////////////////////////////////
  if( counter != 2*total_bonds ){
    cout << " Error: Number of total bonds miscalculated in decomp_hinges." << endl;
    exit(1);
  }

  // 1. If the two sites that define a bond are in different rigid clusters, then
  // the given bond can rotate. Find all rotatable bonds, label them, and store
  // the site numbers in the array "dihedral_angle".
  ////////////////////////////////////////////////////////////////////////////////
  int b = 0;
  for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){

    if( (site_info[current_site].neighbor_list).size() > 1 ){    
      b = -1;

      neighbor_site = (site_info[current_site].neighbor_list).begin();
      while( neighbor_site != (site_info[current_site].neighbor_list).end() ){
	b++;

	if( *neighbor_site > current_site && 
	    site_info[current_site].rigid_label != site_info[*neighbor_site].rigid_label ) {
	  
	  total_hinges++;
	  dihedral_angle[total_hinges][0] = current_site;
	  dihedral_angle[total_hinges][1] = *neighbor_site;

	  hinge_label[total_hinges] = UNLOCKED;
	  
	  labels_array[pointer[current_site]+b] = total_hinges;

	  unsigned int c = 0;
	  do{
	    c++;
	  } while( (neighbors[pointer[*neighbor_site] +(c-1)] != current_site) && c < 20);
	  if( c == 20 ){
	    cout << "couldn't find neighbor in decomp_hinges" << endl;
	    exit(1);
	  }
	  labels_array[pointer[*neighbor_site]+c] = total_hinges;
	  
	}
	neighbor_site++;
      }
    }
  }

  int collective_mode_DOF[total_hinges +1];
  int new_collective_mode_DOF[total_hinges +1];
  int collective_mode_edges[total_hinges +1];
  int collective_mode_bars[total_hinges +1];

  for( unsigned int a = 0; a <= total_hinges; a ++ ){
    collective_mode_DOF[a] = 0;
    new_collective_mode_DOF[a] = 0;
    collective_mode_edges[a] = 0;
    collective_mode_bars[a] = 0;
  }

  // Group dependent rotations in collective motions. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int current_hinge = 1; current_hinge <= total_hinges; current_hinge++ ){
  
    // hinges that belong to a failed pebble search in a previous iteration
    // will no longer be labled UNLOCKED.
    //////////////////////////////////////////////////////////////////////
    if( hinge_label[current_hinge] == UNLOCKED ){
      
      site1 = dihedral_angle[current_hinge][0];
      site2 = dihedral_angle[current_hinge][1];
      bars  = site_info[site1].number_of_bars[site2];

      collective_mode_DOF[current_hinge] = lockHinge( site1, site2, bars );
    }
  }

  // Relabel overlapping collective motions. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int a = 1; a <= total_hinges; a++ ){
    new_label = hinge_label[a];

    if( new_label > 0 )
      hinge_label[a] = hinge_label[new_label];
  }

  // Sort and relabeled collective motions based on size. 
  //////////////////////////////////////////////////////////////////////
  vector<vector<unsigned int> > new_labels(total_hinges +1);

  for( unsigned int a = 1; a <= total_hinges; a++ ){
    if( hinge_label[a] > 0 )
      new_labels[ hinge_label[a] ].push_back(a);
  }

  stable_sort( new_labels.begin(), new_labels.end(), sortOnSize );
  
  for( unsigned int a = 1; a <= total_hinges; a++ ){
    for( unsigned int b = 0; b < new_labels[a-1].size(); b++ )
      hinge_label[ new_labels[a-1][b] ] = a;
  }

  // Collect data for computing the flexibility index.
  //////////////////////////////////////////////////////////////////////
  for( unsigned int a = 1; a <= total_hinges; a++ ){

    site1 = dihedral_angle[a][0];
    site2 = dihedral_angle[a][1];
    bars  = site_info[site1].number_of_bars[site2];
    
    collective_mode_edges[ hinge_label[a] ] += 1;

    //cout << "adding bars " << a << " " << site1 << "-" << site2 << " " << bars << endl;
    //cout << "hinge " << a << " " << new_label << " (" << hinge_label[a] << ") " << collective_mode_DOF[a] << " " << collective_mode_DOF[ hinge_label[a] ] << endl;

    if( hinge_label[a] > 0 ){
      collective_mode_bars[ hinge_label[a] ] += bars;
    }

    if( collective_mode_DOF[a] < 0 )
      new_collective_mode_DOF[ hinge_label[a] ] += collective_mode_DOF[a];
  }

  // Label the bonds.
  //////////////////////////////////////////////////////////////////////
  //  ofstream hinges("cmd.txt");
  int CM_label;

  for( unsigned int a = 1; a <= total_sites; a++ ){

    b = -1;
    neighbor_site = (site_info[a].neighbor_list).begin();
    while( neighbor_site != (site_info[a].neighbor_list).end() ){
      b++;

      if( *neighbor_site > a ){
      
	CM_label = labels_array[ pointer[a]+b ];

	site_info[a].cm_label[*neighbor_site] = CM_label;
	site_info[*neighbor_site].cm_label[a] = CM_label;

	if( CM_label == 0 ){
	  //  hinges << setw(8) << a 
	  // << setw(8) << *neighbor_site
	  // << setw(8) << "-1" << " a " << endl;
	  site_info[a].cm_label[*neighbor_site] = -1;
	  site_info[*neighbor_site].cm_label[a] = -1;
	}
	else if( hinge_label[CM_label] == UNLOCKED ){
	  //hinges << setw(8) << a 
	  //<< setw(8) << *neighbor_site
	  //<< setw(8) << " 0" << " b " << endl;
	  site_info[a].cm_label[*neighbor_site] = 0;
	  site_info[*neighbor_site].cm_label[a] = 0;
	  site_info[a].bond_info[*neighbor_site].flexibility_index = 1.0;
	}
	else{
	  //hinges << setw(8) << a 
	  //<< setw(8) << *neighbor_site
	  //<< setw(8) << hinge_label[CM_label] << " c " << endl;
	  site_info[a].cm_label[*neighbor_site] = hinge_label[CM_label];
	  site_info[*neighbor_site].cm_label[a] = hinge_label[CM_label];
	  
	  numerator = abs(new_collective_mode_DOF[ hinge_label[CM_label] ]);
	  denominator = (6 * collective_mode_edges[ hinge_label[CM_label] ])  - collective_mode_bars[ hinge_label[CM_label] ];

	  //cout << "DOF   " << new_collective_mode_DOF[ hinge_label[CM_label] ] << endl;
	  //cout << "edges " << collective_mode_edges[ hinge_label[CM_label] ] << endl;
	  //cout << "bars  " << collective_mode_bars[ hinge_label[CM_label] ] << endl;
	  site_info[a].bond_info[*neighbor_site].flexibility_index = numerator/denominator;
	}
       
      }
      neighbor_site++;
    }
  }
  //  hinges.close();

  // Determine the collective motion labels according the coloring scheme 
  // outlined in the manual.
  ////////////////////////////////////////////////////////////////////////////////
  int neighbor_site = 0;
  int bond_CM_label = 0;
  int first_CM_label = 0;
  short int CM_counter = 0;
  bool has_rigid_bond = false;
  bool multiple_CM_labels = false;
  vector< vector<int> > CMs_incident_with_this_cluster(total_clusters+1);
  
  for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){

    CM_counter = 0;
    first_CM_label = 0;
    has_rigid_bond = false;
    multiple_CM_labels = false;
    
    for( unsigned int a = 0; a < (site_info[current_site].neighbor_list).size(); a++ ){

      neighbor_site = site_info[current_site].neighbor_list[a];

      bond_CM_label = site_info[current_site].cm_label[neighbor_site];

      if( bond_CM_label >= 1 ){
	if( !multiple_CM_labels ){
	  first_CM_label = bond_CM_label;
	}
	else if( first_CM_label != bond_CM_label )
	  multiple_CM_labels = true;
      }
      else if( bond_CM_label == -1 )
	has_rigid_bond = true;
    }

    if( !multiple_CM_labels ){
      site_info[current_site].coll_mode_label = first_CM_label;
    }
    if( has_rigid_bond &&
	first_CM_label != 0 ){
      CMs_incident_with_this_cluster[site_info[current_site].rigid_label].push_back( first_CM_label );      
      // cout << "pushing back " << first_CM_label << "\t total_clusters = "<<total_clusters<<endl;
    }

  }

  // Total how many collective modes are incident with each rigid 
  // cluster. If a rigid cluster is part of only a single collective
  // mode, label all the atoms are belonging to the collective mode.
  // NOTE: Skip rigid cluster 1. As the largest rigid cluster, it will
  //       provide the baseline.
  //////////////////////////////////////////////////////////////////////
  for( unsigned int a = 1; a <= total_clusters; a++ ){
    sort( CMs_incident_with_this_cluster[a].begin(),
	  CMs_incident_with_this_cluster[a].end() );
  
    vector<int>::iterator pos = unique( CMs_incident_with_this_cluster[a].begin(),
					CMs_incident_with_this_cluster[a].end() );
    CMs_incident_with_this_cluster[a].erase(pos, CMs_incident_with_this_cluster[a].end() );
  }    


  for( unsigned int a = 2; a <= total_clusters; a++ ){
    //    cout << "Rigid cluster " << a << " is incident with " << CMs_incident_with_this_cluster[a].size() << " collective modes" << endl;
    if( CMs_incident_with_this_cluster[a].size() == 1 ){
      int new_label = CMs_incident_with_this_cluster[a].at(0);
      for( unsigned int b = 0; b < molFramework->RC_atom_list[a].size(); b++ ){
	site_info[molFramework->RC_atom_list[a].at(b)].coll_mode_label = new_label;
      }
    }
  }

  // Free alocated arrays.
  //////////////////////////////////////////////////////////////////////
  for( int a = 0; a <= total_bonds; a++ ){
    delete [] dihedral_angle[a];
  }

  delete [] dihedral_angle;
  delete [] hinge_label;
  delete [] pointer;
  delete [] labels_array;
  delete [] neighbors;

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Here we're testing each hinge for independence. If locking a hinge does not 
//   create a rigid cluster, then it is part of set of hinges whose motion depends 
//   on eachother. This group of dependent hinges is called a collective motion. 
//   Because each hinge has five bars covered by pebbles, and we know it is always 
//   possible to find at least 6 pebbles on a single site, we need only to find 
//   one additional pebble on the second site in order to lock the bonds.
//   The previously used function "min_hinge_label" has been inline coded here as
//   a while loop. 
// Parameters:
//   site1 - first atom of a bond.
//   site2 - second atom of a bond.
//   bars - total number of bars placed between site1 and site2
////////////////////////////////////////////////////////////////////////////////
int PebbleGame::lockHinge( const unsigned int site1, const unsigned int site2, const int bars ){

  unsigned int current_site = 0;
  int new_label = 0;
  int hinge_ID;
  int smallest_label = 2000000000;
  int pebbles_needed = 6 - bars;

  collectGroupMaxPebbles( site1 );

  // If we can find free pebbles on site2, use them to cover the 
  // necessary bars. Start with pebble6, as the free pebbles are most 
  // likely to be at the end of the pebble list. 
  //////////////////////////////////////////////////////////////////////
  int current_pebble = PEBBLE_6;
  do{
    if( pebbles[site2][current_pebble] == FREE_PEBBLE ){
      pebbles[site2][current_pebble] = site1;
      pebbles_needed--; 
    }
    current_pebble--;
  } while( pebbles_needed && current_pebble >= PEBBLE_1 );
  
  if( !pebbles_needed ) return( bars-6 );

  // A bond becomes locked when there are six covered bars between them.
  // Here we're going to try and find "6 - bars" pebbles in order to 
  // see if the bond between site1 and site2 is independent or not. If 
  // it is dependent, all the bonds searched in the "find_pebble" 
  // routine are in the same collective motion.
  //////////////////////////////////////////////////////////////////////
  while( pebbles_needed ){
    
    blocked_cutoff++;
    blocked_label[site1] = blocked_cutoff;
    blocked_label[site2] = blocked_cutoff;
    back_track[site2] = -1;

    findPebble( site2 );

    // If the pebble search from site2 did not fail, then there should be a 
    // free pebble on site2 now. Use this free pebble to cover the bond 
    // between site1 and site2. 
    //////////////////////////////////////////////////////////////////////
    if( !failed_pebble_search ){
      for( current_pebble = PEBBLE_1; current_pebble <= PEBBLE_6; current_pebble++ ){
	if( pebbles[site2][current_pebble] == FREE_PEBBLE )
	  pebbles[site2][current_pebble] = site1;
      }
      pebbles_needed--;
    } 
    else{
      // Execute this code if we fail the pebble search
      //////////////////////////////////////////////////////////////////////
      failed_pebble_search = false;
      //int next_site = 0;
      int a = 0;

      sites_to_check[0] = site1;
      sites_to_check[1] = site2;
      smallest_label = 200000000;
      int b;
      for( a = 0; a < total_sites_found; a++ ){
	current_site = sites_to_check[a];
	b = -1;

	neighbor_site = (site_info[current_site].neighbor_list).begin();
	while( neighbor_site != (site_info[current_site].neighbor_list).end() ){
 	  //for( b = 0; b < (site_info[current_site].neighbor_list).size(); b++ ){
	  //next_site = neighbors[pointer[current_site] +b];
	  
	  b++;

	  //hinge_ID = site_info[current_site].cm_label[*neighbor_site];
	  //hinge_ID = site_info[current_site].cm_label[next_site];

	  hinge_ID = labels_array[pointer[current_site] +b];
	    
	  // all rigid bonds have labels_array = hinge_ID = 0, don't check these.
	  // We want to find the smallest hinge_ID of any hinge in the failed pebble search 
	  // graph, and set smallest_label = the smallest hinge_ID. This will be used to label 
	  // all the hinges in the failed pebble search graph, ie. it will be the 
	  // collective motion label. 
	  ////////////////////////////////////////////////////////////////////////////////
	  if( *neighbor_site > current_site && 
	      hinge_ID > 0 &&
	      blocked_label[*neighbor_site] == blocked_cutoff ){
	    
	    new_label = hinge_label[hinge_ID];
	    
	    if( new_label < 0 ){
	      if( hinge_ID < smallest_label )
		smallest_label = hinge_ID;
	    }
	    else{
	      while( hinge_label[new_label] != new_label ){
		new_label = hinge_label[new_label];
	      }
	      if( new_label < smallest_label )
		smallest_label = new_label;
	    }
	  }
	  
	  neighbor_site++;
	}
      }
      
      hinge_label[smallest_label] = smallest_label;
      
      for( a = 0; a < total_sites_found; a++ ){
	current_site = sites_to_check[a];

	b = -1;
	neighbor_site = (site_info[current_site].neighbor_list).begin();
	while( neighbor_site != (site_info[current_site].neighbor_list).end() ){
	  //for( int b = 0; b < (site_info[current_site].neighbor_list).size(); b++ ){
	  //next_site = neighbors[pointer[current_site] +b];

	  b++;

	  //hinge_ID = site_info[current_site].cm_label[*neighbor_site];
	  hinge_ID = labels_array[pointer[current_site] +b];
	  
	  if( *neighbor_site > current_site && 
	      hinge_ID > 0 &&	    
	      blocked_label[*neighbor_site] == blocked_cutoff ){
	      
	    new_label = hinge_label[hinge_ID];
	    
	    if( new_label < 0 )
	      hinge_label[hinge_ID] = smallest_label;
	    else{
	      while( hinge_label[new_label] != new_label ){
		new_label = hinge_label[new_label];
	      }
	      hinge_label[new_label] = smallest_label;
	    }
	  }
	  
	  neighbor_site++;
	}
      }
      
      // we're returning the number of pebbles needed that actually 
      // covered bars before we failed the pebble search.
      // The number is 6 -bars -pebbles_needed.
      //////////////////////////////////////////////////////////////////////
      return( bars +pebbles_needed -6 );
    }

  }
    
  return( bars-6 );
}
////////////////////////////////////////////////////////////////////////////////
