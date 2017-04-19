#include "global_defs.h"
#include "generalUtils.h"
#include "PebbleGame.h" 

extern const Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Rigid_cluster_decomposition - identifies sites belonging to the same rigid 
//   cluster, and assigns them a unique cluster label. If the number of neighbors 
//   of the "initial site" is > than 1, a search for sites belonging to the same 
//   rigid cluster as the "initial site" begins. The search consists of two loops. 
//   The outer search checks all of the topological neighbors of the "initial  
//   site". For each neighbor found, an inner loop commences which will perform 
//   a pebble search from that site. The inner loop pebble game assigns sites as 
//   either RIGID or FLOPPY. Each site that is found to be rigid, is then added 
//   to the topological search, which means a new pebble game search will start 
//   from every neighbor of every site that is found rigid (unless that neighbor 
//   has already been labeled RIGID or FLOPPY). The outer loop exits when no site 
//   that can be found on any branch from the "initial_site" can be found that 
//   hasn't already been labeled RIGID or FLOPPY.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::rigidClusterDecomposition(){ 

  unsigned int this_node_cluster_label = 0;
  int node = 0;

  // output the rigid cluster decomposition for debugging purposes. 
  ////////////////////////////////////////////////////////////////////////////////
  //ofstream rcd( "rcd.txt", ios::out );
  //rcd << "total sites " << total_sites << endl;
  //for( int d = 1; d <= total_sites; d++)
  //rcd << setw(8) << d << setw(8) << site_info[d].rigid_label << endl;
  //rcd.close();

  for( SiteID initial_site = 1; initial_site <= total_sites; initial_site++ ){
    
    if( site_info[initial_site].neighbor_list.size() > 1 ){
      
      if( site_info[initial_site].rigid_label == initial_site ){
	//cout << " initializing pebble search on site: " << initial_site << endl;
	initializeTopologySearchAtSite( initial_site );
	saveThisSiteToSearchList( initial_site );

	// Start searching the topological neighbors of "initial_site"
	////////////////////////////////////////////////////////////////////////////////
	while( sites_checked < total_sites_found ){
	  node = sites_to_check[sites_checked];

	  neighbor_site = (site_info[node].neighbor_list).begin();
	  while( neighbor_site != (site_info[node].neighbor_list).end() ){

	    // Start a new pebble search for each neighbor that hasn't been checked, and 
	    // hasn't been assigned to be rigid or flexible in a previous pebble search
	    ////////////////////////////////////////////////////////////////////////////////
	    if( blocked_label[*neighbor_site] != CHECKED ){
	      this_node_cluster_label = site_info[*neighbor_site].rigid_label;

	      if( blocked_label[this_node_cluster_label] == CHECKED ){
		site_info[*neighbor_site].rigid_label = initial_site;
		saveThisSiteToSearchList( *neighbor_site );
	      }
	      else{

		if( blocked_label[this_node_cluster_label] != is_floppy ){
		  if( site_info[this_node_cluster_label].neighbor_list.size() == 1 ||
		      this_node_cluster_label >= initial_site ){
		    //cout << "starting search from " << this_node_cluster_label << endl;
		    initializePebbleSearchAtSite( *neighbor_site );
		    startPebbleSearchAtSite( initial_site );
		  }
		}

	      }
	    }
	    neighbor_site++;
	  }
	  sites_checked++;
	}
	resetBlockedLabels();
      }
    }
    
    // If the current site has only one neighbor, see if it's in an isolated
    // dimer. 
    //////////////////////////////////////////////////////////////////////
    else if( (site_info[initial_site].neighbor_list).size() == 1 ){
      
      neighbor_site = (site_info[initial_site].neighbor_list).begin();
      if( (site_info[*neighbor_site].neighbor_list).size() == 1 && 
	  initial_site < *neighbor_site ){
	isolated_dimers++;
	site_info[*neighbor_site].rigid_label = site_info[initial_site].rigid_label;
      }
      else if( site_info[*neighbor_site].neighbor_list.size() > 1 ){
	if( site_info[initial_site].number_of_bars[*neighbor_site] >= 1 )
	  site_info[initial_site].rigid_label = site_info[*neighbor_site].rigid_label;
	else
	  total_clusters++;
      }
    }
    
    // The current site has no neighbors, it is an isolated site.
    //////////////////////////////////////////////////////////////////////
    else
      isolated_sites++;
  }
  
  total_clusters += (isolated_dimers + isolated_sites);

  // BMH hack code to generate artificially large rigid cluster
  /*
  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
  if( siteNumber >= 2985 && siteNumber <= 5108 ){
      site_info[siteNumber].rigid_label = 1;
    }
    else{
      site_info[siteNumber].rigid_label++;
    }
  }
  */
   
  // output the rigid cluster decomposition for debugging purposes. 
  ////////////////////////////////////////////////////////////////////////////////
  //ofstream rcd( "rcd.txt", ios::out );
  //rcd << "total sites " << total_sites << endl;
  //for( int d = 1; d <= total_sites; d++)
  //rcd << setw(8) << d << setw(8) << site_info[d].rigid_label << endl;
  //rcd.close();
  
  // Relabel the rigid cluster labels in ascending order according to size. The
  // the largest cluster will be labeled 1. The "new_labels" vector is accessed
  // by rigid cluster label. Will push every atom onto the vector, and then sort
  // the vector based on size. 
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<unsigned int> > new_labels(total_sites+1);
  unsigned int error_check_total_clusters = 0;

  //ofstream test("test.out");
  for(SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    new_labels[site_info[siteNumber].rigid_label].push_back(siteNumber);
    //test << setw(8) << siteNumber << setw(8) << site_info[siteNumber].rigid_label << endl;

    if( site_info[siteNumber].neighbor_list.size() == 1 &&
	site_info[siteNumber].rigid_label != site_info[site_info[siteNumber].neighbor_list[0] ].rigid_label &&
	!molFramework->using_bbg )
      cout << "found bug " << siteNumber << " (" << site_info[siteNumber].rigid_label << ") " 
	   << site_info[siteNumber].neighbor_list[0] << " (" << site_info[site_info[siteNumber].neighbor_list[0]].rigid_label << ")" << endl;
  }
  //test.close();

  stable_sort( new_labels.begin(), new_labels.end(), sortOnSize );
  
  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){

    
    if( new_labels[siteNumber-1].size() > 1 )
      total_clusters_larger_than_one++;

    if( new_labels[siteNumber-1].size() > parameters.min_output_cluster_size )
      total_clusters_larger_than_min_cluster_size++;



    for( unsigned int labelNumber = 0; labelNumber < new_labels[siteNumber-1].size(); labelNumber++ )
      site_info[ new_labels[siteNumber-1].at(labelNumber) ].rigid_label = siteNumber;

    if( new_labels[siteNumber-1].size() > 0 )
      error_check_total_clusters++;
  }

  //total_clusters = error_check_total_clusters;
  
  if( error_check_total_clusters != total_clusters ){
    cout << " Error: Total clusters found by labeling does not equal the total number of" << endl;
    cout << "        clusters found by the pebble game." << endl;
    cout << "        [" << error_check_total_clusters << "] != [" << total_clusters << "]" << endl;
    cout << "        Writing the rigid cluster decomposition to file rcd.txt." << endl;
    ofstream rcd( "rcd.txt", ios::out );
    for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++)
      rcd << setw(8) << siteNumber << setw(8) << site_info[siteNumber].rigid_label << endl;
    rcd.close();
    exit(1);
  }

  // Populate the RC_atom_list vector. This vector is indexed by rigid cluster
  // label. It will return a vector<int> containing all the atoms that belong
  // to the rigid cluster label that was indexed.
  ////////////////////////////////////////////////////////////////////////////////
  molFramework->RC_atom_list.clear();
  molFramework->RC_atom_list.resize( total_clusters+1 );
  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    (molFramework->RC_atom_list[ site_info[siteNumber].rigid_label ]).push_back(siteNumber);
  }


  

  // Tag sites as bulk (connected to only other atoms in the same rigid cluster) or
  // surface (connected to another atom via a flexible bond). Useful in output
  // routines and also for validating certain counts. 
  ////////////////////////////////////////////////////////////////////////////////
  SiteID neighbor = 0;
  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    for( SiteID neighborNumber = 0; neighborNumber < (site_info[siteNumber].neighbor_list).size(); neighborNumber++ ){
      neighbor = site_info[siteNumber].neighbor_list[neighborNumber];
      if( neighbor > siteNumber &&
	  site_info[siteNumber].rigid_label != site_info[neighbor].rigid_label ){
	site_info[siteNumber].surface_site = true;
	site_info[neighbor].surface_site = true;
	break;
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Initialize_new_rigid_cluster_search.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::initializeTopologySearchAtSite( int initial_site ){

      total_clusters++;
      is_floppy = -blocked_cutoff;
      collectGroupMaxPebbles( initial_site ); 

      total_sites_found = 0;
      sites_checked = 0;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Save_this_site_to_search_list.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::saveThisSiteToSearchList( int this_site ){ 

  sites_to_check[total_sites_found] = this_site;
  blocked_label[this_site] = CHECKED;

  total_sites_found++;

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Initialize_pebble_search
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::initializePebbleSearchAtSite( int this_node ){

  blocked_cutoff++;
  pebble_sites_searched = total_sites_found;
  pebbles_followed = total_sites_found;

  blocked_label[this_node] = blocked_cutoff;
  back_track[this_node] = START_NODE;
  sites_to_check[pebble_sites_searched] = this_node;

  pebble_sites_searched++;

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void PebbleGame::startPebbleSearchAtSite( int initial_site ){

  int 
    this_pebble = 0,
    this_pebbles_cluster_label = 0,
    current_node = 0;    

  while( pebbles_followed < pebble_sites_searched ){
   
    current_node = sites_to_check[pebbles_followed];
    
    // Start pebble searching ...
    ////////////////////////////////////////////////////////////////////////////////
    for( int current_pebble = PEBBLE_1; current_pebble <= PEBBLE_6; current_pebble++ ){
      
      this_pebble = pebbles[current_node][current_pebble];

      if( this_pebble == FREE_PEBBLE ){
	while( current_node != START_NODE ){
	  blocked_label[ site_info[current_node].rigid_label ] = is_floppy;
	  current_node = back_track[current_node];
	}
	
	//cout << "found free pebble; returning. " << endl;
	//printPebbleCovering();
	return;
      }
      else if( this_pebble != EMPTY_PEBBLE_SPACE ){
	if( blocked_label[this_pebble] < blocked_cutoff &&
	    blocked_label[this_pebble] != CHECKED ){
	  
	  this_pebbles_cluster_label = site_info[this_pebble].rigid_label;
	  
	  if( this_pebbles_cluster_label < initial_site ||
	      blocked_label[this_pebbles_cluster_label] == is_floppy ){
	    while( current_node != START_NODE ){
	      blocked_label[ site_info[current_node].rigid_label ] = is_floppy;
	      current_node = back_track[current_node];
	    }      
	    return;
	  }
	  else{
	    sites_to_check[pebble_sites_searched] = this_pebble;
	    pebble_sites_searched++;
	    blocked_label[this_pebble] = blocked_cutoff;
	    back_track[this_pebble] = current_node;
	  }
	}
      }
      

    }
    // Finished pebble searching
    pebbles_followed++;
  }
  
  // We failed a pebble search. 
  for( int b = total_sites_found; b < pebble_sites_searched; b++ ){
    site_info[ sites_to_check[b]  ].rigid_label = initial_site;
    blocked_label[ sites_to_check[b] ] = CHECKED;
  }
   
  total_sites_found = pebble_sites_searched;
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void PebbleGame::resetBlockedLabels(){

  for( int a = 0; a < total_sites_found; a++ )
    blocked_label[ sites_to_check[a] ] = 0;

}
////////////////////////////////////////////////////////////////////////////////
