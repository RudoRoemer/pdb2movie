#include "generalUtils.h"
#include "PebbleGame.h" 

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Stressed region labels were assigned in the biuld_network routine 
//   (specifically in the cover_bonds routine). The first loop is necessary, 
//   because when a new new stressed region is found to overlap with a 
//   previously found stressed region (again, in the cover_bonds routine), 
//   only the smallest site in the previous stressed region is relabeled. 
//   Since all the other sites in the previous stressed region point to the 
//   label of the smallest site, we need only check two deep in the labeling.
////////////////////////////////////////////////////////////////////////////////
void PebbleGame::stressedRegionDecomposition(){

  unsigned int current_label = 0;

  blocked_cutoff++;

  // Relabel each site to have the same label as that of the site with the smallest
  // site number in the stressed region. Ensures proper labeling for stressed regions
  // that were found to be overlapping in the "build_network" routine. 
  ////////////////////////////////////////////////////////////////////////////////
  for( unsigned int a = 1; a <= total_sites; a++ )
    site_info[a].stressed_label = site_info[ site_info[a].stressed_label ].stressed_label;
  
  // Count the number of unique stressed regions and keep track of sites in 
  // stressed areas using the "blocked_label" array. Later, all sites that were
  // not blocked, will be labeled 0. 
  ////////////////////////////////////////////////////////////////////////////////
  for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){
    
    current_label = site_info[current_site].stressed_label;
    if( current_label != current_site ){
      if( blocked_label[current_label] < blocked_cutoff ){
	total_stressed_regions++;
	blocked_label[current_label] = blocked_cutoff;
      }
    }
  }

  // Label non-stressed region sites as 0, and output the stressed region label
  // for all sites. 
  ////////////////////////////////////////////////////////////////////////////////
  for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){
    if( blocked_label[ site_info[current_site].stressed_label ] < blocked_cutoff )
      site_info[current_site].stressed_label = 0;
  }

  // BMH debug code. Replaced with *_data.txt file. 2/13/05
  // Output the stressed region label for all sites. 
  ////////////////////////////////////////////////////////////////////////////////
  //ofstream stressed_labels("srd.txt");
  //for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){
  //stressed_labels << setw(8) << current_site 
  //  << setw(8) << site_info[current_site].stressed_label << endl;
  //}
  //stressed_labels.close();

  // Relabel the stressed regions in ascending order according to size. The
  // the largest region will be labeled 1. The "new_labels" vector is accessed
  // by rigid cluster label. Will push every atom onto the vector, and then sort
  // the vector based on size. 
  ////////////////////////////////////////////////////////////////////////////////
  vector<vector<unsigned int> > new_labels(total_sites+1);

  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    if( site_info[siteNumber].stressed_label > 0 )
      new_labels[site_info[siteNumber].stressed_label].push_back(siteNumber);
  }

  stable_sort( new_labels.begin(), new_labels.end(), sortOnSize );

  //for( unsigned int a = 0; a <= total_sites; a++ )
  //cout << a << " " << new_labels[a].size() << endl;

  for( unsigned int a = 0; a <= total_stressed_regions; a++ ){
    for( unsigned int b = 0; b < new_labels[a].size(); b++ ){
      site_info[ new_labels[a].at(b) ].stressed_label = a+1;
    }
  }

  // Error Check
  //////////////////////////////////////////////////////////////////////
  for( unsigned int a = 1; a <= total_sites; a++ ){
    if( site_info[a].stressed_label > total_stressed_regions ){
      cout << " Error in labeling stressed regions. " 
	   << a << " " << site_info[a].stressed_label << " > " << total_stressed_regions << endl;
      exit(0);
    }
  }
  
  // Populate the SR_atom_list vector. This vector is indexed by stressed region
  // label. It will return a vector<int> containing all the atoms that belong
  // to the stressed region that was indexed.
  ////////////////////////////////////////////////////////////////////////////////
  molFramework->SR_atom_list.resize( total_stressed_regions+1 );
  for( int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    molFramework->SR_atom_list[ site_info[siteNumber].stressed_label ].push_back(siteNumber);
  }

  // Compute the flexibility index for each stressed region. The default
  // value for flexibility index is zero, so rigid-isostatic regions do 
  // not need to be evaluated, provided we handle flexible and stressed 
  // regions, which we do. 
  //////////////////////////////////////////////////////////////////////
  float isostatic_value = 0.0;
  float flex_index = 0.0;
  int bars_in_fully_connected_graph = 0;
  int SR_bars[total_stressed_regions +1];
  int SR_size = 0;
  int current_site = 0;
  int neighbor_site = 0;


  // compute the flexibility index for each region. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int current_stressed_region = 1; current_stressed_region <= total_stressed_regions; current_stressed_region++ ){

    SR_bars[current_stressed_region] = 0;
    
    for( unsigned int a = 0; a < molFramework->SR_atom_list[current_stressed_region].size(); a++ ){
      current_site = molFramework->SR_atom_list[current_stressed_region][a];

      for( unsigned int b = 0; b < site_info[current_site].neighbor_list.size(); b++ ){
	neighbor_site = site_info[current_site].neighbor_list[b]; 

	if( neighbor_site > current_site &&
	    site_info[current_site].stressed_label == site_info[neighbor_site].stressed_label ){
	  SR_bars[current_stressed_region] += site_info[current_site].number_of_bars[neighbor_site];
	}
      }
    }
  }

  for( unsigned int current_stressed_region = 1; current_stressed_region <= total_stressed_regions; current_stressed_region++ ){
    
    SR_size = molFramework->SR_atom_list[current_stressed_region].size();
    isostatic_value = (6 * SR_size) -6;
    bars_in_fully_connected_graph = 3 * (SR_size) * (SR_size -1);
    flex_index = (SR_bars[current_stressed_region] -isostatic_value) / (bars_in_fully_connected_graph -isostatic_value);

    //cout << "Stressed region " << current_stressed_region << " " << SR_size << " " << SR_bars[current_stressed_region] << endl;
    //cout << "    " << isostatic_value << " " << bars_in_fully_connected_graph << " " << flex_index << endl;

    for( int a = 0; a < SR_size; a++ ){
      current_site = molFramework->SR_atom_list[current_stressed_region][a];

      for( unsigned int b = 0; b < site_info[current_site].neighbor_list.size(); b++ ){
	neighbor_site = site_info[current_site].neighbor_list[b]; 

	if( neighbor_site > current_site &&
	    site_info[current_site].stressed_label == site_info[neighbor_site].stressed_label ){
	  site_info[current_site].bond_info[neighbor_site].flexibility_index = -flex_index;
	  site_info[neighbor_site].bond_info[current_site].flexibility_index = -flex_index;
	}
      }
    }
  }

  // TODO add error check that all atoms are in SR_atom_list (non-stressed atoms 
  // are in the zeroth index.

  return;
}
////////////////////////////////////////////////////////////////////////////////
