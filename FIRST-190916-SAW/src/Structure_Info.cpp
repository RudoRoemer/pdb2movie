#include "global_defs.h"
#include "Structure_Info.h"
#include "hybrid_36_c.h"
#include "generalUtils.h"

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Program warnings are given when an operation may affect the results of the
//   program, but the program will still run. Three warning levels have been 
//   arbitrarily assigned, with warning level 1 being the most severe warngin. 
////////////////////////////////////////////////////////////////////////////////
void Structure_Info::add_warning( int warning_level, string warning ){
  
  if( warning_level == 1 )
    warning_level_1.push_back( warning );
  else if( warning_level == 2 )
    warning_level_2.push_back( warning );
  else if( warning_level == 3 )
    warning_level_3.push_back( warning );
  
  else{
    cout << " Error: Incorrect Warning Level used." << endl;
    exit(1);
  }
  
  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   A macro to populate the "site_info[]" array.
////////////////////////////////////////////////////////////////////////////////
void Structure_Info::add_to_site_info_array( SiteID site_1, SiteID site_2, int bars ){

  // Warning check:
  // Make sure these two atom numbers actually exist in the system
  //////////////////////////////////////////////////////////////////////
  
    
  site_info[site_1].neighbor_list.push_back( site_2 );
  site_info[site_1].number_of_bars.insert( pair<SiteID,int> (site_2, bars) );
  site_info[site_2].neighbor_list.push_back( site_1 );
  site_info[site_2].number_of_bars.insert( pair<SiteID,int> (site_1, bars) );
  
  site_info[site_1].cm_label.insert( pair<SiteID,int> (site_2, 0));
  site_info[site_2].cm_label.insert( pair<SiteID,int> (site_1, 0));
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void Structure_Info::remove_from_site_info_array( SiteID site_1, SiteID site_2 ){
  
  // 1. Erase site_2 from site_1's neighbor list.
  //////////////////////////////////////////////////////////////////////
  vector<SiteID>::iterator found_site_2 = (site_info[site_1].neighbor_list).begin();
  while( found_site_2 != (site_info[site_1].neighbor_list).end() ){
    if( *found_site_2 == site_2 ){
      (site_info[site_1].neighbor_list).erase( found_site_2 );
    }
    else
      found_site_2++;    
  }
  
  // 2. Erase the bond to site_2 from site_1's number of bars mapping.
  //////////////////////////////////////////////////////////////////////
  map<int,int>::iterator found_bond_2 = (site_info[site_1].number_of_bars).find(site_2);
  if( found_bond_2 != (site_info[site_1].number_of_bars).end() ){
    (site_info[site_1].number_of_bars).erase( found_bond_2 );
  }
  else{
    cout << "Error: couldn't find " << site_2 << " in " << site_1 << " map" << endl;
    exit(1);
  }
  
  //3. Erase site_1 from site_2's neighbor list.
  //////////////////////////////////////////////////////////////////////
  vector<SiteID>::iterator found_site_1 = (site_info[site_2].neighbor_list).begin();
  while( found_site_1 != (site_info[site_2].neighbor_list).end() ){
    if( *found_site_1 == site_1 ){
      (site_info[site_2].neighbor_list).erase( found_site_1 );
    }
    else
      found_site_1++;
  }
  
  // 4. Erase the bond to site_2 from site_1's number of bars mapping.
  //////////////////////////////////////////////////////////////////////
  map<int,int>::iterator found_bond_1 = (site_info[site_2].number_of_bars).find(site_1);
  if( found_bond_1 != (site_info[site_2].number_of_bars).end() ){
    (site_info[site_2].number_of_bars).erase( found_bond_1 );
  }
  else{
    cout << "Error: couldn't find " << site_2 << " in " << site_1 << " map" << endl;
    exit(1);
  }
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
void Structure_Info::reset_site_info_labels(){
  
  for( SiteID siteNumber = 0; siteNumber <= total_sites; siteNumber++ ){
    site_info[siteNumber].rigid_label = siteNumber;
    site_info[siteNumber].stressed_label = siteNumber;
  }
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void Structure_Info::reset_pruned_label(){
  
  for( unsigned int siteNumber = 0; siteNumber <= total_sites; siteNumber++ )
    site_info[siteNumber].pruned = false;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
string Structure_Info::get_orig_atom_number( SiteID atom_1 ){
  string str_atom_number;
  char write_atom_number[6];
  int  int_atom_number;

  str_atom_number= site_info[atom_1].orig_atom_number;
  if( isNumber( str_atom_number ) )
    {
      int_atom_number = atoi(  str_atom_number.c_str() );
      if(  int_atom_number > 99999){
	std::stringstream ss;
	hy36encode(5, int_atom_number, write_atom_number );
	ss << write_atom_number;
	str_atom_number = ss.str();
      }
    }

  return( str_atom_number );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void Structure_Info::freeGrid(){

  //size_t total = 0;

  for( unsigned int a = 0; a < coordinate_grid.size(); a++ ){
    //for( unsigned int b = 0; b < coordinate_grid[a].size(); b++ ){
    //cout << "Grid " << a << ":" << b << " = " << coordinate_grid[a].at(b) << " (" << sizeof(coordinate_grid[a].at(b)) << ")" <<endl;
    //total += sizeof(coordinate_grid[a].at(b) );
    //}
    //total += sizeof( coordinate_grid[a] );
    coordinate_grid[a].clear();
  }

  coordinate_grid.clear();
  //total += sizeof( coordinate_grid );
  //cout << "Current Grid Size: " << total << endl;
}
