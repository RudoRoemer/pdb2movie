#ifndef _SITE_INFO_
#define _SITE_INFO_

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "SiteID.h"

using namespace std;

class MolFramework;

struct bond_data{

  float flexibility_index;

  bond_data(){
    flexibility_index = 0.0;
  }
  
};

////////////////////////////////////////////////////////////////////////////////
struct Site_Info{
  MolFramework *molFramework;
  
  vector<SiteID> neighbor_list;
  int nCovalentNeighbors;
  int nHydrogens;
  map<unsigned int,int> cm_label;
  map<int,int> number_of_bars;

  map<int,bond_data> bond_info;

  unsigned int rigid_label;
  unsigned int stressed_label; 
  int coll_mode_label;

  float flexibility;
  float flexibility_index;

  bool excluded;
  bool pruned;
  bool checked;
  bool surface_site;
  bool isExposed;

  string record_name;
  //KGS 03/04/08 Changed orig_atom_number from int to string
  //int orig_atom_number;
  string orig_atom_number;
  int FIRST_number;
  string atom_name;
  string element_name;
  string alt_location;
  string residue_name;
  string insert_code;
  unsigned int chain_ID;
  SiteID site_number;
  int hbond_status;

  int FIRST_chain_ID;
  int seq_number;

  int FIRST_group_ID; // FIRST-defined field, columns 82-87 of PDB file

  int grid_X;
  int grid_Y;
  int grid_Z;

  bool ambiguous_atom_type;
  bool no_missing_neighbors;

  float coords[3];
  float occupancy;
  float temp_factor;
  float vdw_radius;

  string seg_id;
  float  charge;

  Site_Info(){
    initialize();
  };
  
  Site_Info(MolFramework * molFramework){
    initialize();
    setMolFramework(molFramework);
  };
  
  void setMolFramework(MolFramework* molFramework) {
    this->molFramework = molFramework;
  };
  
  MolFramework* getMolFramework() {
    return molFramework;
  };
  
  void initialize() {
    molFramework = NULL;
    excluded = pruned = checked = false;
    surface_site = false;
    isExposed = false;
    rigid_label = 0;
    stressed_label = 0;
    coll_mode_label = 0;
    //KGS 03/04/08 Changed orig_atom_number from int to string
    //orig_atom_number = 0;
    flexibility = 0.0;
    flexibility_index = 1.0;
    FIRST_number = 0;
    chain_ID = 0;
    FIRST_chain_ID = 0;
    seq_number = 0;
    FIRST_group_ID = 0;
    coords[0] = 0.0;
    coords[1] = 0.0;
    coords[2] = 0.0;
    grid_X = grid_Y = grid_Z = 0;
    hbond_status = 0;
    no_missing_neighbors = true;
    occupancy = 0.0;
    temp_factor = 0.0;
    vdw_radius = 0.0;
    ambiguous_atom_type = false;
    charge = 0.0;
    element_name = "  ";
    nCovalentNeighbors = 0;
    nHydrogens = 0;
  };

  ~Site_Info(){};

  void print(){
    cout << setw(10) << orig_atom_number 
	 << setw(10) << FIRST_number 
	 << setw(10) << atom_name
	 << setw(10) << element_name
	 << setw(2)  << chain_ID
	 << setw(2)  << FIRST_chain_ID
	 << setw(10) << seq_number
       	 << setw(15) << coords[0]
	 << setw(15) << coords[1]
	 << setw(15) << coords[2] 
	 << setw(5)  << excluded
	 << endl;
  }

};

#endif
