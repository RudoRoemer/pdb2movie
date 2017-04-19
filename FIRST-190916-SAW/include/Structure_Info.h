#ifndef _STRUCTURE_INFO_
#define _STRUCTURE_INFO_

#include "global_defs.h"
#include "Site_Info.h"

class Structure_Info{

 public:

  void add_warning( int warning_level, string warning );
  void add_to_site_info_array( unsigned int site_1, unsigned int site_2, int bars );
  void freeGrid();
  void remove_from_site_info_array( unsigned int site_1, unsigned int site_2 );
  void reset_site_info_labels();
  void reset_pruned_label();

  string get_orig_atom_number( unsigned int atom_1 );

 public:

  bool using_pdb;
  bool using_ctf;
  bool using_mol2;
  bool using_mmcif;

  bool includeHydrogenbonds;
  bool includeHydrophobicTethers;
  
  bool using_dataset;
  bool using_bbg;

  bool using_axyz;
  bool unrecognized_file_format;  

  bool using_harlem;
  bool is_nmr;
  short int model_number;

  string infile_name;
  string base_name;
  string path_name;
  string ctf_file_name;
  map<int,string> harlem_names;

  int x_cells;
  int y_cells;
  int z_cells;

  float min_coords[3];
  float max_coords[3];
  vector< vector<unsigned int> > coordinate_grid;

  unsigned int total_residues;
  unsigned int floppy_modes;
  unsigned int total_sites;
  unsigned int total_excluded_sites;
  unsigned int total_bonds;
  unsigned int total_clusters;
  unsigned int total_stressed_regions;
  unsigned int total_clusters_larger_than_one;
  unsigned int total_clusters_larger_than_min_cluster_size;
  unsigned int isolated_sites;
  unsigned int isolated_dimers;
  
  float mean_coordination; 
  float pruned_mean_coordination;

  Site_Info *site_info;

  vector<string> warning_level_1;
  vector<string> warning_level_2;
  vector<string> warning_level_3;

  vector<new_bonds> dilution_list;
  vector< vector<int> > RC_atom_list;
  vector< vector<int> > SR_atom_list;

  list<unsigned int> chainTermini;

 public: 

  Structure_Info(){

    using_pdb = false;
    using_ctf = false;
    using_mol2 = false;
    using_mmcif = false;
    using_dataset = false;
    using_bbg = false;
    using_axyz = false;
    using_harlem = false;
    is_nmr = false;

    includeHydrogenbonds = true;
    includeHydrophobicTethers = true;

    model_number = 0;
    total_residues = 0;
    total_sites = 0;
    total_excluded_sites = 0;
    total_bonds = 0;
    floppy_modes = 0;
    total_clusters = 0;
    total_stressed_regions = 0;
    total_clusters_larger_than_one = 0;
    total_clusters_larger_than_min_cluster_size = 0;
    isolated_sites = 0;
    isolated_dimers = 0;
    x_cells = 0;
    y_cells = 0;
    z_cells = 0;

    min_coords[0] = 1000.0;
    min_coords[1] = 1000.0;
    min_coords[2] = 1000.0;
    max_coords[0] = -1000.0;
    max_coords[1] = -1000.0;
    max_coords[2] = -1000.0;
    mean_coordination = 0.0;
    pruned_mean_coordination = 0.0;
  }

  ~Structure_Info(){

  }

};
#endif
