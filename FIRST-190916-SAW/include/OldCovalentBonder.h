#ifndef OLDCOVALENTBONDER_H_
#define OLDCOVALENTBONDER_H_

#include "MolFramework.h"
#include "Site_Info.h"
#include "global_defs.h"

using namespace std;

class MolFramework;
class Site_Info;

class OldCovalentBonder
{
public:
  MolFramework* molFramework;
  Site_Info* site_info;
  vector<template_data> standard_template;
  vector<string> unknown_residue;
  map<string, float> covalent_bond_table;
  map<string, list<string> >::iterator atom_connections;
  vector<unsigned int>::iterator neighbor_atom;
  vector<int>::iterator current_grid;
  list<string>::iterator known_neighbor; 
  map<string, float>::iterator distance_cutoff;
    
  OldCovalentBonder(MolFramework* molFramework);
  ~OldCovalentBonder();
  
  void autobond();
  void identifyCovalentBonds();
  bool find_bonds( int current_template, unsigned int current_atom );
  bool find_bonds( unsigned int current_atom );
  bool is_disulfide_bond( unsigned int atom_1, unsigned int atom_2 );
  bool is_internucleic_bond( unsigned int atom_1, unsigned int atom_2 );
  bool is_known_neighbor( unsigned int current_atom, unsigned int neighbor, 
             list<string> list_of_neighbors, int *bars,
             int current_template );
  bool is_locked( int site_1, int mult_site_1,
             int site_2, int mult_site_2 );
  bool is_peptide_bond( unsigned int atom_1, unsigned int atom_2 );
  int meets_atom_atom_distance_cutoff( unsigned int atom_1, unsigned int atom_2 );
  void read_distance_cutoff_table();
  void read_standard_residue_library();
  bool search_het_group_dictionary( string residue_name );
  void store_residue_template_data( fstream *file_ptr, string residue_name );
};

#endif /*OLDCOVALENTBONDER_H_*/
