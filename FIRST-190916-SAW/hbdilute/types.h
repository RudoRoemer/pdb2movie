/* Structure declarations for hbdilute.c */

#ifndef _TYPES_H_
#define _TYPES_H_

typedef struct chem_file{
  int  atom_or_hetatm;
  int  atom_number;
  char atom_type[4];
  int  residue;
  char chain_ID;
  int  insertion_space;
  struct chem_file *next_atom;
  struct chem_file *last_atom;
} atom_list;

typedef struct rcl{
  int  start_num;
  char start_type[4];
  int  end_num;
  char end_type[4];
  int  set_number;
  int  subset_number;
} rclist[30];

typedef struct n_c {
  int   atom_number;
  int   residue_number;
  int   insertion_space;
  char  atom_type[4];
  char  chain_ID;
  int   cluster_color;
  struct n_c *next_element;
} clusters;

#endif
