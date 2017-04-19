#ifndef _HDFUNCTIONS_
#define _HDFUNCTIONS_

#include "types.h"

void read_chem_file( atom_list **, FILE *, int[], int[], char[], int *, int *, int[], int[], char[] );
int  compute_number_of_clusters( int [], atom_list *, int [], int );
void first_decomp( int, int, int[], atom_list *, clusters *[], int [], int [30] );
void set_new_decomp_info( int, int[], int, clusters *[], atom_list *, int [] );
void compute_new_cluster_list( clusters *[], int, clusters *[], int, int [], int [] );

int  compare_old_and_new_decomps( clusters *[120], clusters *[120], int, int[], int );

void set_temp_files( FILE *[], int *, int * );
void clean_up( FILE *, int);

#endif


