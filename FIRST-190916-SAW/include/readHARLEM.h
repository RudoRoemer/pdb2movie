#ifndef _READ_HARLEM_
#define _READ_HARLEM_

#include "global_defs.h"
#include "MolFramework.h"

void read_harlem_file( MolFramework &structure );
int  scan_harlem_file( MolFramework &structure );
void read_harlem_data( MolFramework &structure );
string get_harlem_element_name( string &atom_name );

#endif
