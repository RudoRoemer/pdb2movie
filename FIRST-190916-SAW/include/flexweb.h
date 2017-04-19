#ifndef _FLEXWEB_
#define _FLEXWEB_

#include "global_defs.h"
#include "MolFramework.h"

void outputStatus();

void output_Jmol_script(ofstream &jmol_script, 
                        MolFramework &structure );
void output_Jmol_script( MolFramework &structure );
void output_FRODA_TIMME_Jmol_script( MolFramework &structure );
void output_RCD_XML( MolFramework &structure );
void output_PID( string path );
void output_details_flexweb();
void setTaskNames();

#endif
