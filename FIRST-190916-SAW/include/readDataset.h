#ifndef _READ_DATASET_
#define _READ_DATASET_

#include "global_defs.h"
#include "MolFramework.h"

void readDatasetFile( MolFramework &structure );
void readCF_AndTF_Constraints( MolFramework &structure );
void readHB_AndPH_Constraints( MolFramework &structure );

#endif
