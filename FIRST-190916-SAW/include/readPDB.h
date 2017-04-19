#ifndef _READ_PDB_
#define _READ_PDB_

#include "global_defs.h"
#include "MolFramework.h"

void readPDB_File( MolFramework &structure );

int  scanPDB_File( MolFramework &structure );

int  getModelNumber( string infile_name );

int  listModels( string infile_name );

bool isCorrectModel( string current_line, int model_number );

void readPDB_Data( MolFramework &structure );

void readPDB_Data( MolFramework &structure, int model_number );

void readPDB_Data( vector<MolFramework*> &structures, 
		   vector<string> pdb_file_names);

void readPDB_Data( vector<MolFramework*> &structures, 
		   string pdb_file_name);

void readATOM_OrHETATM_Line( string &current_line, 
			     MolFramework &structure, 
			     int counter, 
			     int &rescount);

void readCONECT_Record( const string &current_line, 
			MolFramework &structure );

#endif
