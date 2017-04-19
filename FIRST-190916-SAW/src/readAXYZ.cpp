#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "readAXYZ.h"

extern const Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Analyze the flexibility in a PDB file. This function is defined as virtual
//   in the base class "topology_object", and is redefined here. 
////////////////////////////////////////////////////////////////////////////////
void readAXYZ_File( MolFramework &structure ){

  ////////////////////////////////////////////////////////////////////////////////
  int total_sites = scanFile( structure );

  ////////////////////////////////////////////////////////////////////////////////
  structure.initialize( total_sites );

  ////////////////////////////////////////////////////////////////////////////////
  readAXYZ_Data( structure );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
int scanFile( MolFramework &structure ){

  int total_sites = 0;
  string linebuf;
  fstream axyz_file( structure.infile_name.c_str() );

  while( !axyz_file.eof() ){
    getline(axyz_file, linebuf);
    if( linebuf.find_first_of(" \n", 0 ) != string::npos )
      total_sites++;
  }
  axyz_file.close();

  return( total_sites );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void readAXYZ_Data( MolFramework &structure ){

  int field_start = 0;
  int field_end = 0;
  int counter = 0;
  string number;
  string linebuf;

  // 1. Open the file. 
  //////////////////////////////////////////////////////////////////////
  fstream infile( structure.infile_name.c_str() );
  if( !infile ){
    cout << "Error: Could not open file (" << structure.infile_name << "). " << endl;
    exit(1);
  }

  std::stringstream ss;

  while( !infile.eof() ){
    getline(infile, linebuf);
    
    if( linebuf.find_first_of(" \n", 0 ) != string::npos ){
      counter++;

      // Read the atom name
      //////////////////////////////////////////////////////////////////////
      field_start = linebuf.find_first_not_of(" \n", 0);
      field_end   = linebuf.find_first_of(" \n", field_start); 
      structure.site_info[counter].atom_name = linebuf.substr( field_start, field_end-field_start );
      structure.site_info[counter].element_name = structure.getElementName( structure.site_info[counter].atom_name );
      ss<< counter;
      structure.site_info[counter].orig_atom_number = ss.str();
      structure.site_info[counter].FIRST_number = counter;
      
      // Read the X-coord.
      //////////////////////////////////////////////////////////////////////
      field_start = linebuf.find_first_not_of(" \n", field_end);
      field_end   = linebuf.find_first_of(" \n", field_start);
      number = linebuf.substr( field_start, field_end-field_start );
      structure.site_info[counter].coords[X] = atof( number.c_str() );
      //  Compare to min/max in X direction. 
      if( structure.site_info[counter].coords[X] < structure.min_coords[X] )
	structure.min_coords[X] = structure.site_info[counter].coords[X];
      if( structure.site_info[counter].coords[X] > structure.max_coords[X] )
	structure.max_coords[X] = structure.site_info[counter].coords[X];
      
      // Read the Y-coord.
      //////////////////////////////////////////////////////////////////////
      field_start = linebuf.find_first_not_of(" \n", field_end);
      field_end   = linebuf.find_first_of(" \n", field_start);
      number = linebuf.substr( field_start, field_end-field_start );
      structure.site_info[counter].coords[Y] = atof( number.c_str() );
      //  Compare to min/max in Y direction. 
      if( structure.site_info[counter].coords[Y] < structure.min_coords[Y] )
	structure.min_coords[Y] = structure.site_info[counter].coords[Y];
      if( structure.site_info[counter].coords[Y] > structure.max_coords[Y] )
	structure.max_coords[Y] = structure.site_info[counter].coords[Y];
      
      // Read the Z-coord.
      //////////////////////////////////////////////////////////////////////
      field_start = linebuf.find_first_not_of(" \n", field_end);
      field_end   = linebuf.find_first_of(" \n", field_start);
      number = linebuf.substr( field_start, field_end-field_start );
      structure.site_info[counter].coords[Z] = atof( number.c_str() );
      //  Compare to min/max in Z direction. 
      if( structure.site_info[counter].coords[Z] < structure.min_coords[Z] )
	structure.min_coords[Z] = structure.site_info[counter].coords[Z];
      if( structure.site_info[counter].coords[Z] > structure.max_coords[Z] )
	structure.max_coords[Z] = structure.site_info[counter].coords[Z];
      
    }
  }

  infile.close();
}
////////////////////////////////////////////////////////////////////////////////

