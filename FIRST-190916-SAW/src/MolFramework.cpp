#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "MolFramework.h"
#include "CovalentBonder.h"
#include "OldCovalentBonder.h"
#include "SiteID.h"
#include "hybrid_36_c.h"

#include <cstring>

extern const Parameters parameters;

unsigned int new_bonds::total_in_network;

const Vector NULL_VEC = Vector(0,0,0);

inline float rad2deg( float a ){
  return(a*(180.0/PI));
}
inline float deg2rad( float a ){
  return( a*(PI/180.0) );
}
////////////////////////////////////////////////////////////////////////////////


MolFramework::~MolFramework(){
  //delete [] site_info; FIXME - figure out why this segfaults for Adam's cell data in TIMME
  
  orig_2_FIRST_atom_number.clear();
};


////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Read a list of user-defined constraints.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::addUserDefinedConstraints(){
  
  if( !parameters.read_user_defined_constraint_file )
    return;

  if( parameters.verbose > 0 )
    cout << "Reading user-defined constraints from file..." << endl;
  
  string filename = "user_def.list";
  string linebuf;
  ifstream user_def_list;

  // If we are in interactive mode, prompt the user for the name of 
  // the file listing their hydrogen bonds. If we're in non-interactive
  // mode, look for a file named hbond.list. If this file does not 
  // exist, exit. 
  //////////////////////////////////////////////////////////////////////
  user_def_list.open( filename.c_str() );
  if( !user_def_list ){
    clear_screen;
    cout << endl;
    cout << " The file [" << filename << "] was not found in this directory." << endl << endl;
    do{
      user_def_list.clear();
      cout << " Please enter the name of the file that contains the list of user-defined constraints." << endl;
      cout << " Filename = ";
      getline(cin, filename);
      
      user_def_list.open( filename.c_str(), ios::in );
      if( !user_def_list || did_press_enter(filename) ){
	clear_screen;
	cout << endl << "File not found." << endl << endl;
      }
      
    } while( !user_def_list || did_press_enter(filename) );
  }
  
  // Read each line of the Hphobic list file. The first number should
  // correspond to a hydrogen atom, and the second atom to the accpt
  // atom. An optional third number can be used to force a user-defined
  // energy for the given bond. If no energy is given, it will be 
  // computed. 
  //////////////////////////////////////////////////////////////////////
  unsigned int atom_1 = 0;
  unsigned int atom_2 = 0;
  int bars   = parameters.hp_bars;
  size_t field_start = 0;
  size_t field_end = 0;
  string number;

  while( !user_def_list.eof() ){
    getline(user_def_list, linebuf);

    // Skip blank lines.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos )
      break;

    // Read in site 1 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    atom_1 = atoi( number.c_str() );
    
    // Read in site 2 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    atom_2 = atoi( number.c_str() );
    
    // OPTIONAL Read the number of bars to use in the constraint, if it is listed. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      if( number.find("123456789") == string::npos)
	bars = parameters.ud_bars;
      else
	bars = atoi( number.c_str() );

      if( bars < 0 ||
	  bars > 6 ){
	cout << " ERROR: Found a user-defined constraint with an invalid number of bars while reading" << endl;
	cout << "        the file " << filename << "." << endl;
	cout << "        INVALID LINE: " << linebuf << endl;
	exit(1);
      }
    }

		if (!parameters.use_first_numbering) {
			// Map the original numbers into the numbers used internally by FIRST.
			//////////////////////////////////////////////////////////////////////    
			atom_1 = orig_2_FIRST_atom_number[atom_1];
			atom_2 = orig_2_FIRST_atom_number[atom_2];
		} 

    if( atom_1 <= 0 || atom_1 > total_sites ||
	atom_2 <= 0 || atom_2 > total_sites ){
      cout << " ERROR: An atom in the user-defined constraint file [" << filename << "] does not exist." << endl;
      cout << "        The error occurred with the following line:" << endl;
      cout << "        " << linebuf << endl << endl;
      exit(1);
    }

    //cout << atom_1 << " " << atom_2 << " " << bars << endl;
    user_defined_constraints.push_back( new_bonds(atom_1, atom_2, bars) );
  }

  if( parameters.verbose )
    cout << "Read " << user_defined_constraints.size() << " user-defined constraints from external file." << endl;
  user_def_list.close();

}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
// Description:
//   Function to map the X-, Y-, and Z-grid positions of an atom to the array
//   position of the 1-dimensional "coordinate_grid" array. 
// Parameters:
//   grid_X - The x-coordinate of the 3D grid that contains the atom coordinates.
//   grid_Y - The y-coordinate of the 3D grid that contains the atom coordinates.
//   grid_Z - The z-coordinate of the 3D grid that contains the atom coordinates.
// Return Value List:
//   The position in the 1-dimensional array that represents the coordinat grid.
////////////////////////////////////////////////////////////////////////////
int MolFramework::gridPosition( int grid_X, int grid_Y, int grid_Z ){

  return( ( grid_Z * x_cells * y_cells) + 
	  ( grid_Y * x_cells ) +
	  ( grid_X ) );
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   The 3D coordinate grid is stored in a one-dimensional vector of integers. This 
//   function maps the the 3D position to the 1D postion.
// Parameters:
//   current_atom - The FIRST number of an atom.
// Return Value List:
//   int: The index of the array "coordinate_grid" to which current_atom belongs.
////////////////////////////////////////////////////////////////////////////
int MolFramework::gridPosition( int current_atom ){

  return( ( site_info[current_atom].grid_Z * x_cells * y_cells) + 
	  ( site_info[current_atom].grid_Y * x_cells ) +
	  ( site_info[current_atom].grid_X ) );
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Create a list of the 27 coordinate grids that form a 3-dimensional cube surrounding
//   the current atom's grid, which is in the center of the cube.
// Parameters:
//   current_atom - FIRST_number of an atom in the molecule. 
// Return Value List:
//   A vector of 27 integers that correspond to the array postions in "coordinat_grid"
//   that will be searched when finding neighbors of the current_atom.
////////////////////////////////////////////////////////////////////////////
//inline
vector<int> MolFramework::getGrids( int current_atom ){

  int current_X = site_info[current_atom].grid_X;
  int current_Y = site_info[current_atom].grid_Y;
  int current_Z = site_info[current_atom].grid_Z;
  
  vector<int> grids_to_check;

  for( int a = -1; a <= 1; a++ ){
    for( int b = -1; b <= 1; b++ ){
      for( int c = -1; c <= 1; c++ ){
	grids_to_check.push_back( gridPosition(current_X+a, current_Y+b, current_Z+c) );
      }
    }
  }

  return( grids_to_check );
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Create a list of the 27 coordinate grids that form a 3-dimensional cube surrounding
//   the current atom's grid, which is in the center of the cube.
// Parameters:
//   current_atom - FIRST_number of an atom in the molecule. 
// Return Value List:
//   A vector of 27 integers that correspond to the array postions in "coordinat_grid"
//   that will be searched when finding neighbors of the current_atom.
////////////////////////////////////////////////////////////////////////////
void MolFramework::getGrids( int current_atom, int grids[27] ){

  int current_X = site_info[current_atom].grid_X;
  int current_Y = site_info[current_atom].grid_Y;
  int current_Z = site_info[current_atom].grid_Z;
  int counter = 0;

  //if( current_atom == 357 ){
  //cout << "Grid position of atom 357" << endl;
  //cout << "\t" << site_info[current_atom].grid_X << "\t" << site_info[current_atom].grid_Y << "\t" << site_info[current_atom].grid_Z << endl;
  //}

  for( int a = -1; a <= 1; a++ ){
    for( int b = -1; b <= 1; b++ ){
      for( int c = -1; c <= 1; c++ ){
	grids[counter] = gridPosition(current_X+a, current_Y+b, current_Z+c);
	counter++;
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Compute the Cartesian distance between two atoms.
// Parameters:
//   atom_1 - FIRST_number of a given atom.
//   atom_2 - FIRST_number of a given atom.
// Return Value List:
//   distance - The Cartesian distance between atom_1 and atom_2.
////////////////////////////////////////////////////////////////////////////
float MolFramework::getDistance( SiteID atom_1, SiteID atom_2 ){
	
  float X1 = site_info[atom_1].coords[X];
  float Y1 = site_info[atom_1].coords[Y];
  float Z1 = site_info[atom_1].coords[Z];

  float X2 = site_info[atom_2].coords[X];
  float Y2 = site_info[atom_2].coords[Y];
  float Z2 = site_info[atom_2].coords[Z];

  return sqrt( (X1-X2)*(X1-X2) + (Y1-Y2)*(Y1-Y2) + (Z1-Z2)*(Z1-Z2) );
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Determine the atom_1 - atom_2 -atom_3 bond angle using the law of 
//   cosines. The return value is in radians.
// Parameters:
//   atom_1 - point A of angle ABC. 
//   atom_2 - point B of angle ABC.
//   atom_3 - point C of angle ABC.
// Return Value:
//   acos(angle) - The angle in radians formed by atoms 1, 2, and 3. 
////////////////////////////////////////////////////////////////////////////
float MolFramework::computeAngle( SiteID atom_1, SiteID atom_2, SiteID atom_3 ){

  float dist_1  = getDistance( atom_1, atom_2 );
  float dist_2  = getDistance( atom_2, atom_3 );
  float dist_3  = getDistance( atom_1, atom_3 );

  float angle = ( -pow(dist_3,2) +pow(dist_1,2) +pow(dist_2,2)) / ( 2 * dist_1 * dist_2 );

  return( acos(angle) );
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
// Description:
//   Compute the vector product of (site_1-site_2) X (site_1-site_3).
////////////////////////////////////////////////////////////////////////////
void MolFramework::getUnitNormalToPlane( float *unit_normal, int site_1, 
						       int site_2, int site_3 ){

  float vector_1[3];
  float vector_2[3];

  for( int a = 0; a < 3; a++ ){
    vector_1[a] = site_info[site_2].coords[a] - site_info[site_1].coords[a];
    vector_2[a] = site_info[site_3].coords[a] - site_info[site_1].coords[a];
  }

  unit_normal[X] = ( vector_1[Y] * vector_2[Z] - 
		     vector_1[Z] * vector_2[Y] );

  unit_normal[Y] = ( vector_1[Z] * vector_2[X] - 
		     vector_1[X] * vector_2[Z] );
  
  unit_normal[Z] = ( vector_1[X] * vector_2[Y] - 
		     vector_1[Y] * vector_2[X] );

  float magnitude = sqrt( pow(unit_normal[X],2) + pow(unit_normal[Y],2) + pow(unit_normal[Z],2) );

  unit_normal[X] /= magnitude;
  unit_normal[Y] /= magnitude;
  unit_normal[Z] /= magnitude;

  return;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
float dotProduct( float *vector_1, float *vector_2 ){

  float sum = 0.0;
  for( int a = 0; a < 3; a++ )
    sum += vector_1[a] * vector_2[a];

  return( acos(sum) );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This function is used when computing the hydrogen bond energy when both 
//   the donor and acceptor are sp2 hybridized. The neighbors of sp2 hybridized 
//   atoms lie in a plane. The function calculates the normals from each sp2 
//   plane, and returns the angle between the normals.
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computeOutOfPlaneAngle( angle_triple *donor_list,
							  angle_triple *accpt_list ){

  float normal_1[3];
  float normal_2[3];
  float gamma = 0.0;

  getUnitNormalToPlane( normal_1, donor_list->atom_1, donor_list->atom_2, donor_list->atom_3 );
  getUnitNormalToPlane( normal_2, accpt_list->atom_1, accpt_list->atom_2, accpt_list->atom_3 );

  gamma = dotProduct(normal_1, normal_2);

  if( gamma < PI/2.0 ) 
    gamma = PI - gamma;

  return( gamma );
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isIsolated( SiteID atom_1 ){

  if( (site_info[atom_1].neighbor_list).size() )
    return(false);

  return(true);
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Given two atoms numbers, check that the atoms are in different residues. 
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isDifferentResidue( SiteID atom_1, SiteID atom_2 ){

  return ( site_info[atom_1].FIRST_chain_ID != site_info[atom_2].FIRST_chain_ID ||
	   site_info[atom_1].seq_number != site_info[atom_2].seq_number ||
	   site_info[atom_1].insert_code != site_info[atom_2].insert_code );

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Read in a list of hydrogen bonds. The format should be three white-space 
//   delimited columns listing the hydrogen atom number, the acceptor atom 
//   number, and the energy of the bond. The energy column may be blank, in 
//   which case the energy will be computed.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::readHydrogenBondList(){

  string filename = path_name + "hbonds.in";
  string linebuf;
  ifstream Hbond_list;
  ifstream old_file_name;

  // Try to open a file named "hbond.in". If the file can't be found, 
  // prompt the user for the name of the file listing their hydrogen bonds. 
  //////////////////////////////////////////////////////////////////////
  cout << "reading " << filename << endl;
  Hbond_list.open( filename.c_str() );
  
  // Give a warning if they are using hbonds.list
  //////////////////////////////////////////////////////////////////////
  if( !Hbond_list ){
    old_file_name.open("hbonds.list");
    if( !old_file_name ){
      cout << " The input file name for hbonds has changed. You should use \"hbonds.in\" NOT hbonds.list." << endl;
      exit(1);
    }
  }
 
  if( !Hbond_list ){

    clear_screen;
    cout << endl;
    cout << " The file [" << filename << "] was not found in this directory." << endl << endl;
    do{
      Hbond_list.clear();
      cout << " Please enter the name of the file that contains the list of hydrogen bonds." << endl;
      cout << " Filename = ";
      getline(cin, filename);
      
      Hbond_list.open( filename.c_str(), ios::in );
      if( !Hbond_list || did_press_enter(filename) ){
	clear_screen;
	cout << endl << "File not found." << endl << endl;
      }
      
    } while( !Hbond_list || did_press_enter(filename) );
  }

  // Read each line of the Hbond list file. The first number should
  // correspond to a hydrogen atom, and the second atom to the acceptor
  // atom. An optional third number can be used to force a user-defined
  // energy for the given bond. If no energy is given, it will be 
  // computed. 
  //////////////////////////////////////////////////////////////////////
  //  unsigned int donor = 0; // FIXME - warning: unused variable 'donor'
  unsigned int hydrogen = 0;
  unsigned int acceptor = 0;
  int bars = parameters.hb_bars;
  float energy = -23.4567;
  size_t field_start = 0;
  size_t field_end = 0;
  string number;

  while( !Hbond_list.eof() ){
    getline(Hbond_list, linebuf);
    field_start = field_end = 0;

    // Skip blank lines and lines that begin with '#'.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
	isComment(linebuf) )
      break;

    // Read in the hydrogen atom number
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    hydrogen = atoi( number.c_str() );
    
    // Read in the acceptor atom number
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    acceptor = atoi( number.c_str() );

    // OPTIONAL Read the energy of the given bond, if it is listed. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      energy = atof( number.c_str() );
    }
    
    // OPTIONAL Read in the number of bars to use for this hydrogen bond
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      bars = atoi( number.c_str() );

      if( bars < 0 ||
	  bars > 6 ){
	cout << " ERROR: Found a hydrogen bond with an invalid number of bars while reading" << endl;
	cout << "        the file " << filename << "." << endl;
	cout << "        INVALID LINE: " << linebuf << endl;
	exit(1);
      }
    }

    if (!parameters.use_first_numbering) {
			// Map the original numbers into the numbers used internally by FIRST.
			//////////////////////////////////////////////////////////////////////    
			hydrogen = orig_2_FIRST_atom_number[hydrogen];
			acceptor = orig_2_FIRST_atom_number[acceptor];
		} 

    if( hydrogen <= 0 || hydrogen > total_sites ||
	acceptor <= 0 || acceptor > total_sites ){
      cout << " ERROR: An atom in the file [" << filename << "] does not exist." << endl;
      cout << "        The error occurred with the following line:" << endl;
      cout << "        " << linebuf << endl << endl;
      exit(1);
    }

    hydrogen_bonds.push_back( new_bonds(hydrogen, acceptor, bars, energy) );
  }

  if( parameters.verbose )
    cout << "Read " << hydrogen_bonds.size() << " hydrogen bonds from an external file." << endl;
  Hbond_list.close();
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::readHydrophobicTetherList(){

  string filename = path_name + "hphobes.in";
  string linebuf;
  ifstream Hphobic_list;
  ifstream old_file_name;

  // If we are in interactive mode, prompt the user for the name of 
  // the file listing their hydrogen bonds. If we're in non-interactive
  // mode, look for a file named hphobes.in. If this file does not 
  // exist, exit. 
  //////////////////////////////////////////////////////////////////////
  Hphobic_list.open( filename.c_str() );

  // Give a warning if they are using hbonds.list
  //////////////////////////////////////////////////////////////////////
  if( !Hphobic_list ){
    old_file_name.open("hphobes.list");
    if( !old_file_name ){
      cout << " The input file name for hphobes has changed. You should use \"hphobes.in\" NOT hbonds.list." << endl;
      exit(1);
    }
  }

  if( !Hphobic_list ){
    clear_screen;
    cout << endl;
    cout << " The file " << filename << " was not found in this directory." << endl << endl;
    do{
      Hphobic_list.clear();
      cout << " Please enter the name of the file that contains the list of hydrophobic tethers." << endl;
      cout << " Filename = ";
      getline(cin, filename);
      
      Hphobic_list.open( filename.c_str(), ios::in );
      if( !Hphobic_list || did_press_enter(filename) ){
	clear_screen;
	cout << endl << "File not found." << endl << endl;
      }
      
    } while( !Hphobic_list || did_press_enter(filename) );
  }
  
  // Read each line of the Hphobic list file. The first number should
  // correspond to a hydrogen atom, and the second atom to the accpt
  // atom. An optional third number can be used to force a user-defined
  // energy for the given bond. If no energy is given, it will be 
  // computed. 
  //////////////////////////////////////////////////////////////////////
  unsigned int atom_1 = 0;
  unsigned int atom_2 = 0;
  int bars   = parameters.hp_bars;
  size_t field_start = 0;
  size_t field_end = 0;
  string number;

  while( !Hphobic_list.eof() ){
    getline(Hphobic_list, linebuf);

    // Skip blank lines.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
	isComment(linebuf) )
      break;

    // Read in site 1 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    atom_1 = atoi( number.c_str() );
    
    // Read in site 2 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    atom_2 = atoi( number.c_str() );
    
    // OPTIONAL Read the number of bars to use in the constraint, if it is listed. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      if( number.find("123456789") == string::npos)
	bars = atoi( number.c_str() );	
      else
	bars = parameters.hp_bars;

      if( bars < 0 ||
	  bars > 6 ){
	cout << " ERROR: Found a hydrophobic tether with an invalid number of bars while reading" << endl;
	cout << "        the file " << filename << "." << endl;
	cout << "        INVALID LINE: " << linebuf << endl;
	exit(1);
      }
    }

		if (!parameters.use_first_numbering) {
      // Map the original numbers into the numbers used internally by FIRST.
      //////////////////////////////////////////////////////////////////////    
      atom_1 = orig_2_FIRST_atom_number[atom_1];
      atom_2 = orig_2_FIRST_atom_number[atom_2];
		} 

    if( atom_1 <= 0 || atom_1 > total_sites ||
	atom_2 <= 0 || atom_2 > total_sites ){
      cout << " ERROR: An atom in the file [" << filename << "] does not exist." << endl;
      cout << "        The error occurred with the following line:" << endl;
      cout << "        " << linebuf << endl << endl;
      exit(1);
    }

    hydrophobic_tethers.push_back( new_bonds(atom_1, atom_2, bars) );
  }

  if( parameters.verbose )
    cout << "Read " << hydrophobic_tethers.size() << " hydrophobic tethers from external file." << endl;
  Hphobic_list.close();

  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::readCovalentBondList(){

  string filename = path_name + "cov.in";
  string linebuf;
  ifstream cov_list;

  cov_list.open( filename.c_str() );

  //////////////////////////////////////////////////////////////////////
  if( !cov_list ){
    if( parameters.interactive ){
      clear_screen;
      cout << endl;
      cout << " The file " << filename << " was not found in this directory." << endl << endl;
      do{
	cov_list.clear();
	cout << " Please enter the name of the file that contains the list of covalent bonds." << endl;
	cout << " Filename = ";
	getline(cin, filename);
	
	cov_list.open( filename.c_str(), ios::in );
	if( !cov_list || did_press_enter(filename) ){
	  clear_screen;
	  cout << endl << "File not found." << endl << endl;
	}
	
      } while( !cov_list || did_press_enter(filename) );
    }
    else{
      cout << " File " << filename << " not found." << endl;
      exit(1);
    }
  }
  
  // Read each line of the Hphobic list file. The first number should
  // correspond to a hydrogen atom, and the second atom to the accpt
  // atom. An optional third number can be used to force a user-defined
  // energy for the given bond. If no energy is given, it will be 
  // computed. 
  //////////////////////////////////////////////////////////////////////
  SiteID atom_1 = 0;
  SiteID atom_2 = 0;
  int bars   = parameters.cov_bars;
  size_t field_start = 0;
  size_t field_end = 0;
  int counter = 0;
  string number;

  while( !cov_list.eof() ){
    getline(cov_list, linebuf);

    // Skip blank lines.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
	isComment(linebuf) )
      break;

    // Read in site 1 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", 0);
    field_end   = linebuf.find_first_of(" \t\n", field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    atom_1 = atoi( number.c_str() );
    
    // Read in site 2 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    atom_2 = atoi( number.c_str() );
    
    // OPTIONAL Read the number of bars to use in the constraint, if it is listed. 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of(" \t\n", field_end);
    field_end   = linebuf.find_first_of(" \t\n", field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      if( number.find_first_of("123456789") == string::npos)
	bars = parameters.cov_bars;
      else
	bars = atoi( number.c_str() );

      if( bars < 0 ||
	  bars > 6 ){
	cout << " ERROR: Found a covalent bond with an invalid number of bars while reading" << endl;
	cout << "        the file " << filename << "." << endl;
	cout << "        INVALID LINE: " << linebuf << endl;
	exit(1);
      }
    }

    if (!parameters.use_first_numbering) {
			// Find the internally used FIRST atom number that maps to the
			// original atom number. 
			//////////////////////////////////////////////////////////////////////    
			atom_1 = orig_2_FIRST_atom_number[atom_1];
			atom_2 = orig_2_FIRST_atom_number[atom_2];
		} 

    if( atom_1 <= 0 || atom_1 > total_sites ||
	atom_2 <= 0 || atom_2 > total_sites ){
      cout << " ERROR: An atom in the file [" << filename << "] does not exist." << endl;
      cout << "        The error occurred with the following line:" << endl;
      cout << "        " << linebuf << endl << endl;
      exit(1);
    }

    add_to_site_info_array( atom_1, atom_2, bars );
    counter++;
  }

  if( parameters.verbose )
    cout << "Read " << counter << " covalent bonds from external file named " << filename << "." << endl;
  cov_list.close();

  return;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Assign each site that will be included in the network to a unique cell 
//   of a 3-dimensional square grid that encompasses the protein. The cell 
//   length is defined in "parameters". Assigning sites to grid cell minimizes 
//   the pairwise distance search for a given atom to those atoms found in 
//   adjacent "cells" Keep an empty shell around the protein to facilitate 
//   easier searching of edge sites.
// Parameters:
//   grid_length - Edge length of the 3D grid the atom_data coordinates are
//     stored in. This value is set in the file global_defs.h. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::assignAtomsTo3DGrid(){
  
  float grid_length = parameters.grid_length;

  float x_length = max_coords[X] - min_coords[X];
  float y_length = max_coords[Y] - min_coords[Y];
  float z_length = max_coords[Z] - min_coords[Z];

  //cout << "x length = " << x_length << endl;
  //cout << "y length = " << y_length << endl;
  //cout << "z length = " << z_length << endl;

  // 1. Determine how many cells are needed in each direction (x,y,z) to 
  //    encompass the entire structure. Then create an array of integer 
  //    lists to store the atoms that are found in each cell. 
  ////////////////////////////////////////////////////////////////////////////////
  x_cells = (int) ceil( x_length/grid_length );
  y_cells = (int) ceil( y_length/grid_length );
  z_cells = (int) ceil( z_length/grid_length );

  // If the coordinates are constant in one dimension, allocate one grid 
  // for that dimension. 
  //////////////////////////////////////////////////////////////////////
  if( x_cells == 0 ) 
    x_cells = 1;
  if( y_cells == 0 ) 
    y_cells = 1;
  if( z_cells == 0 )
    z_cells = 1;

  x_cells += 2; // Allocate an empty shell around the
  y_cells += 2; // structure. Makes it easier when
  z_cells += 2; // searching grids at the edge of the box. 

  coordinate_grid.resize( x_cells * y_cells * z_cells );
  
  // 2. Determine the 3D cell to which "current_atom" belongs, and add it to the  
  //    list of atoms ( of type vector<int> ) in that cell. 
  ////////////////////////////////////////////////////////////////////////////////
  for( unsigned int current_atom = 1; current_atom <= total_sites; current_atom++ ){

    // a. Find the current_atom's grid position.
    //////////////////////////////////////////////////////////////////////
    site_info[current_atom].grid_X = (int) ( (site_info[current_atom].coords[X] - min_coords[X])/grid_length );
    site_info[current_atom].grid_Y = (int) ( (site_info[current_atom].coords[Y] - min_coords[Y])/grid_length );
    site_info[current_atom].grid_Z = (int) ( (site_info[current_atom].coords[Z] - min_coords[Z])/grid_length );
    
    // b. Move each atom into the cell +1 in the x-, y- and z- directions, 
    //    thus making a 1-cell empty shell around the structure.
    //////////////////////////////////////////////////////////////////////
    site_info[current_atom].grid_X++;
    site_info[current_atom].grid_Y++;
    site_info[current_atom].grid_Z++;
    
    // c. Error check.
    //////////////////////////////////////////////////////////////////////
    if( site_info[current_atom].grid_X > x_cells || 
	site_info[current_atom].grid_Y > y_cells || 
	site_info[current_atom].grid_Z > z_cells ){
      cout << " ERROR: Coordinates of atom not in 3D grid." << endl;
      exit(1);
    }
    
    // d. Add the "current_atom" to the list of atoms for the corresponding 
    //    grid position.
    //////////////////////////////////////////////////////////////////////
    coordinate_grid[ gridPosition(current_atom) ].push_back( current_atom );
  }

  return;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::identifyCovalentBonds(){

  if( parameters.verbose )
    cout << "Identifying covalent bonds..." << endl;

  // Read covalent bonds from a file named "cov.in". Do NOT identify 
  // new bonds...
  //////////////////////////////////////////////////////////////////////
  if( parameters.covin ){
    readCovalentBondList();
    return;
  }

  OldCovalentBonder *oldcovalentbonder;
  CovalentBonder *covalentbonder;
  
  if( parameters.bondnew ) {
    //if( parameters.verbose ){
    //cout << endl;
    //cout << "*NOTE: Modifications were made to the default behaviour of the" << endl;
    //cout << "       identifyCovalentBonds functions. In FIRST v5.2, covalent" << endl;
    //cout << "       bonds that were not found in template files were NOT bonded" << endl;
    //cout << "       by default. Now, they ARE BONDED by default, and more sophisticated" << endl;
    //cout << "       heuristics have been added to improve bonding. In most cases" << endl;
    //cout << "       this should not affect your results, however, if you require" << endl;
    //cout << "       the previous covalent-bonding behaviour, please use the" << endl;
    //cout << "       command-line flag: \"-bondold\". This message will remain in" << endl;
    //cout << "       the source code until 12/31/2006." << endl;
    //}
    covalentbonder = new CovalentBonder(this);
    covalentbonder->autobond();
    delete covalentbonder;
  }
  else{
    oldcovalentbonder = new OldCovalentBonder(this);
    oldcovalentbonder->autobond();
    delete oldcovalentbonder;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::identifyHydrogenBonds(){

  if( parameters.skip_identify_hydrogen_bonds )
    return;

  if( parameters.verbose )
    cout << "Identifying hydrogen bonds..." << endl;

  // Read hydrogen bonds from a file. Do not look for them.
  //////////////////////////////////////////////////////////////////////
  if( parameters.hbin ){
    readHydrogenBondList();
    return;
  }

  int donor = 0;
  // int prev_donor = 0; // FIXME - warning: unused variable 'prev_donor'
  float SB_energy = 0.0;
  float HB_energy = 0.0;
  vector<int> grids_to_check;
  bool exposedA = false;
  bool exposedB = false;

  // 2. Look for potential hydrogen bonds. 
  //////////////////////////////////////////////////////////////////////
  for( SiteID current_atom = 1; current_atom <= total_sites; current_atom++ ){
    exposedA = site_info[current_atom].isExposed;
    exposedB = false;
    if( isDonorHydrogen( current_atom, &donor ) ){

      //if we're using the exposed surface routine with exposedBonds0
      //then skip any exposed atom 
      if ( parameters.exposedBonds == 0 ) {
        if ( exposedA ) continue;
      }

      // a. Create a list of local coordinate grids to check.
      //////////////////////////////////////////////////////////////////////
      grids_to_check = getGrids( current_atom );
      
      // b. Check each local coordinate grid.
      //////////////////////////////////////////////////////////////////////
      current_grid = grids_to_check.begin();
      while( current_grid != grids_to_check.end() ){
	
	// c. Get the list of the atoms in the current grid. 
	//////////////////////////////////////////////////////////////////////
	neighbor_atom = coordinate_grid[*current_grid].begin();
	while( neighbor_atom != coordinate_grid[*current_grid].end() ){
          exposedB = site_info[*neighbor_atom].isExposed;

          bool exposedGood = true;
          if ( parameters.exposedBonds == 0) {
            if ( exposedB ) {
              exposedGood = false;
            } 
          }
          else if( parameters.exposedBonds == 1 ) {
            if ( exposedB && exposedA ) {
              exposedGood = false;
            } 
          }

	  // d. Check the distance between current_atom and neighbor_atom. Only check
	  //    atoms with larger site numbers, to prevent double counting. 
	  //    Also check for hydrogen bond versus salt bridge, and compute the 
	  //    appropriate energy.
	  ////////////////////////////////////////////////////////////////////////////////
	  if( getDistance(current_atom, *neighbor_atom) < 5.0 &&
	      isHydrogenBondAcceptor( *neighbor_atom ) && 
	      exposedGood &&
	      !NthNearestNeighbor( current_atom, *neighbor_atom, 3) ){

	    if( isSaltBridge(donor, current_atom, *neighbor_atom) ){
	    	       
	      SB_energy = saltBridgeEnergy(current_atom, *neighbor_atom);
	      hydrogen_bonds.push_back( new_bonds(current_atom, *neighbor_atom, parameters.hb_bars, SB_energy, "SB" ) );
	      checkForPreviousBond( donor, current_atom, *neighbor_atom );
	    }
	    else if( isHydrogenBond(current_atom, *neighbor_atom) ){

	      HB_energy = hydrogenBondEnergy(current_atom, *neighbor_atom);
	      
	      if( HB_energy < 1000.0 ){
		hydrogen_bonds.push_back( new_bonds(current_atom, *neighbor_atom, parameters.hb_bars, HB_energy) );
		checkForPreviousBond( donor, current_atom, *neighbor_atom );
	      }

	    }

	  }
	  neighbor_atom++;
	}
	current_grid++;
      }
      
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Hydrophobic tethers can be identified in several ways. The specific 
//   function that we're going to use is stored in the parameter 
//   "hphobe_fxn".
////////////////////////////////////////////////////////////////////////////////
void MolFramework::identifyHydrophobicTethers(){

  if( parameters.skip_identify_hydrophobic_tethers )
    return;

  if( parameters.verbose )
    cout << "Identifying hydrophobic interactions by using function " << parameters.hphobe_fxn << "..." << endl;

  if( parameters.phin ){
    readHydrophobicTetherList();
    return;
  }
  
  // 1. If the user selected hydrophobic function 0, then they don't want
  //    to include tethers. Return to running the file. 
  ////////////////////////////////////////////////////////////////////////////////
  if( parameters.hphobe_fxn == 0 )
    return;

  // 2. Search for hydrophobic tethers. 
  ////////////////////////////////////////////////////////////////////////////////
  vector<int> grids_to_check;
  bool exposedA = false;
  bool exposedB = false;

  for( SiteID current_atom = 1; current_atom <= total_sites; current_atom++ ){
    exposedA = site_info[current_atom].isExposed;
    exposedB = false;
    if( isHydrophobicAtom(current_atom) ){
      if ( parameters.exposedBonds == 0 ) {
        if ( exposedA ) {
          continue;
        }
      }

      // a. Create a list of coordinate grids to check. Check each grid. 
      //////////////////////////////////////////////////////////////////////
      grids_to_check = getGrids( current_atom );
      current_grid = grids_to_check.begin();
      while( current_grid != grids_to_check.end() ){
	
	// b. Check each atom in the current grid. 
	//////////////////////////////////////////////////////////////////////
	neighbor_atom = coordinate_grid[*current_grid].begin(); // FIXME - replace with unsigned int
	while( neighbor_atom != coordinate_grid[*current_grid].end() ){
          exposedB = site_info[*neighbor_atom].isExposed;

          bool exposedGood = true;
          if ( parameters.exposedBonds == 0) {
            if ( exposedB ) {
              exposedGood = false;
            } 
          }
          else if( parameters.exposedBonds == 1 ) {
            if ( exposedB && exposedA ) {
              exposedGood = false;
            } 
          }

	  // c. Check the distance between current_atom and neighbor_atom. Only 
	  //    check atoms with larger site numbers, to prevent double counting. 
	  //////////////////////////////////////////////////////////////////////
	  if( *neighbor_atom > current_atom &&
	      isHydrophobicAtom(*neighbor_atom) &&
              exposedGood &&
	      !NthNearestNeighbor(current_atom, *neighbor_atom, 3 ) &&
	      isWithinHydrophobicCutoff(current_atom, *neighbor_atom)
	      && !in_ring_stacked_pair(current_atom, *neighbor_atom) ){

	    //For Nucleic acids. Hydrophobic Tethers may exist within the same residue.
	    if( isNucleicAcidAtom(current_atom) && isNucleicAcidAtom( *neighbor_atom ) ){
	    // Use any Hydrophobic identification function
	    //////////////////////////////////////////////////////////////////////
	      if( parameters.hphobe_fxn == 1 ||  parameters.hphobe_fxn == 2 ||  parameters.hphobe_fxn == 3){
		hydrophobic_tethers.push_back( new_bonds(current_atom, *neighbor_atom, parameters.hp_bars) );
	      }
	    // Error check if the user asked for an undefined function.
	    //////////////////////////////////////////////////////////////////////
	    else{
	      cout << " ERROR: The hydrophobic identification function you designated does not exist." << endl;
	      exit(1);
	    }
	    // For proteins. Atoms belong to different residues.
	    } else{ 
	      if(  isDifferentResidue(current_atom, *neighbor_atom) ){
		// Use  Hydrophobic identification function #1
		//////////////////////////////////////////////////////////////////////
		if( parameters.hphobe_fxn == 1 ){
		  hydrophobic_tethers.push_back( new_bonds(current_atom, *neighbor_atom, parameters.hp_bars) );
		}	    
		
		// Use Hydrophobic identification function #2
		//////////////////////////////////////////////////////////////////////
		else if( parameters.hphobe_fxn == 2 ){
		  if( isConnectedToOnlyHydrophobicAtoms(current_atom) &&
		      isConnectedToOnlyHydrophobicAtoms(*neighbor_atom) ){
		    hydrophobic_tethers.push_back( new_bonds(current_atom, *neighbor_atom, parameters.hp_bars) );
		  }
		}
		
		// Use hydrophobic identification function #3
		//////////////////////////////////////////////////////////////////////
		else if( parameters.hphobe_fxn == 3 ){
		  if( isConnectedToOnlyHydrophobicAtoms(current_atom) &&
		      isConnectedToOnlyHydrophobicAtoms(*neighbor_atom) &&
		      oneInterresidueHydrophobicTetherPerAtom(current_atom, *neighbor_atom) ){
		    hydrophobic_tethers.push_back( new_bonds(current_atom, *neighbor_atom, parameters.hp_bars) );
		  }
		}
		// Error check if the user asked for an undefined function.
		//////////////////////////////////////////////////////////////////////
		else{
		  cout << " ERROR: The hydrophobic identification function you designated does not exist." << endl;
		  exit(1);
		}
	      } //End of if(isDifferentResidue)
	    }
	    
	  }
	  neighbor_atom++;
	}
	current_grid++;
      }
    }
  }

  // For nucleic acid chain check if neighboring bases have only one tether,
  // paper of S.Fulle and H.Gohlke
  // "Analyzing the Flexibility of RNA Structures by Constraint Counting" Biophys. J., 2008, 94,4202-4219
  SiteID atom_1 = 0;
  SiteID atom_2 = 0;
  typedef vector< pair <int,int> > residueNumberMap;
  residueNumberMap residue_check_map ;
  residueNumberMap::iterator res_iter ;
  vector<new_bonds> hydrophobic_tethers_toBeErased;
  vector<new_bonds>::iterator current_tether = hydrophobic_tethers.begin();
  while( current_tether != hydrophobic_tethers.end() ){
    atom_1 = current_tether->site_1;
    atom_2 = current_tether->site_2;
    if (is_nucleobase_atom (atom_1) && is_nucleobase_atom (atom_2) ){
      Site_Info site1 = this->site_info[atom_1];
      Site_Info site2 = this->site_info[atom_2];      

      stringstream residueIDString1;
      stringstream residueIDString2;
      residueIDString1 << site1.seq_number << ";" << site1.insert_code << ";"<< site1.chain_ID;
      int res_id1 = this->unique_res_id[residueIDString1.str()];
      residueIDString2 << site2.seq_number << ";" << site2.insert_code << ";"<< site2.chain_ID;
      int res_id2 = this->unique_res_id[residueIDString2.str()];
      
      for (res_iter = residue_check_map.begin(); res_iter != residue_check_map.end(); res_iter++){
	int res_id1_saved = res_iter->first;
	int res_id2_saved = res_iter->second;
	if (res_id1_saved == res_id1 && res_id2_saved == res_id2){
	  hydrophobic_tethers_toBeErased.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
	  break;
	}
      }
	res_iter = residue_check_map.begin();
	residue_check_map.insert(res_iter, make_pair(res_id1, res_id2) );

    }
    current_tether++;
  }

  vector<new_bonds>::iterator current_tether1;

  SiteID atom_1_erase;
  SiteID atom_2_erase;
  current_tether1 = hydrophobic_tethers_toBeErased.begin();
  while( current_tether1 != hydrophobic_tethers_toBeErased.end() ){
    atom_1_erase = current_tether1->site_1;
    atom_2_erase = current_tether1->site_2;
    current_tether = hydrophobic_tethers.begin();
    while( current_tether != hydrophobic_tethers.end() ){
      atom_1 = current_tether->site_1;
      atom_2 = current_tether->site_2;
      if (atom_1 == atom_1_erase && atom_2 == atom_2_erase){
	hydrophobic_tethers.erase( current_tether );
	break;
      }
      current_tether ++;
    }
    current_tether1 ++;
  }

  //For nucleic acid chain. Now pass hydrophobic tethers that are from nuclebase pairs to stacked_rings. 
  //They are effectivelly the same that are identified in MolFramework::is_stacked_ring_pair.
  //This separation on hydrophobic_tethers and stacked_rings is done only for nucleic acid pairs.
  // Nucleic acid-protein or protein-protein pairs are identified normally, in MolFramework::is_stacked_ring_pair.
  hydrophobic_tethers_toBeErased.clear();
  current_tether = hydrophobic_tethers.begin();
  while( current_tether != hydrophobic_tethers.end() ){
    atom_1 = current_tether->site_1;
    atom_2 = current_tether->site_2;
    if ( is_nucleobase_atom(atom_1) && is_nucleobase_atom(atom_2)){
      if( !parameters.srin ) stacked_rings.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
      hydrophobic_tethers_toBeErased.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
    }
    else if ( is_nucleobase_atom(atom_1) && is_ring_atom(atom_2)){
      if( !parameters.srin ) stacked_rings.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
      hydrophobic_tethers_toBeErased.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
    }
    else if ( is_nucleobase_atom(atom_2) && is_ring_atom(atom_1)){
      if( !parameters.srin ) stacked_rings.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
      hydrophobic_tethers_toBeErased.push_back( new_bonds(atom_1, atom_2, parameters.hp_bars) );
    }
      current_tether++;
  }
  current_tether1 = hydrophobic_tethers_toBeErased.begin();
  while( current_tether1 != hydrophobic_tethers_toBeErased.end() ){
    atom_1_erase = current_tether1->site_1;
    atom_2_erase = current_tether1->site_2;
    current_tether = hydrophobic_tethers.begin();
    while( current_tether != hydrophobic_tethers.end() ){
      atom_1 = current_tether->site_1;
      atom_2 = current_tether->site_2;
      if (atom_1 == atom_1_erase && atom_2 == atom_2_erase){
	hydrophobic_tethers.erase( current_tether );
	break;
      }
      current_tether ++;
    }
    current_tether1 ++;
  }



  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Choose which atoms can partake in hydrophobic interactions. Currently, 
//   only carbon and sulfur are allowed to form hydrophobic tethers.
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isHydrophobicAtom( unsigned int this_atom ){

  if( site_info[this_atom].element_name == "C " ||
      site_info[this_atom].element_name == "S " ){
    return(true);
  }
  return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Check each atom connected to "atom_1". If any of these atoms is not a 
//   carbon "C " or sulfur "S ", return 0. To add additional atoms to check 
//   against, add an additional or (||) to the if statement below.
// Parameters:
//   atom_1 - FIRST_number of a given atom.
//////////////////////////////////////////////////////////////////////////////// 
bool MolFramework::isConnectedToOnlyHydrophobicAtoms( unsigned int atom_1 ){

  int neighbor_atom = 0;
  string neighbor_element;

  for( unsigned int neighborSiteNumber = 0; neighborSiteNumber < (site_info[atom_1].neighbor_list).size(); neighborSiteNumber++ ){
    
    neighbor_atom = site_info[atom_1].neighbor_list[neighborSiteNumber];

    if( site_info[neighbor_atom].element_name != "C " && 
	site_info[neighbor_atom].element_name != "H " &&
	site_info[neighbor_atom].element_name != "S " ){
      return(false);
    }
  }

  return(true);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isWithinHydrophobicCutoff( SiteID atom_1, SiteID atom_2 ){

  float radii_1 = 0.0;
  float radii_2 = 0.0;
  float hydrophobic_cutoff  = 0.0;

  // Modification of cutoff distance to take into account parameters in the paper of S.Fulle and H.Gohlke 
  // "Analyzing the Flexibility of RNA Structures by Constraint Counting" Biophys. J., 2008, 94,4202-4219
  if ( isNucleicAcidAtom(atom_1) && isNucleicAcidAtom(atom_2 ) ){
    hydrophobic_cutoff = parameters.PH_nucleic_cutoff;
  }else{
    hydrophobic_cutoff = parameters.PH_cutoff;
  }

  if( site_info[atom_1].element_name == "C " )
    radii_1 = 1.7;
  if( site_info[atom_2].element_name == "C " )
    radii_2 = 1.7;
  if( site_info[atom_1].element_name == "S " )
    radii_1 = 1.8;
  if( site_info[atom_2].element_name == "S " )
    radii_2 = 1.8;

  if( getDistance(atom_1, atom_2) <= (radii_1 + radii_2 + hydrophobic_cutoff) )
    return(true);
  else
    return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework:: oneInterresidueHydrophobicTetherPerAtom( SiteID current_atom, SiteID neighbor ){

  SiteID atom_1 = 0;
  SiteID atom_2 = 0;

  bool current_flag =1;

  vector<new_bonds>::iterator current_tether = hydrophobic_tethers.begin();
  while( current_tether != hydrophobic_tethers.end() ){
    atom_1 = current_tether->site_1;
    atom_2 = current_tether->site_2;

    if( ( current_atom == atom_1 &&
	  !isDifferentResidue(atom_2, neighbor) ) ||
	( neighbor == atom_2 &&
	  !isDifferentResidue(atom_1, current_atom) ) ){



      if( getDistance(current_atom, neighbor) < getDistance(atom_1, atom_2) ){
	hydrophobic_tethers.erase( current_tether );
	current_flag =1;
	break;
      }      
      else
	{
	  current_flag =0;
	  break;
	}
    }
    current_tether++;
  }
  
  // To remove dependence on how the grid is positioned on a molecule
  // we have to repeat the loop
  current_tether = hydrophobic_tethers.begin(); 
  while( current_tether != hydrophobic_tethers.end() ){
    atom_1 = current_tether->site_1;
    atom_2 = current_tether->site_2;

    if( ( current_atom == atom_1 &&
	  !isDifferentResidue(atom_2, neighbor) ) ||
	( neighbor == atom_2 &&
	  !isDifferentResidue(atom_1, current_atom) ) ){

      if( getDistance(current_atom, neighbor) < getDistance(atom_1, atom_2) ){
	hydrophobic_tethers.erase( current_tether );
	current_flag =1;
	break;
      }      
      else
      {
      current_flag =0;
      break;
      }
    }
    current_tether++;
  }
  
  return current_flag;
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Determines if an atom is part of the protein backbone (EXCLUDES carbonyl 
//   oxygen). The function relies on the name of the atom, and the names of the
//   neighboring atoms. For example, if the input atom is named "CA", and it has
//   a nearest neighbor that is named "C" and a nearest neighbor that is named "N"
//   this function will return the integer 1, indicating that the input atom was
//   a main chain alpha carbon. 
// Parameters:
//   atom_1 - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
short int MolFramework::isBackbone( SiteID atom_1 ){

  if( site_info[atom_1].atom_name == "CA" &&
      NthNearestNeighbor( atom_1, "N", 1 ) &&
      NthNearestNeighbor( atom_1, "C", 1 ) )
    return(2);
  
  else if( site_info[atom_1].atom_name == "N" && 
	   NthNearestNeighbor( atom_1, "CA", 1 ) &&
	   NthNearestNeighbor( atom_1, "C", 2 ) )
    return(1);
  
  else if( site_info[atom_1].atom_name == "C" && 
	   NthNearestNeighbor( atom_1, "CA", 1 ) &&
	   NthNearestNeighbor( atom_1, "O", 1 ) )
    return(3);

  else
    return(0);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Determines if an atom is part of the main-chain (INCLUDES carbonyl
//   oxygen. The function relies on the name of the atom, and the names of the
//   neighboring atoms. For example, if the input atom is named "CA", and it has
//   a nearest neighbor that is named "C" and a nearest neighbor that is named "N"
//   this function will return the integer 1, indicating that the input atom was
//   a main chain alpha carbon. 
//  Parameters:
//   atom_1 - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
short int MolFramework::isMainchain( SiteID atom_1 ){

  short int backbone = isBackbone( atom_1 );
  if( backbone )
    return( backbone );

  else if( site_info[atom_1].atom_name == "O" && 
	   NthNearestNeighbor( atom_1, "C", 1 ) &&
	   NthNearestNeighbor( atom_1, "CA", 2 ) )
    return(4);

  else if( site_info[atom_1].atom_name == "H" &&
	   NthNearestNeighbor( atom_1, "N", 1 ) &&
	   NthNearestNeighbor( atom_1, "CA", 2 ) )
    return(5);

  // DNA, RNA structures
  else if( site_info[atom_1].atom_name == "P" &&
	   NthNearestNeighbor( atom_1, "O5'", 1) )
    return(101);
  

  return(0);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Check some standard double bond pairs to see if this is really a double bond.
//   The most important cases will be for those atoms that can hydrogen bonds, as
//   the bond order can effect the assigned hybridization. 
// Parameters:
//   atom_1 - 
//   atom_2 -
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isDoubleBond( SiteID atom_1, SiteID atom_2 ){

  if( site_info[atom_1].element_name == "N " ){
    
    if( site_info[atom_2].element_name == "C "){
      if( getDistance( atom_1, atom_2 ) <= 1.4 ) 
	return( true );
      else
	return( false );
    }
  }

  return( false );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
// Parameters:
//   atom_1 - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
int MolFramework::isChargedResidue( SiteID atom_1 ){

  if( site_info[atom_1].residue_name == "ARG" &&
      site_info[atom_1].element_name == "N " )
    return( DONOR_CHARGED_SP2 );

  else if( site_info[atom_1].residue_name == "LYS" && 
	   site_info[atom_1].element_name == "N " )
    return( DONOR_CHARGED_SP3 );
  
  else if( (site_info[atom_1].residue_name == "ASP" ||
	    site_info[atom_1].residue_name == "GLU") &&
	   site_info[atom_1].element_name == "O " )
    return( ACCEPTOR_CHARGED_SP2 );
  
  // Check if atom_1 is a terminal main-chain carboxyl oxygen
  //////////////////////////////////////////////////////////////////////
  else if( site_info[atom_1].element_name == "O " && 
	   !NthNearestNeighbor( atom_1, "CB", 1 ) &&
 	    NthNearestNeighbor( atom_1, "CA", 2 ) &&
	   !NthNearestNeighbor( atom_1, "N", 2 ) ){
    return( ACCEPTOR_CHARGED_SP2 );
  }

  // Check if atom_1 is a terminal main-chain amide. 
  //////////////////////////////////////////////////////////////////////
  else if( site_info[atom_1].element_name == "N " && 
	   NthNearestNeighbor( atom_1, "CA", 1 ) &&
	  !NthNearestNeighbor( atom_1, "C", 1 ) )
    return( DONOR_CHARGED_SP3 );

  else
    return(0);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
// Parameters:
//   atom_1 - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
int MolFramework::getHistidineProtonationState( SiteID atom_1 ){

  int protonation_state = 0;
  int CE_atom = 0;
  SiteID neighbor_atom = 0;
  int other_nitrogen = 0;
  int nitrogen_neighbor = 0;

  // 1. Check to see if there's a proton on atom_1, and find the atom number
  //    of the " CE " side-chain atom. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int neighborSiteNumber = 0; neighborSiteNumber < (site_info[atom_1].neighbor_list).size(); neighborSiteNumber++ ){

    neighbor_atom = site_info[atom_1].neighbor_list[neighborSiteNumber];
   
    if( site_info[neighbor_atom].element_name == "H " )
      protonation_state++;

    if( site_info[neighbor_atom].atom_name == "CE1" )
      CE_atom = neighbor_atom;
  } 

  // 2. Start searching from the " CE " atom in the HIS side chain for
  //    the nitrogen that wasn't found above. Once found, check its 
  //    protonation state. 
  //////////////////////////////////////////////////////////////////////
  for( SiteID ceNeighborSiteNumber = 0; ceNeighborSiteNumber < (site_info[CE_atom].neighbor_list).size(); ceNeighborSiteNumber++ ){
  
    neighbor_atom = site_info[CE_atom].neighbor_list[ceNeighborSiteNumber];
    if( site_info[neighbor_atom].element_name == "N " &&
	neighbor_atom != atom_1 ){
      
      other_nitrogen = neighbor_atom;
      for( unsigned int otherNitrogenNeighborSiteNumber = 0; otherNitrogenNeighborSiteNumber < (site_info[other_nitrogen].neighbor_list).size(); otherNitrogenNeighborSiteNumber++ ){

	nitrogen_neighbor = site_info[other_nitrogen].neighbor_list[otherNitrogenNeighborSiteNumber];
	if( site_info[nitrogen_neighbor].element_name == "H " )
	  protonation_state += 2;
      }
    } 

  } 

  // 3. Return hbond status label.
  //////////////////////////////////////////////////////////////////////
  return( protonation_state );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
// Parameters:
//   atom_1 - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
int MolFramework::determineNitrogenHydrogenBondStatus( SiteID atom_1 ){

  // Handle Asparagine and Glutamine side chains.
  //////////////////////////////////////////////////////////////////////
  if( site_info[atom_1].residue_name == "ASN" ||
      site_info[atom_1].residue_name == "GLN" )
    return( DONOR_SP2 );

  // Handle Histidine side chains.
  //////////////////////////////////////////////////////////////////////
  else if( site_info[atom_1].residue_name == "HIS" ){

    if( getHistidineProtonationState(atom_1) == 1 )
      return( DONOR_CHARGED_SP2 );

    else if( getHistidineProtonationState(atom_1) == 2 )
      return( ACCEPTOR_CHARGED_SP2 );

    else if( getHistidineProtonationState(atom_1) == 3 )
      return( DONOR_CHARGED_SP2 );

    else
      return( ACCEPTOR_CHARGED_SP2 );
  }

  // Handle Tryptophan side chains.
  //////////////////////////////////////////////////////////////////////
  else if( site_info[atom_1].residue_name == "TRP" )
    return( DONOR_SP2 );

  // If the nitrogen does not belong to a standard amino acid, set the 
  // hydrogen bond status based on the total number of bonded neighbors 
  // and the total number of neighbors that are nitrogen. 
  // BMH 10.26.04 These default rules check only the valency to assign
  //   hybridization. Really should check the bond order of each bond
  //   to accurately determine the hybridization. 
  //////////////////////////////////////////////////////////////////////
  else{
    if( (site_info[atom_1].neighbor_list).size() == 2 ){
      if( isDoubleBond(atom_1, site_info[atom_1].neighbor_list[0] ) )
	return( ACCEPTOR_SP2 );
      else if( isDoubleBond(atom_1, site_info[atom_1].neighbor_list[1] ) )
	return( ACCEPTOR_SP2 );
      else
	return( ACCEPTOR_CHARGED_SP2 );
    }
    else if( (site_info[atom_1].neighbor_list).size() == 3 )
      return( DONOR_SP2 );
    else if( (site_info[atom_1].neighbor_list).size() == 4 ){
      if( parameters.energyFxnA )
	return( DONOR_CHARGED_SP2 );
      else
	return( DONOR_CHARGED_SP3 );
    }
    else if( (site_info[atom_1].neighbor_list).size() >= 5 ){
      // char write_atom_number[6];
      //hy36encode(5,site_info[atom_1].orig_atom_number, write_atom_number);
      cout << "Nitrogen " << site_info[atom_1].orig_atom_number << " has 5 or more covalently bonded neighbors." << endl;
      exit(2);
    }
  }

  // Default. Return Donor SP2.
  //////////////////////////////////////////////////////////////////////
  return( DONOR_SP2 );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isTerminalNitrogen( SiteID nitrogen ){
  
  if( site_info[nitrogen].neighbor_list.size() == 4 )
    return true;

  return false;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
// Parameters:
//   oxygen - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
int MolFramework::determineOxygenHydrogenBondStatus( SiteID oxygen ){

  //cout << "    Determining hydrogen-bond status for oxygen: " << oxygen << endl;
  // Handle Serine and Threonine side chains.
  //////////////////////////////////////////////////////////////////////
  if( site_info[oxygen].residue_name == "SER" ||
      site_info[oxygen].residue_name == "THR" )
    return( DONOR_AND_ACCEPTOR_SP3 );
  
  // Handle Tyrosine side chains.
  //////////////////////////////////////////////////////////////////////
  else if( site_info[oxygen].residue_name == "TYR" )
    return( DONOR_AND_ACCEPTOR_SP2 );

  // Handle Asparagine and Glutamine side chains.
  //////////////////////////////////////////////////////////////////////
  else if( site_info[oxygen].residue_name == "ASN" ||
	   site_info[oxygen].residue_name == "GLN" )
    return( ACCEPTOR_SP2 );

  // Handle an oxygen from a non-standard amino acid.
  //////////////////////////////////////////////////////////////////////  
  else{

    // if the oxygen has only one covalent bond, see if it's part of a 
    // carboxylate, or not. 
    //////////////////////////////////////////////////////////////////////
    if( (site_info[oxygen].neighbor_list).size() == 1 ){              // If the current oxygen has one neighbor, check

      if( NthNearestNeighbor( oxygen, "O ", 2 ) ){                  // if there is another oxygen within two bonds
	int nearest_oxygen =  NthNearestNeighbor( oxygen, "O ", 2 );

	if( (site_info[nearest_oxygen].neighbor_list).size() == 1 )   // if there is, then this
	  return( ACCEPTOR_CHARGED_SP2 );                             // is a charged carboxylate moetiy.
	else
	  return( ACCEPTOR_SP2 );                                     // Otherwise it's an sp2 hybridized acceptor. 
      }
      else
	return( ACCEPTOR_SP2 );
    }
    // if the oxygen has two neighbors, one of which is a hydrogen, set it 
    // to donor and acceptor.
    //////////////////////////////////////////////////////////////////////
    else if( (site_info[oxygen].neighbor_list).size() == 2 ){

      if( getElementName(site_info[oxygen].neighbor_list[0]) == "H " ||
	  getElementName(site_info[oxygen].neighbor_list[1]) == "H " )
	return( DONOR_AND_ACCEPTOR_SP3 ); 
      else // If the oxygen has two bonds, neither of which is hydrogen, it's an SP3 acceptor.
	return( ACCEPTOR_SP3 );
    }
    // If the oxygen is isolated.
    //////////////////////////////////////////////////////////////////////
    else if( (site_info[oxygen].neighbor_list).size() == 0 ){
      return(0);
    }
    // If the oxygen has 3 covalent bonds, something went wrong. 
    //////////////////////////////////////////////////////////////////////
    else{ 
      stringstream warning;
      warning << "WARNING: Ternary Oxygen found at site " << oxygen << " (original number = " 
	      <<  site_info[oxygen].orig_atom_number<< ")" << endl;
      warning << " " << site_info[site_info[oxygen].neighbor_list[0]].FIRST_number <<" "
	      << site_info[site_info[oxygen].neighbor_list[0]].atom_name << " " 
	      << getDistance( oxygen,site_info[oxygen].neighbor_list[0] ) <<  endl << " " 
	      << site_info[site_info[oxygen].neighbor_list[1]].FIRST_number <<" "
	      << site_info[site_info[oxygen].neighbor_list[1]].atom_name << " " 
	      << getDistance( oxygen,site_info[oxygen].neighbor_list[1] ) << endl << " " 
	      << site_info[site_info[oxygen].neighbor_list[2]].FIRST_number << " "
	      << site_info[site_info[oxygen].neighbor_list[2]].atom_name << " " 
	      << getDistance( oxygen,site_info[oxygen].neighbor_list[2] ) << endl
	      << endl;
      add_warning( 3, warning.str() );
      return(0);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
// Parameters:
//   sulfur - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
int MolFramework::determineSulfurHydrogenBondStatus( SiteID sulfur ){

  if( site_info[sulfur].residue_name == "MET" )
    return( ACCEPTOR_SP3 );
  
  // If it's a Cysteine residue, mark as an sp3 donor-acceptor. If it's
  // a cystine involved in a disulfide, mark it as an sp3 acceptor only. 
  //////////////////////////////////////////////////////////////////////
  else if( site_info[sulfur].residue_name == "CYS" ||
	   site_info[sulfur].residue_name == "CYX" ){
    if( NthNearestNeighbor( sulfur, "H ", 1 ) ){
      return( DONOR_AND_ACCEPTOR_SP3 );
    }
    else{
      return( DONOR_SP3 );
      //return( ACCEPTOR_SP3 ); 
    }
  }

  // Process sulfur atom in an unknown group. 
  //////////////////////////////////////////////////////////////////////
  else if( (site_info[sulfur].neighbor_list).size() == 1 )
    return( ACCEPTOR_CHARGED_SP2 );
  else if( (site_info[sulfur].neighbor_list).size() == 2 ){
    if( NthNearestNeighbor( sulfur, "H ", 1 ) )
      return( DONOR_AND_ACCEPTOR_SP3 );
    else
      return( ACCEPTOR_SP3 );
  }
  else if( (site_info[sulfur].neighbor_list).size() == 3 )
    return( DONOR_CHARGED_SP3 );
  else{
    stringstream warning;
    warning << "Warning: Set hyrdogen-bond status for sulfur atom " << site_info[sulfur].orig_atom_number << " to 0 (no hydrogen bonding)." << endl;
    add_warning( 1, warning.str() );
    return(0);
  }
    
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   count the covalent neighbors of a site
//   and how many of them are hydrogens
//   these affect parameters for surface areas
////////////////////////////////////////////////////////////////////////////////
void MolFramework::countNeighbors(){
  string element;
  //element = site_info[current_atom].element_name;
  SiteID other_atom;
  int nHydrogen;
  int nNeighbors;
  for( SiteID current_atom = 1; current_atom <= total_sites; current_atom++ ){
    nHydrogen = 0;
    nNeighbors = site_info[current_atom].neighbor_list.size() ;
    for ( unsigned int whichN = 0; whichN < site_info[current_atom].neighbor_list.size();
          whichN++ ) {
      other_atom = site_info[current_atom].neighbor_list.at(whichN);
      element = site_info[other_atom].element_name;
      if( element == "H " ){
        nHydrogen++;
      }
    }
    site_info[current_atom].nCovalentNeighbors = nNeighbors;
    site_info[current_atom].nHydrogens = nHydrogen;
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Check the hydrogen bonding status of each atom that is being included 
//   in the analysis. Main-chain and side-chain donors and acceptors are 
//   assigned explicitly. Unknown groups are assigned based on topology.
// Parameters:
//   atom_1 - FIRST_number of an atom in the molecule.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::assignHydrogenBondStatus(){

  string element;

  // 1. Assign a hydrogen bond status (donor - acceptor - both; charged - not charged), 
  //    and a hybridization (sp2 - sp3) to each nitrogen, oxygen, and sulfur.
  //////////////////////////////////////////////////////////////////////
  for( SiteID current_atom = 1; current_atom <= total_sites; current_atom++ ){
    if( site_info[current_atom].neighbor_list.size() ){
      
      element = site_info[current_atom].element_name;
      
      // Process nitrogens
      //////////////////////////////////////////////////////////////////////
      if( element == "N " ){
	//cout << "checking N " << current_atom;
	if( isMainchain( current_atom) ) {
	  if( parameters.energyFxnA ){
	    site_info[current_atom].hbond_status = DONOR_SP2; 
	  }
	  else{
	    if( isTerminalNitrogen(current_atom) )
	      site_info[current_atom].hbond_status = DONOR_CHARGED_SP3; 
	    else
	      site_info[current_atom].hbond_status = DONOR_SP2; 
	  }
	}
	else if( isChargedResidue( current_atom) )
	  site_info[current_atom].hbond_status = isChargedResidue( current_atom);
	else
	  site_info[current_atom].hbond_status = determineNitrogenHydrogenBondStatus( current_atom );
      }
      
      // Process oxygens
      //////////////////////////////////////////////////////////////////////
      else if( element == "O " ){
	//cout << "checking O " << current_atom << " " << site_info[current_atom].atom_name;
	if( isMainchain(current_atom) ) 
	  site_info[current_atom].hbond_status = ACCEPTOR_SP2; 
	else if( isChargedResidue(current_atom) )
	  site_info[current_atom].hbond_status = isChargedResidue(current_atom);
	else
	  site_info[current_atom].hbond_status = determineOxygenHydrogenBondStatus( current_atom );
      }
      
      // Process sulfurs
      //////////////////////////////////////////////////////////////////////
      else if( element == "S " ){
	//cout << "checking S " << current_atom;
	site_info[current_atom].hbond_status = determineSulfurHydrogenBondStatus( current_atom );
      }
      
      //cout << "status " << current_atom << " " << site_info[current_atom].hbond_status <<  endl;
      //cout << setw(5) << site_info[current_atom].orig_atom_number << setw(8) << site_info[current_atom].hbond_status << endl;
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   For "atom_1", first check to see if it is a hydrogen, then check to 
//   see if the atom it is connected to has been assigned a non-zero
//    hbond_status. If the hbond_status is zero, it can not form a hydrogen 
//    bond.
// Parameters:
//   hydrogen - FIRST_number of a hydrogen atom in the molecule.
//   donor - FIRST_number of the donor atom to which hydrogen is connected.
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isDonorHydrogen( int hydrogen, int *donor ){

  if( site_info[hydrogen].element_name != "H " ||
      isIsolated(hydrogen) )
    return(false);
  
  *donor = site_info[hydrogen].neighbor_list[0];
  
  if( site_info[*donor].hbond_status == DONOR_SP2 ||
      site_info[*donor].hbond_status == DONOR_SP3 ||
      site_info[*donor].hbond_status == DONOR_AND_ACCEPTOR_SP2 ||
      site_info[*donor].hbond_status == DONOR_AND_ACCEPTOR_SP3 ||
      site_info[*donor].hbond_status == DONOR_CHARGED_SP2 ||
      site_info[*donor].hbond_status == DONOR_CHARGED_SP3 )
    return(true);
  else
    return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isHydrogenBondAcceptor( unsigned int accpt ){

  if( site_info[accpt].hbond_status == ACCEPTOR_SP2 ||
      site_info[accpt].hbond_status == ACCEPTOR_SP3 ||
      site_info[accpt].hbond_status == DONOR_AND_ACCEPTOR_SP2 ||
      site_info[accpt].hbond_status == DONOR_AND_ACCEPTOR_SP3 ||
      site_info[accpt].hbond_status == ACCEPTOR_CHARGED_SP2 ||
      site_info[accpt].hbond_status == ACCEPTOR_CHARGED_SP3 )
    return(true);
  else
    return(false);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isSaltBridge( SiteID donor, SiteID hydrogen, SiteID acceptor ){

  // 1. See if both the donor and acceptor are charged. The number 16 refers to the smallest sum
  //    of the #define'd constants for CHARGED donors and acceptors in **global_defs.h**
  //////////////////////////////////////////////////////////////////////
  if( (site_info[donor].hbond_status + site_info[acceptor].hbond_status) >= 24 ){
    
    // 2. Make sure the pair satisfy the geometric cutoff.
    //////////////////////////////////////////////////////////////////////
    if( getDistance(hydrogen, acceptor) <= parameters.cutoff_SB_hyd_accpt_dist  &&
	getDistance(donor, acceptor)    <= parameters.cutoff_SB_donor_accpt_dist  &&
	(rad2deg( computeAngle(donor, hydrogen, acceptor) ) >= parameters.cutoff_SB_donor_hyd_accpt_angle) ){
     
      for( SiteID a = 0; a < site_info[acceptor].neighbor_list.size(); a++ ){

	SiteID baseAtom = site_info[acceptor].neighbor_list[a];
	if( rad2deg( computeAngle(hydrogen, acceptor, baseAtom) ) <= parameters.cutoff_SB_hyd_accpt_base_angle ){
	  return(false);
	}
      }
      
      return(true);
    }
    
  }  

  return(false);
}

////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isHydrogenBond( SiteID hydrogen, SiteID acceptor ){

  // 1. Determine the donor atom number (assumes hydrogen only have one 
  //    covalent bond)
  //////////////////////////////////////////////////////////////////////
  int donor = site_info[hydrogen].neighbor_list[0];


  // 2. Make sure the pair satisfy the geometric cutoff.
  //////////////////////////////////////////////////////////////////////
  if( getDistance(hydrogen, acceptor) <= parameters.cutoff_HB_hyd_accpt_dist  &&
      getDistance(donor, acceptor)    <= parameters.cutoff_HB_donor_accpt_dist  &&
      rad2deg(computeAngle(donor, hydrogen, acceptor)) >= parameters.cutoff_HB_donor_hyd_accpt_angle )
    return(true);
  else
    return(false); 
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Compute the energy for salt bridges.
//
//   The salt bridge function has the following form: (L-J 10-12 type)
//                                                                    
//                 {   ( R_s )12      ( R_s )10 }                       
//     E_sb = 10.0*{ 5*(-----)   -  6*(-----)   }                       
//                 (   (R + a)        (R + a)   }                       
//                                                                    
//   where: R = donor-acceptor distance.                                
//          The prefactor 10.0 is the well depth.                       
//          R_s is the donor-acceptor distance at the minimum.          
////////////////////////////////////////////////////////////////////////////////
float MolFramework::saltBridgeEnergy( int hydrogen, int acceptor ){

  float R_s = 3.2;
  float a = 0.375;

  // 1. Determine the donor atom number (assumes hydrogen only have one 
  //    covalent bond)
  //////////////////////////////////////////////////////////////////////
  int donor = site_info[hydrogen].neighbor_list[0];

  hbond_data << "\t ------------------------------------------------------------" << endl;
  hbond_data << "\t COMPUTING ENERGY FOR SALT BRIDGE (D-H - A): " 
	     << site_info[donor].orig_atom_number << " " 
	     << site_info[hydrogen].orig_atom_number << " " 
	     << site_info[acceptor].orig_atom_number << endl;

  // 2. Measure the energy.
  //////////////////////////////////////////////////////////////////////
  hbond_data << "\t ANGLE TERMS" << endl << endl;  
  hbond_data << "\t\t COMPUTING THETA ANGLE: (donor-hydrogen-acceptor angle)" << endl;
  hbond_data << "\t\t ------------------------------" << endl;
  float theta = computeAngle( donor, hydrogen, acceptor );
  hbond_data << "\t\t THETA = " << rad2deg(theta) << endl << endl;
  

  float DA_distance = getDistance(donor, acceptor);

  float ratio = (R_s / ( DA_distance + a) );
  float ratio_squared = pow(ratio,2);

  float energy = 10.0 * ( 5*pow(ratio_squared,6) - 6*pow(ratio_squared,5) );

  hbond_data << "\t DISTANCE TERMS " << endl << endl;
  hbond_data << "\t\t DONOR-ACCEPTOR DISTANCE = " << DA_distance << endl << endl;

  hbond_data << "\t TOTAL ENERGY = " <<  energy << endl;
  hbond_data << "\t ------------------------------------------------------------" << endl << endl;

  return( energy );
}

////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isOnlyBondedToHydrogen( SiteID atom_1 ){

  if( site_info[atom_1].residue_name == "HOH" )
    return( true );
  else
    return( false );

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Find a set of three consecutively bonded atoms centered on the input atom.
//   These three atoms will form an angle, which is used in the energy functions. 
//   If the input atom is singly coordinated, search for all triples in which the
//   the input atom is a terminal vertex on the two-edge chain. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::findAllAtomTriples( int site1, vector<angle_triple> *atom_list ){

  int site2 = 0;
  int site3 = 0;

  if( (site_info[site1].neighbor_list).size() == 1 ){

    site2 = site_info[site1].neighbor_list[0];
    for( unsigned int neighborSiteNumber = 0; neighborSiteNumber < (site_info[site2].neighbor_list).size(); neighborSiteNumber++ ){
      site3 = site_info[site2].neighbor_list[neighborSiteNumber];
      
      if( site3 != site1 ){
	if( isOnlyBondedToHydrogen( site1 ) )
	  atom_list->push_back( angle_triple(site1, site2, site3) );
	else if( site_info[site2].element_name != "H " )
	  atom_list->push_back( angle_triple(site1, site2, site3) );
      }

    }
  } 
  else{
    for( unsigned int secondNeighborSiteNumber = 0; secondNeighborSiteNumber < (site_info[site1].neighbor_list).size(); secondNeighborSiteNumber++ ){
      site2 = site_info[site1].neighbor_list[secondNeighborSiteNumber];
      for( unsigned int thirdNeighborSiteNumber = 0; thirdNeighborSiteNumber < (site_info[site1].neighbor_list).size(); thirdNeighborSiteNumber++ ){ 
	site3 = site_info[site1].neighbor_list[thirdNeighborSiteNumber];
	
	if( site2 != site3 ) {
	  if( isOnlyBondedToHydrogen( site1 ) )
	    atom_list->push_back( angle_triple(site1, site2, site3) );
	  else if( site_info[site2].element_name != "H " )
	    atom_list->push_back( angle_triple(site1, site2, site3) );
	}

      }
    }
  }

  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The hydrogen bond energy consists of a distance and an angle term. The 
//   angle term is a 10-12 Leonnard-Jones type potential, like the one used 
//   for salt bridges. However, there are four different dis. functions, 
//   one for each combination of SP2 or SP3 hybridization of donor and 
//   acceptor atoms. Please see manual for details.
//
//   For acceptors, compute the energy, as necessary, for all the incident 
//   base atoms. Return the lowest energy. This shouldn't occur too often, 
//   as most of the hydrogen bond acceptors are carbonyl oxygens, however, 
//   some cases will occur (such as the alcohol group of serine amino acids).
//
//   Hydrogen bonds are not included in site_info[].neighbor_list. 
//
//   BMH: TODO - in the case of sp3-sp2 bonds, we only need the acceptor
//        and base atom, but we cycle over all acceptor-base1-base2 triples.
//        This is fine, we always get the same answer, but may be a 
//        problem if we start to use this function many times. 
////////////////////////////////////////////////////////////////////////////////
float MolFramework::hydrogenBondEnergy( int hydrogen, int acceptor ){

  bool D_SP2 = false;
  bool D_SP3 = false;
  bool A_SP2 = false;
  bool A_SP3 = false;
  
  vector<angle_triple> donor_list;
  vector<angle_triple> accpt_list;

  // 1. Find the donor atom number (assumes hydrogen only have one 
  //    covalent bond, which should be true.)
  //////////////////////////////////////////////////////////////////////
  int donor = site_info[hydrogen].neighbor_list[0];
  hbond_data << "\t ------------------------------------------------------------" << endl;
  hbond_data << "\t COMPUTING ENERGY FOR HYDROGEN BOND (D-H - A): " 
	     << site_info[donor].orig_atom_number << " " 
	     << site_info[hydrogen].orig_atom_number << " " 
	     << site_info[acceptor].orig_atom_number << endl;

  // 2. Store all sets of three atoms that form an angle and can possibly 
  //    be used to compute the hydrogen bond energy. 
  ////////////////////////////////////////////////////////////////////////////////
  findAllAtomTriples( hydrogen, &donor_list );
  findAllAtomTriples( acceptor, &accpt_list );

  // 3. Compute the distance dependent term of the total energy.
  //////////////////////////////////////////////////////////////////////
  float R_s = 2.8;

  float DA_distance = getDistance(donor, acceptor);
  float ratio = (R_s / DA_distance );
  float ratio_squared = pow(ratio,2);
  float E_distance = 8.0 * ( 5*pow(ratio_squared,6) - 6*pow(ratio_squared,5) );

  hbond_data << "\t DISTANCE TERMS " << endl << endl;
  hbond_data << "\t\t DONOR-ACCEPTOR DISTANCE = " << DA_distance << endl;
  hbond_data << "\t\t E_distance = " << E_distance << endl << endl;

  // 4. Compute the angular prefactor [cos^2(theta) * e^(-(pi - theta)^6)].
  //    All of the hybridization-dependent angular terms contain this prefactor.
  //////////////////////////////////////////////////////////////////////
  float theta = computeAngle( donor, hydrogen, acceptor );
  float prefactor = pow(cos(theta), 2) * exp( -pow((PI -theta),6) );

  // 5. Determine the M.O. hybridization of the donor and acceptor atoms.
  ////////////////////////////////////////////////////////////////////////////////
  if( site_info[donor].hbond_status == DONOR_SP2 ||
      site_info[donor].hbond_status == DONOR_AND_ACCEPTOR_SP2 ||
      site_info[donor].hbond_status == DONOR_CHARGED_SP2 )
    D_SP2 = true;
  else
    D_SP3 = true;

  if( site_info[acceptor].hbond_status == ACCEPTOR_SP2 ||
      site_info[acceptor].hbond_status == DONOR_AND_ACCEPTOR_SP2 ||
      site_info[acceptor].hbond_status == ACCEPTOR_CHARGED_SP2 )
    A_SP2 = true;
  else
    A_SP3 = true;

  // 6. Compute the additional angular terms, as necessary, for the 4 donor-
  //    acceptor hybridization pairs. 
  ////////////////////////////////////////////////////////////////////////////////

  float E_angular = 0.0;    
  float E_ang_best = 0.0;
  float phi = 0.0;
  float gamma = 0.0;
  bool  found_unacceptable_phi = false;
  bool  found_acceptable_phi = false;

  // SP2 donor - SP2 acceptor hydrogen bond.
  //////////////////////////////////////////////////////////////////////
  if( D_SP2 && A_SP2 ){

    hbond_data << "\t ANGLE TERMS: HYBRIDIZATION: SP2-SP2" << endl << endl;  
    hbond_data << "\t\t COMPUTING THETA ANGLE: (donor-hydrogen-acceptor angle)" << endl;
    hbond_data << "\t\t ------------------------------" << endl;
    hbond_data << "\t\t THETA = " << rad2deg(theta) << endl << endl;
    hbond_data << "\t\t PREFACTOR = " << prefactor << endl << endl;

    hbond_data << "\t\t COMPUTING PHI AND GAMMA ANGLES:" << endl;
    hbond_data << "\t\t ------------------------------" << endl;

    int count = 1;
    for( unsigned int donorSiteNumber = 0; donorSiteNumber < donor_list.size(); donorSiteNumber++ ){
      for( unsigned int acceptorSiteNumber = 0; acceptorSiteNumber < accpt_list.size(); acceptorSiteNumber++ ){

	hbond_data << "\t\t Measurement " << count++ << endl;

	phi   = computeAngle( hydrogen, accpt_list[acceptorSiteNumber].atom_1, accpt_list[acceptorSiteNumber].atom_2 );
	gamma = computeOutOfPlaneAngle( &donor_list[donorSiteNumber], &accpt_list[acceptorSiteNumber] );

	if( parameters.energyFxnA ){
	  if( rad2deg(phi) > 90.0 )
	    found_acceptable_phi = true;
	}
	else{
	  if( rad2deg(phi) <= 90.0 )
	    found_unacceptable_phi = true;
	}

	if( phi > gamma ){
	  E_angular = phi;
	  hbond_data << "\t\t *PHI   = " << rad2deg(phi) << endl;
	  hbond_data << "\t\t  GAMMA = " << rad2deg(gamma) << endl << endl;
	}
	else{
	  E_angular = gamma;
	  hbond_data << "\t\t  PHI   = " << rad2deg(phi) << endl;
	  hbond_data << "\t\t *GAMMA = " << rad2deg(gamma) << endl << endl;
	}

	if( E_angular > E_ang_best )
	  E_ang_best = E_angular;
      }
    }

    if( parameters.energyFxnA ){
      if( !found_acceptable_phi ){
	hbond_data << "\t\t DISCARDING THE POTENTIAL HYDROGEN BOND." << endl;
	hbond_data << "\t\t NO PHI ANGLE LESS THAN 90 DEGREES WAS FOUND." << endl;
	return(1001.0);
      }
    }
    else{
      if( found_unacceptable_phi ){
	hbond_data << "\t\t DISCARDING THE POTENTIAL HYDROGEN BOND." << endl;
	hbond_data << "\t\t NO PHI ANGLE LESS THAN 90 DEGREES WAS FOUND." << endl;
	return(1001.0);
      }
    }

    E_angular = cos(E_ang_best)*cos(E_ang_best) * prefactor;
    hbond_data << "\t\t E_angular = " << cos(E_ang_best)*cos(E_ang_best) * prefactor << endl << endl;
  }  

  // SP3 donor - SP2 acceptor hydrogen bond.
  ////////////////////////////////////////////////////////////////////////////////
  if( D_SP3 && A_SP2 ){
    
    float best_phi = 0.0;

    hbond_data << "\t ANGLE TERMS: HYBRIDIZATION: SP3-SP2" << endl << endl;  
    hbond_data << "\t\t COMPUTING THETA ANGLE: (donor-hydrogen-acceptor angle)" << endl;
    hbond_data << "\t\t ------------------------------" << endl;
    hbond_data << "\t\t THETA = " << rad2deg(theta) << endl << endl;

    hbond_data << "\t\t COMPUTING PHI ANGLE: (hydrogen-acceptor-base atom angle)" << endl;
    hbond_data << "\t\t ------------------------------" << endl;

    int count = 1;
    hbond_data << "\t\t " << donor_list.size() << " Donor Atoms Found" << endl;
    hbond_data << "\t\t " << accpt_list.size() << " Acceptor Atoms Found" << endl << endl;

    for( unsigned int acceptorNumber = 0; acceptorNumber < accpt_list.size(); acceptorNumber++ ){
      phi   = computeAngle( hydrogen, accpt_list[acceptorNumber].atom_1, accpt_list[acceptorNumber].atom_2 );

      if( parameters.energyFxnA ){
	if( rad2deg(phi) > 90.0 )
	  found_acceptable_phi = true;
      }
      else{
	if( rad2deg(phi) <= 90.0 )
	  found_unacceptable_phi = true;
      }

      hbond_data << "\t\t Measurement " << count++ <<  " " << hydrogen << " " << accpt_list[acceptorNumber].atom_1 << " " << accpt_list[acceptorNumber].atom_2 << endl;
      hbond_data << "\t\t PHI   = " << rad2deg(phi) << endl << endl;
      hbond_data << "\t\t Energy = " << prefactor * pow( cos(phi), 2) * E_distance << endl << endl;
      if( phi > best_phi )
	best_phi = phi;
    }

    if( parameters.energyFxnA ){
      if( !found_acceptable_phi){
	hbond_data << "\t\t DISCARDING THE POTENTIAL HYDROGEN BOND." << endl;
	hbond_data << "\t\t NO PHI ANGLE LESS THAN 90 DEGREES WAS FOUND." << endl;
	return(1001.0);
      }
    }
    else{
      if( found_unacceptable_phi){
	hbond_data << "\t\t DISCARDING THE POTENTIAL HYDROGEN BOND." << endl;
	hbond_data << "\t\t PHI ANGLE LESS THAN 90 DEGREES WAS FOUND." << endl;
	return(1001.0);
      }
    }
        
    hbond_data << endl << "\t\t PHI   = " << rad2deg(best_phi) << endl;
    E_angular = prefactor * pow( cos(best_phi), 2);

    hbond_data << "\t\t E_angular = " << E_angular << endl << endl;
  }

  // SP2 donor - SP3 acceptor hydrogen bond.
  ////////////////////////////////////////////////////////////////////////////////
  if( D_SP2 && A_SP3 ){
    
    hbond_data << "\t ANGLE TERMS: HYBRIDIZATION: SP2-SP3" << endl << endl;  
    hbond_data << "\t\t COMPUTING THETA ANGLE: (donor-hydrogen-acceptor angle)" << endl;
    hbond_data << "\t\t ------------------------------" << endl;
    hbond_data << "\t\t THETA = " << rad2deg(theta) << endl << endl;
    E_angular = prefactor * prefactor;

    hbond_data << "\t\t E_angular = " << E_angular << endl << endl;
  }

  // SP3 donor - SP3 acceptor hydrogen bond.
  ////////////////////////////////////////////////////////////////////////////////
  if( D_SP3 && A_SP3 ){

    const float angleDiff = deg2rad(109.50);
    float best_phi = 100.0;

    hbond_data << "\t ANGLE TERMS: HYBRIDIZATION: SP3-SP3" << endl << endl;  
    hbond_data << "\t\t COMPUTING THETA ANGLE: (donor-hydrogen-acceptor angle)" << endl;
    hbond_data << "\t\t ------------------------------" << endl;
    hbond_data << "\t\t THETA = " << rad2deg(theta) << endl << endl;
    hbond_data << "\t\t PREFACTOR = " << prefactor << endl << endl;

    hbond_data << "\t\t COMPUTING PHI ANGLE: (hydrogen-acceptor-base atom angle)" << endl;
    hbond_data << "\t\t ------------------------------" << endl;

    int count = 1;
    for( unsigned int donorNumber = 0; donorNumber < donor_list.size(); donorNumber++ ){
      for( unsigned int acceptorNumber = 0; acceptorNumber < accpt_list.size(); acceptorNumber++ ){
	
	phi = computeAngle( hydrogen, accpt_list[acceptorNumber].atom_1, accpt_list[acceptorNumber].atom_2 );
	
	if( !parameters.energyFxnA ){
	  if( (phi -angleDiff) > deg2rad(90.0) )
	    found_unacceptable_phi = true;
	}
	
	hbond_data << "\t\t Measurement " << count++ << endl;
	hbond_data << "\t\t PHI   = " << rad2deg(phi) << endl;
	hbond_data << "\t\t Hbond Energy = " <<  prefactor * pow( cos(phi - angleDiff), 2 ) * E_distance << endl << endl;

	if( parameters.energyFxnA ){
	  if( fabs(phi-angleDiff) < fabs(best_phi-angleDiff) )
	    best_phi = phi;
	}
	else{
	  if( phi < best_phi )
	    best_phi = phi;
	}
      }  
    }

    if( !parameters.energyFxnA ){
      if( found_unacceptable_phi ){
	hbond_data << "\t\t DISCARDING THE POTENTIAL HYDROGEN BOND." << endl;
	hbond_data << "\t\t PHI ANGLE LESS THAN 90 DEGREES WAS FOUND." << endl;
	return(1001.0);
      }
    }

    E_angular = prefactor * pow( cos(best_phi - angleDiff), 2 );
    hbond_data << "\t\t E_angular = " << E_angular << endl << endl;
  }
  ////////////////////////////////////////////////////////////////////////////////

  hbond_data << "\t TOTAL ENERGY = " <<  E_distance * E_angular << endl;
  hbond_data << "\t ------------------------------------------------------------" << endl << endl;

  if( parameters.output_level < 3 )
    hbond_data.seekp(0, ios::beg);

  return( E_distance * E_angular );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Search the list of hydrogen bonds to see if there is already a bond
//   between this donor and accpetor. If so, delete the bond with the higher 
//   energy.
// Parameters:
//   this_donor - FIRST_number of a donor atom for a given hydrogen bond.
//   this_hydro - FIRST_number of a hydrogen atom for a given hydrogen bond.
//   this_accpt - FIRST_number of a acceptor atom for a given hydrogen bond.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::checkForPreviousBond( int this_donor, int this_hydro, int this_accpt ){

  int prev_donor = 0;
  int prev_hydro = 0;
  int prev_accpt = 0;

  vector<new_bonds>::iterator prev_bond = hydrogen_bonds.begin();

  while( prev_bond < hydrogen_bonds.end() ){
    prev_donor = site_info[prev_bond->site_1].neighbor_list[0];
    prev_hydro = prev_bond->site_1;
    prev_accpt = prev_bond->site_2;
    
    if( this_donor == prev_donor &&
	this_hydro != prev_hydro &&
	this_accpt == prev_accpt ){

      if( prev_bond->bond_ranking >= hydrogen_bonds.back().energy ){
	hydrogen_bonds.erase( prev_bond );
      }
      else{
	prev_bond = hydrogen_bonds.end();
	hydrogen_bonds.erase( prev_bond );
      }
    }

    prev_bond++;
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   A depth-first search from atom_1 through the covalent-bond network to find 
//   atom_2. The funtion searches recursively through the number of bonds set by
//   the variable degree. No record of the sites that have been checked is kept,
//   which can cause a performance drop when searching over loops, as it will 
//   recheck sites. However, for degree < 5 (5 is the size of the smallest common 
//   ring found in proteins), it works faster than if information about checked sites
//   is kept. 
//   The atoms that are found in the search are not necessarily those with 
//   exactly the degree specified -- any atom with a separation less than or
//   equal to the given degree will match.
//   Also note that atom_1 = atom_2 will not return a positive match for degree 
//   < 2, but will for degree >= 2, since atoms are considered to be their own
//   second-nearest neighbor.
////////////////////////////////////////////////////////////////////////////////
unsigned int MolFramework::NthNearestNeighbor( SiteID atom_1, SiteID atom_2, 
                                                  int degree, SiteID root_node ){

//  int atom = 0; // FIXME - warning: unused variable 'atom'

  if( degree == 0 ) 
    return(0);
  
  vector<unsigned int>::iterator next_atom = (site_info[atom_1].neighbor_list).begin();
  while( next_atom != (site_info[atom_1].neighbor_list).end() ){

    if( *next_atom == atom_2 )
      return(*next_atom);
    else if( *next_atom == root_node )
      next_atom++;
    else if( NthNearestNeighbor(*next_atom, atom_2, degree-1, atom_1) )
      return( atom_2 );
    else
      next_atom++;
  }
  
  return(0);
}

////////////////////////////////////////////////////////////////////////////////
unsigned int MolFramework::NthNearestNeighbor( SiteID atom_1, string atom_2,
                                                  int degree, SiteID root_node ){

  unsigned int atom = 0;

  if( degree == 0 ) 
    return(0);
  
  vector<unsigned int>::iterator next_atom = (site_info[atom_1].neighbor_list).begin();
  while( next_atom != (site_info[atom_1].neighbor_list).end() ){

    if( site_info[*next_atom].atom_name == atom_2 )
      return(*next_atom);
    else if( *next_atom == root_node )
      next_atom++;
    else{
      atom = NthNearestNeighbor( *next_atom, atom_2, degree-1, atom_1 );
      if(atom)
	return(atom);
      else
	next_atom++;
    }

  }
  
  return(0);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Special function to determine the element type of atoms in PDB formatted 
//   files. The function takes advantage of the 4-column width atom-name field
//   in which the element name occupies the first two columns. Special checks
//   are coded for hydrogen atoms, which are an exception to the element field
//   rules in PDB files. 
////////////////////////////////////////////////////////////////////////////////
string MolFramework::getElementNamePDBFormat( string &atom_name, bool ambiguous_atom_type ){

  string element;

  size_t start1 = atom_name.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");

  if( start1 == string::npos ){
    cout << " ERROR: Empty atom name string read from file." << endl;
    exit(1);
  }

  if( start1 == 0 )
    element = atom_name.substr( 0,2 );
  else if( start1 == 1 )
    element = atom_name.substr( 1,1 ) + " ";
  else{
    element = atom_name.substr( start1,1 ) + " ";
    stringstream warning;
    warning << "Atom name " << atom_name << " is not in PDB file format." << endl;
    add_warning( 2, warning.str() );
    warning.clear();
  }

  // By default, assign any atom name that has an H in the first column 
  // to be a hydrogen atom. 
  //////////////////////////////////////////////////////////////////////
  if( atom_name[0] == 'H' )
    element = "H ";
  
  // Look for potential two-character element names. Flag any ambiguous 
  // atoms that may need a topology check. 
  //////////////////////////////////////////////////////////////////////
  
  // Check FE
  //////////////////////////////////////////////////////////////////////
  if( element == "F " && 
      atom_name.find("FE") != string::npos )
    element = "FE";

  // Check MG
  //////////////////////////////////////////////////////////////////////
  else if( element == "M " ){
    if( atom_name.find("MG") != string::npos )
      element = "MG";
    if( atom_name.find("MN") != string::npos )
      element = "MN";
  }

  // Check LI
  //////////////////////////////////////////////////////////////////////
  else if( element == "L " && 
      atom_name.find("LI") != string::npos )
    element = "LI";
  
  // Check ZN
  //////////////////////////////////////////////////////////////////////
  else if( element == "Z " && 
      atom_name.find("ZN") != string::npos )
    element = "ZN";
  
  //ofstream atom_to_element_list("atom2element.txt", ios::app);
  //atom_to_element_list << setw(8) << "[" << atom_name << "]" << setw(8) << "[" << element << "]" << endl;
  //atom_to_element_list.close();
  
  return( toUpper(element) );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This function returns a 2-character string containing the name of the element 
//   corresponding to the input atom name. If the element contains only one 
//   character, the character is stored in the first position of the element string.
//   For example, if the atom name CE1 is given, element = "C ". 
// Parameters:
//   atom_name - The atom_name data element of an atom_data object.
////////////////////////////////////////////////////////////////////////////////
string MolFramework::getElementName( string &atom_name, bool ambiguous_atom_type ){
  
  string element;

  // remove leading whitespace and numbers. Read up to the next white
  // space of number (we want to read "N A" as "N"). 
  //////////////////////////////////////////////////////////////////////
  size_t start = atom_name.find_first_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
  if( start == string::npos ){
    cout << "Error: atom_name [" << atom_name << "] empty in getElementName()" << endl;
    exit(1);
  }
  size_t end = atom_name.find_first_of(" 1234567890\n", start);
  if( end == string::npos ) 
    end = atom_name.size();

  string chars_only = atom_name.substr( start, end-start );

  // For atom name longer than one character, set the element as the 
  // first character. This should cover, C, H, N, O, P, and S. Then
  // do a check for possible two letter elements.
  //////////////////////////////////////////////////////////////////////
  if( chars_only.size() == 1 ){
    element = chars_only + " ";

    //ofstream atom_to_element_list("atom2element.txt", ios::app);
    //atom_to_element_list << setw(8) << atom_name << setw(8) << element << endl;
    //atom_to_element_list.close();
    
    return( toUpper( element) );
  }

  else{
    element = chars_only.at(0);
    element += " ";
  }

  // Special check for hydrogen atoms. 
  //////////////////////////////////////////////////////////////////////
  if( chars_only.find("H") != string::npos )
    ambiguous_atom_type = true;

  // Look for potential two-character element names. Flag any ambiguous 
  // atoms that may need a topology check. 
  //////////////////////////////////////////////////////////////////////
  
  // Check FE
  //////////////////////////////////////////////////////////////////////
  if( element == "F " && 
      atom_name.find("FE") != string::npos )
    element = "FE";

  // Check MG
  //////////////////////////////////////////////////////////////////////
  else if( element == "M " ){
    if( atom_name.find("MG") != string::npos )
      element = "MG";
    if( atom_name.find("MN") != string::npos )
      element = "MN";
  }

  // Check LI
  //////////////////////////////////////////////////////////////////////
  else if( element == "L " && 
      atom_name.find("LI") != string::npos )
    element = "LI";
  
  // Check ZN
  //////////////////////////////////////////////////////////////////////
  else if( element == "Z " && 
      atom_name.find("ZN") != string::npos )
    element = "ZN";


  //ofstream atom_to_element_list("atom2element.txt", ios::app);
  //atom_to_element_list << setw(8) << atom_name << setw(8) << element << endl;
  //atom_to_element_list.close();

  return( toUpper(element) );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This function attempts to return the name of the element when given the 
//   "atom_name" of an atom_data object. The string is converted to all upper-
//   case letters before returned.
// Parameters:
//   atom_1 - FIRST_number of an atom. 
////////////////////////////////////////////////////////////////////////////////
string MolFramework::getElementName( unsigned int atom_1 ){
  
  string element;
  string atom_name = site_info[atom_1].atom_name;
  return( getElementName( atom_name ) );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
int MolFramework::getFIRSTNumber( int  orig_atom_number){
  return( orig_2_FIRST_atom_number[orig_atom_number] );
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Look for a bond between atom_1 and atom_2 in the bond_list passed to this 
//   function. If we find the bond, remove it from the list.
////////////////////////////////////////////////////////////////////////////////
int MolFramework::removeFromBondList( SiteID atom_1, SiteID atom_2, vector<new_bonds> &bond_list ){

  vector<new_bonds>::iterator current_bond = bond_list.begin();

  while( current_bond != bond_list.end() ){

    //cout << current_bond->site_1 << " " << atom_1 << endl;
    //cout << current_bond->site_2 << " " << atom_2 << endl;
    if( (current_bond->site_1 == atom_1 && current_bond->site_2 == atom_2 ) ||
	(current_bond->site_1 == atom_2 && current_bond->site_2 == atom_1 ) ){
      bond_list.erase( current_bond );
      return(0);
    }
    else
      current_bond++;
  }

  return(1);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   There is a public variable named dilution_list defined in the topology_object
//   class. This variable is of type vector<new_bonds>, same as hydrogen_bonds,
//   hydrophobic_tethers, and user_defined_bonds. Any hydrogen_bond, 
//   hydrophobic_tether or user_defined_bond that shall be diluted from the 
//   network is added to the dilution_list. By default, only hydrogen bonds that
//   mee tthe energy_cutoff criteria are added. When all bonds have been added
//   to dilution_list, it is sorted. The sort is from weakest to strongest, but
//   can be reversed by uncommenting the second sort option. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::prepareFrameworkForAnalysis(){

  // user-defined constraints are not diluted from the network. 
  ////////////////////////////////////////////////////////////////////////////////
  for( unsigned int constraintNumber = 0; constraintNumber < user_defined_constraints.size(); constraintNumber++ ){
    int site_1 = user_defined_constraints[constraintNumber].site_1;
    int site_2 = user_defined_constraints[constraintNumber].site_2;
    int bars   = user_defined_constraints[constraintNumber].bars;
    add_to_site_info_array( site_1, site_2, bars );
  } 

  // In interactive mode, prompt the user for a "bond_ranking" cutoff. By
  // default this will be an energy, unless the program has been modified.
  ////////////////////////////////////////////////////////////////////////////////
  if( !parameters.energy_set_on_command_line && 
      parameters.interactive &&
      (using_pdb || using_dataset || using_harlem) )
    userInteractionEnergyCutoff();

  for( unsigned int hydrogenBondNumber = 0; hydrogenBondNumber < hydrogen_bonds.size(); hydrogenBondNumber++ ){
    if( hydrogen_bonds[hydrogenBondNumber].energy <= parameters.energy_cutoff ){

      int site_1 = hydrogen_bonds[hydrogenBondNumber].site_1;
      int site_2 = hydrogen_bonds[hydrogenBondNumber].site_2;
      int bars   = hydrogen_bonds[hydrogenBondNumber].bars;

      add_to_site_info_array( site_1, site_2, bars );
      total_hbonds++;
      dilution_list.push_back( hydrogen_bonds[hydrogenBondNumber] );
      included_hbonds << setw(8) << site_info[site_1].orig_atom_number
		      << setw(8) << site_info[site_2].orig_atom_number
		      << setw(10) << hydrogen_bonds[hydrogenBondNumber].energy
		      << setw(5)  << bars
		      << endl;
    }
  } 

  for( unsigned int hydrophobicTetherNumber = 0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++ ){
    int site_1 = hydrophobic_tethers[hydrophobicTetherNumber].site_1;
    int site_2 = hydrophobic_tethers[hydrophobicTetherNumber].site_2;
    int bars   = hydrophobic_tethers[hydrophobicTetherNumber].bars;
    
    add_to_site_info_array( site_1, site_2, bars );
    total_hphobes++;
    included_hphobes << setw(8) << site_info[site_1].orig_atom_number
		     << setw(8) << site_info[site_2].orig_atom_number
		     << setw(5)  << bars
		     << endl;
  } 
  
  for( unsigned int stackedRingNumber = 0; stackedRingNumber < stacked_rings.size(); stackedRingNumber++ ){
    int site_1 = stacked_rings[stackedRingNumber].site_1;
    int site_2 = stacked_rings[stackedRingNumber].site_2;
    int bars   = stacked_rings[stackedRingNumber].bars;

    add_to_site_info_array( site_1, site_2, bars );
    total_stacked_rings++;
    included_aromats << setw(8) << site_info[site_1].orig_atom_number
		     << setw(8) << site_info[site_2].orig_atom_number
		     << setw(5)  << bars
		     << endl;
  }  
  // Sort from most positive to most negative. Uses bond_ranking member
  // of dilution_list to perform sort. 
  //////////////////////////////////////////////////////////////////////
  stable_sort( dilution_list.begin(), dilution_list.end(), greater<new_bonds>() );
  
  // Sort from most negative to most positive. Uses bond_ranking member
  // of dilution_list to perform sort. 
  //////////////////////////////////////////////////////////////////////
  //sort( dilution_list.begin(), dilution_list.end() )

  if( dilution_list.size() )
    dilution_list[0].total_in_network = dilution_list.size();

  // Delete some variables no longer required to free up memory for the
  // pebble game. 
  //////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////
void MolFramework::prepareFrameworkForFreeEnergyDilution(){

  // user-defined constraints are not diluted from the network. 
  ////////////////////////////////////////////////////////////////////////////////
  for( unsigned int userDefinedConstraintNumber = 0; userDefinedConstraintNumber < user_defined_constraints.size(); userDefinedConstraintNumber++ ){
    int site_1 = user_defined_constraints[userDefinedConstraintNumber].site_1;
    int site_2 = user_defined_constraints[userDefinedConstraintNumber].site_2;
    int bars   = user_defined_constraints[userDefinedConstraintNumber].bars;
    add_to_site_info_array( site_1, site_2, bars );
  } 

  // In interactive mode, prompt the user for a "bond_ranking" cutoff. By
  // default this will be an energy, unless the program has been modified.
  ////////////////////////////////////////////////////////////////////////////////
  if( parameters.interactive )
    userInteractionEnergyCutoff();

  for( unsigned int hydrogenBondNumber = 0; hydrogenBondNumber < hydrogen_bonds.size(); hydrogenBondNumber++ ){
    if( hydrogen_bonds[hydrogenBondNumber].bond_ranking <= parameters.energy_cutoff ){

      int site_1 = hydrogen_bonds[hydrogenBondNumber].site_1;
      int site_2 = hydrogen_bonds[hydrogenBondNumber].site_2;
      int bars   = hydrogen_bonds[hydrogenBondNumber].bars;
      add_to_site_info_array( site_1, site_2, bars );

      dilution_list.push_back( hydrogen_bonds[hydrogenBondNumber] );
    }
  } 

  for( unsigned int hydrophobicTetherNumber = 0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++ ){

    int site_1 = hydrophobic_tethers[hydrophobicTetherNumber].site_1;
    int site_2 = hydrophobic_tethers[hydrophobicTetherNumber].site_2;
    int bars   = hydrophobic_tethers[hydrophobicTetherNumber].bars;
    add_to_site_info_array( site_1, site_2, bars );

    dilution_list.push_back( hydrophobic_tethers[hydrophobicTetherNumber] );
  } 

  // Sort from most positive to most negative. Uses bond_ranking member
  // of dilution_list to perform sort. 
  //////////////////////////////////////////////////////////////////////
  sort( dilution_list.begin(), dilution_list.end(), greater<new_bonds>() );

  // Sort from most negative to most positive. Uses bond_ranking member
  // of dilution_list to perform sort. 
  //////////////////////////////////////////////////////////////////////
  //sort( dilution_list.begin(), dilution_list.end() )

  if( dilution_list.size() )
    dilution_list[0].total_in_network = dilution_list.size();

  // Delete some variables no longer required to free up memory for the
  // pebble game. 
  //////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   The following hueristic routines are used to check the topology of the 
//   structure. If the user believes that the input structure is OK, the validation
//   routine can be turned off with the command-line flag "-noval". Validations
//   that must occur, such as ensuring between 1 and 6 bars per bond, are placed
//   elsewhere in the code. Validate structure is run before any of the non-
//   covalent bond identification routines. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::validateStructure(){

  if( parameters.verbose ){
    cout << "Performing structure validation..." << endl;
  }

  // Look for redundant atom numbers in input file.
  //////////////////////////////////////////////////////////////////////
  string temp1;
  const char* errmsg;
  int int_atom_number;
  int slen;

  set<int> atomNumberRedundancyCheck;
  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){

    if( site_info[siteNumber].FIRST_number <= 99999  && !parameters.flexweb) { // PDB file format v2.2 upper bound..
      temp1 =  site_info[siteNumber].orig_atom_number;
      slen = strlen ( temp1.c_str());
      errmsg = hy36decode(5, temp1.c_str(), slen, &int_atom_number);
      if ( errmsg )
	{
	  if (isNumber( temp1 )){
	    int_atom_number = atoi( temp1.c_str() );
	    atomNumberRedundancyCheck.insert( int_atom_number );       
	    if( atomNumberRedundancyCheck.size() != siteNumber ){
	      cout << "ERROR: Duplicate atom number [" 
		   << site_info[siteNumber].orig_atom_number 
		   << "] found in the input file." << endl;
	      cout << "       NOTE: only atom numbers <= 99,999 are checked." << endl;
	      exit(1);
	    }
	  }
	  else{
	    cout << "ERROR: Atom [" << temp1 <<"] has non-numeric identifier" << endl;
	    exit(1);
	  }
	}
      else{
	atomNumberRedundancyCheck.insert( int_atom_number );
	if( atomNumberRedundancyCheck.size() != siteNumber ){
	  cout << "ERROR: Duplicate atom number [" 
	       << site_info[siteNumber].orig_atom_number 
	       << "] found in the input file." << endl;
	  cout << "       NOTE: only atom numbers <= 99,999 are checked." << endl;
	  exit(1);
	}
      }
    }
  }

  // Make sure hydrogen atoms have a multiplicity of one.
  //////////////////////////////////////////////////////////////////////
  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
  
    if( site_info[siteNumber].element_name == "H " &&
	site_info[siteNumber].neighbor_list.size() >= 2 ){
      cout << "WARNING: Hydrogen atom " << site_info[siteNumber].orig_atom_number 
	   << " is bonded to " << site_info[siteNumber].neighbor_list.size() << " atoms" << endl;
      cout << "         Keeping the shortest bond length neighbor, labeled [X], in the network:" << endl;
      
      float shortest = 1000.0;
      SiteID kept = 0;
      for( unsigned int neighbor = 0; neighbor < site_info[siteNumber].neighbor_list.size(); neighbor++ ){

	if( getDistance( siteNumber, neighbor ) < shortest ){
	  shortest = getDistance( siteNumber, site_info[siteNumber].neighbor_list[neighbor] );
	  kept = neighbor;
	}
      }
      for( unsigned int neighbor = 0; neighbor < site_info[siteNumber].neighbor_list.size(); neighbor++ ){

	if( neighbor == kept )
	  cout << "         [X] ";
	else
	  cout << "         [ ] ";
	cout << site_info[siteNumber].orig_atom_number << " "
	     << site_info[siteNumber].atom_name << " " 
	     << " -- " 
	     << site_info[site_info[siteNumber].neighbor_list[neighbor]].orig_atom_number << " " 
	     << site_info[site_info[siteNumber].neighbor_list[neighbor]].atom_name
	     << " (" << getDistance( siteNumber, site_info[siteNumber].neighbor_list[neighbor] ) << " Angstroms)" << endl;
      }
    }
  }

  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){

    // Store ID of isolated sites.
    //////////////////////////////////////////////////////////////////////
    if( site_info[siteNumber].neighbor_list.size() == 0 ){

      stringstream warning;
      warning << "Site " << site_info[siteNumber].orig_atom_number << " has no connections." << endl;
      add_warning( 2, warning.str() );
      warning.clear();
    }

    // Check for a dangling end with 6 bars. Even if the dangling end should
    // be 6-bars, ie. missing atoms, reset it so the count comes up correct. 
    //////////////////////////////////////////////////////////////////////
    if( site_info[siteNumber].neighbor_list.size() == 1 &&
	site_info[siteNumber].number_of_bars[ site_info[siteNumber].neighbor_list[0] ] == 6 ){

      int neighbor = site_info[siteNumber].neighbor_list[0];
      
      cout << "  Warning: Bond with multiplicity of 1 has six bars. Resetting to 5 bars. " << endl; 
      cout << "          " << site_info[siteNumber].orig_atom_number << " " << site_info[neighbor].orig_atom_number; 
      cout << "\t " << site_info[siteNumber].atom_name << " (" << site_info[siteNumber].element_name << ") " 
	   << site_info[neighbor].atom_name << " (" << site_info[neighbor].element_name << ") " << endl;

      site_info[siteNumber].number_of_bars[neighbor] = 5;
      site_info[neighbor].number_of_bars[siteNumber] = 5;
    }
    //////////////////////////////////////////////////////////////////////
   
    // Check all hydrogen atoms and atoms tagged as potentially mislabeled.
    // If they're properly labeled, but have extended bond lengths, 
    // generate a warning.
    //////////////////////////////////////////////////////////////////////
    if( site_info[siteNumber].element_name == "H " &&
	site_info[siteNumber].neighbor_list.size() != 0 &&
 	getDistance( siteNumber, site_info[siteNumber].neighbor_list[0] ) > 1.2 ){

      stringstream warning;
      warning << "Poor covalent bond length of " 
	      << getDistance( siteNumber, site_info[siteNumber].neighbor_list[0] ) << " for hydrogen atom " 
	      << site_info[siteNumber].orig_atom_number << " - " 
	      << site_info[siteNumber].neighbor_list[0] << endl;
      add_warning( 2, warning.str() );

      if( site_info[siteNumber].neighbor_list.size() > 1 ){
	if( parameters.verbose ){
	  cout << "Warning: Potentially mislabeled hydrogen element type for atom " << site_info[siteNumber].orig_atom_number << endl;
	  cout << "         The covalent bond distance between " << site_info[siteNumber].orig_atom_number 
	       << " and " << site_info[site_info[siteNumber].neighbor_list[0]].orig_atom_number << " is " 
	       << getDistance( siteNumber, site_info[siteNumber].neighbor_list[0] ) << " Angstroms." << endl;
	  cout << "          Reevaluating element type... ";
	}
	if( using_pdb )
	  site_info[siteNumber].element_name = getElementNamePDBFormat( site_info[siteNumber].atom_name, false );
	if( parameters.verbose ){
	  cout << "[H ] -> [" << site_info[siteNumber].element_name << "]" << endl;
	}
      }
    }
    //////////////////////////////////////////////////////////////////////
  }  

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Compare the flexibility results from two different instances of the network. 
//   The default is to compare only backbone atoms, although any comparison is 
//   possible. All the information from the pdb file is available through the 
//   pdb_data array.                           
////////////////////////////////////////////////////////////////////////////////
void MolFramework::compareTopologies(){
  
  // Compare the backbone rigidity information between two instances of 
  // the protein network. 
  //////////////////////////////////////////////////////////////////////
  for( unsigned int siteNumber = 0; siteNumber <= total_sites; siteNumber++ ){
    
    if( isBackbone(siteNumber) ){
      if( prev_topology[siteNumber].rigid_label != site_info[siteNumber].rigid_label ){
	cout << endl << "backbone flexibilities differ" << endl;
	return;
      }
    }
      
  } 

  cout << endl << "backbone flexibilities are identical" << endl;

}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::saveTopology( unsigned int hbond_index ){

  //static int counter = 1;
  //stringstream filename;
  //filename << "decomp_list_" << (dilution_list.size() -hbond_index -1);
  //fstream decomp_list( filename.str().c_str(), ios::app );

  string output_file = path_name + "decomp_list";
  ofstream decomp_list( output_file.c_str(), ios::app );
  
  // If this is the first time we are saving the topology, then prev_topology
  // will be NULL. Create a new topology array.
  //////////////////////////////////////////////////////////////////////
  if( prev_topology == NULL ){
    prev_topology = new Site_Info[total_sites+1];
    
    decomp_list << "HEADER" 
		<< setw(8) << total_residues 
		<< setw(8) << total_sites 
		<< "     "
		<< base_name
		<< endl;
  }

  // Copy site_info into prev_topology.
  //////////////////////////////////////////////////////////////////////
  for( unsigned int siteNumber = 1; siteNumber <= total_sites; siteNumber++ )
    prev_topology[siteNumber] = site_info[siteNumber];


  // Create "decomp_list" like file similar to old version of FIRST. 
  ////////////////////////////////////////////////////////////////////////////////
  if( hbond_index != dilution_list.size() ){
    int hbond = dilution_list.size() -hbond_index -1;
    int hydro = dilution_list[hbond_index].site_1;
    int donor = site_info[hydro].neighbor_list[0];
    int accpt = dilution_list[hbond_index].site_2;
    float energy = dilution_list[hbond_index].energy;

    float mean_c = 0.00;
    float frac_f = float(floppy_modes)/(3.0*total_sites);
    float frac_lrc = float(RC_atom_list[1].size())/float(total_sites);

    if( parameters.use_unpruned_mean_coord_in_stripy_plot )
      mean_c = mean_coordination;
    else
      mean_c = pruned_mean_coordination;

    decomp_list << "A:" 
		<< setw(8)  << hbond
		<< setw(10) << showpoint << setprecision(4) << energy 
		<< setw(8)  << total_residues 
		<< setw(8)  << site_info[donor].seq_number 
		<< setw(4)  << site_info[donor].atom_name 
		<< setw(12) << site_info[accpt].seq_number 
		<< setw(4)  << site_info[accpt].atom_name 
		<< setw(15) << showpoint << setprecision(5) << setiosflags(ios::fixed) << mean_c
		<< endl;
    decomp_list << "B:" 
		<< setw(10) << mean_c
		<< setw(13) << frac_lrc
		<< setw(13) << frac_f
		<< setw(11) << showpoint << setprecision(4) << energy
		<< endl;
  }
  else{
    decomp_list << "A:" << setw(8) << hbond_index << "  -0.00000                                                 " << endl;
    decomp_list << "B:                                                                   " << endl;
  }

  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){
    prev_topology[siteNumber] = site_info[siteNumber];

    decomp_list << setw(5) << site_info[siteNumber].rigid_label;
    if( !(siteNumber%20) )
      decomp_list << endl;
    else
      decomp_list << ":";
  } 

  for( unsigned int columnNumber = (total_sites-1)%20; columnNumber < 20; columnNumber++ ){
    decomp_list << "    0";
    if( !(columnNumber%19) )
      decomp_list << endl;
    else
      decomp_list << ":";
  }

  decomp_list << "END" << endl;
  decomp_list.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::noSideChainBonds(){

  int current_atom = 0;
  int neighbor = 0; 
  int this_residue = 0;
  int this_chain = 0;

  for( unsigned int a = 0; a < side_chain_list.size(); a++ ){
    current_atom = side_chain_list[a];
    
    this_residue = site_info[current_atom].seq_number;
    this_chain   = site_info[current_atom].FIRST_chain_ID;

    for( unsigned int b = 0; b < (site_info[current_atom].neighbor_list).size(); b++ ){
      
      neighbor = site_info[current_atom].neighbor_list[b];
      if( isDifferentResidue( current_atom, neighbor) )
	return(false);
      
    }
  }

  return(true);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
inline void MolFramework::setListAsPruned( vector<SiteID> atom_list ){

  for( unsigned int siteNumber = 0; siteNumber < atom_list.size(); siteNumber++ )
    site_info[ atom_list[siteNumber] ].pruned = true;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a depth-first recursive search through the site_info topology 
//   starting on site "current_atom". Only connections to sites with the same 
//   residue number and chain ID are followed. If the site is a side-chain atom, 
//   it is added to "side_chain_list" which is of type vector<unsigned int>. A list of the
//   sites that have been checked is kept to avoid cycles. 
// Parameters:
//   current_atom - 
////////////////////////////////////////////////////////////////////////////////////
void MolFramework::listSideChainAtoms( int current_atom, int current_residue, 
						    int current_chain_ID ){

  search_list.push_back( current_atom );
  site_info[current_atom].checked = true;

  if( !isMainchain(current_atom) ){
    side_chain_list.push_back( current_atom );
  }
  
  int neighbor = 0;

  for( unsigned int neighborNumber = 0; neighborNumber < (site_info[current_atom].neighbor_list).size(); neighborNumber++ ){
    neighbor = site_info[current_atom].neighbor_list[neighborNumber];
    
    if( site_info[neighbor].seq_number  == current_residue &&
	site_info[neighbor].FIRST_chain_ID == current_chain_ID &&
	!site_info[neighbor].checked ){
      listSideChainAtoms( neighbor, current_residue, current_chain_ID );
    }
  } 

}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
void MolFramework::resetSearchList(){
  
  for( SiteID siteNumber = 0; siteNumber < search_list.size(); siteNumber++ )
    site_info[ search_list[siteNumber] ].checked = false;
  
  search_list.clear();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::pruneSideChains(){

  /*int current_single_site = 0; // FIXME - unsed variables
  int total_pruned_neighbors = 0;
  int total_neighbors = 0;
  int next_site = 0;
  int neighbor = 0;*/
  int current_residue = -999;
  int current_chain_ID = 0;

  for( SiteID current_atom = 1; current_atom <= total_sites; current_atom++ ){
    
    if( !site_info[current_atom].pruned ){
      
      if( site_info[current_atom].residue_name == "PHE" ||
	  site_info[current_atom].residue_name == "HIS" ||
	  site_info[current_atom].residue_name == "TYR" ||
	  site_info[current_atom].residue_name == "TRP" ){
	
	current_residue  = site_info[current_atom].seq_number;
	current_chain_ID = site_info[current_atom].FIRST_chain_ID;

	listSideChainAtoms( current_atom, current_residue, current_chain_ID );
	if( noSideChainBonds() ){
	  setListAsPruned( side_chain_list );
	}
	side_chain_list.clear();
	resetSearchList();
      }

    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::pruneDanglingEnds(){

  int total_pruned_neighbors = 0;
  unsigned int total_neighbors = 0;
  int next_site = 0;
  int neighbor = 0;

  for( SiteID siteNumber = 1; siteNumber <= total_sites; siteNumber++ ){

    if( !site_info[siteNumber].pruned ){
       
      int current_site = siteNumber;
      do{
	total_pruned_neighbors = next_site = 0;
	total_neighbors = (site_info[current_site].neighbor_list).size();

	for( SiteID neighborNumber = 0; neighborNumber < total_neighbors; neighborNumber++ ){
	  neighbor = site_info[current_site].neighbor_list[neighborNumber];
	  if( site_info[neighbor].pruned )
	    total_pruned_neighbors++;
	  else
	    next_site = neighbor;
	}

	if( total_neighbors - total_pruned_neighbors <= 1 ){
	  site_info[current_site].pruned = true;
	  current_site = next_site;
	}
      } while( total_neighbors - total_pruned_neighbors <= 1 &&
	       total_neighbors > 0 );
      
    }
  }

}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   MolFramework member function used for debugging. Allows for printing of
//   private member functions. 
////////////////////////////////////////////////////////////////////////////////
void MolFramework::printInfo(){
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Remove the given bond from the site_info array.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::removeBond( int this_bond ){

  int site_1 = dilution_list[this_bond].site_1;
  int site_2 = dilution_list[this_bond].site_2;

  remove_from_site_info_array( site_1, site_2 );
  
  dilution_list[this_bond].in_network = false;
  dilution_list[this_bond].total_in_network--;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Put the given bond back into the network.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::addBond( int this_bond ){

  int site_1 = dilution_list[this_bond].site_1;
  int site_2 = dilution_list[this_bond].site_2;
  int bars   = dilution_list[this_bond].bars;

  add_to_site_info_array( site_1, site_2, bars );

  dilution_list[this_bond].in_network = true;
  dilution_list[this_bond].total_in_network++;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Compute the mean coordination for the network. Computed as the number of bonds
//   divided by the number of sites. Sites with the prune member set to "true" are not 
//   included in the calculation.
///////////////////////////////////////////////////////////////////////////////
void MolFramework::computeMeanCoordination( bool pruned ){

  int total_bonds = 0;
  unsigned int neighbor = 0;
  int locked_dihedral = 0;
  int total_included_sites = 0;
  mean_coordination = 0.0;

  //ofstream fort("fort14", ios::app );
  //ofstream pruned_sites("pruned_list.txt", ios::app );
  
  int temp = 0;
  for( unsigned int current_site = 1; current_site <= total_sites; current_site++ ){
    temp = 0;
    if( !(site_info[current_site].pruned) ){

      for( unsigned int b = 0; b < (site_info[current_site].neighbor_list).size(); b++ ){
	neighbor = site_info[current_site].neighbor_list[b];

	if( !(site_info[neighbor].pruned) )
	  temp++;

	if( neighbor > current_site &&
	    !(site_info[neighbor].pruned) ){
	  total_bonds++;

	  // Look for locked dihedral angles. The total number of these bonds is
	  // is used to correct the total mean coordination for better comparison
	  // to network glasses.
	  //////////////////////////////////////////////////////////////////////
	  if( site_info[current_site].number_of_bars[neighbor] == 6 )
	    locked_dihedral++;
	}
      }
      total_included_sites++;

    }
    //pruned_sites << " " << current_site << " " << temp << endl;
  }

  // Compute the mean coordination for the network, and factor in a correction 
  // for locked dihedral angles (bonds with 6 bars).
  //////////////////////////////////////////////////////////////////////

  float locked_dihedral_correction = (float) locked_dihedral / total_included_sites;
  // float corrected_ave_R = 0.0; // FIXME - warning: unused variable 'corrected_ave_R'

  if( parameters.correction_for_2_4 ){
    total_included_sites += 3*total_hphobes;
    locked_dihedral_correction = (float) locked_dihedral / total_included_sites;
    total_bonds += 3*total_hphobes;
  }

  if( pruned ){
    pruned_mean_coordination = ( 2.0*total_bonds ) / (float) ( total_included_sites );

    if( parameters.correction_for_2_4 )
      pruned_mean_coordination = (pruned_mean_coordination + 10.0*locked_dihedral_correction)/(1.0 +4.0*locked_dihedral_correction);
  }
  else{
    mean_coordination = ( 2.0*total_bonds ) / (float) ( total_included_sites );
    
    if( parameters.correction_for_2_4 )
      mean_coordination = (mean_coordination + 10.0*locked_dihedral_correction)/(1.0 +4.0*locked_dihedral_correction);
  }

  /*
  fort.setf(ios::fixed);
  fort << showpoint;
  fort << setw(6) << total_bonds*2 + total_hphobes*3
  << setw(6) << total_included_sites
  << setw(10) << setprecision(6) << pruned_mean_coordination 
  << setw(10) << setprecision(6) << locked_dihedral_correction 
  << setw(10) << setprecision(6) << corrected_ave_R
  << endl;
  fort.close();
  */

  //pruned_sites.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::initialize( int sites ){
  
  total_sites = sites;
  site_info = new Site_Info[total_sites +1];

  reset_site_info_labels();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void MolFramework::resetCheckedLabel() {

  for( unsigned int a = 0; a <= total_sites; a++ )
    site_info[a].checked = false;
}

////////////////////////////////////////////////////////////////////////////////
// Description: compute and return the Mean Squared Deviation between 
//              a pair of MolFrameworks 
//              for now, assumes that both molFrameworks have identical structure
// Returns: Mean Squared Deviation or a negative number if an error occurs
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computePairMSD(MolFramework & firstMolFramework,
                                    MolFramework & secondMolFramework,
                                    SiteSelector  & siteSelector) {
  size_t totalSites = firstMolFramework.total_sites; 
  if (totalSites != firstMolFramework.total_sites) {
    // firstMolFramework and firstMolFramework have a different number of sites
    
    return -1.0; 
  }
  
  if (totalSites < 1) {
    return -1.0;
  }
  
  float sumOfSquaredDeviations = 0.0;
  
  size_t totalMatchingSites = 0;
  
  for (SiteID siteID = 1; siteID <= totalSites; siteID++) {
    if (!siteSelector.includeSite(siteID)) {
      continue;
    }
        
    Site_Info firstSiteInfo = firstMolFramework.site_info[siteID];
    float * firstSiteCoordinates = firstSiteInfo.coords;
    
    Site_Info secondSiteInfo = secondMolFramework.site_info[siteID];
    float * secondSiteCoordinates = secondSiteInfo.coords;
    
    totalMatchingSites ++;
    
    sumOfSquaredDeviations += 
      VectorAlgebra::distanceSquared(firstSiteCoordinates,
                                     secondSiteCoordinates);
    
  }
  
  float meanSquaredDeviation = 
    sumOfSquaredDeviations / (float)totalMatchingSites;
  
  return meanSquaredDeviation;
}

////////////////////////////////////////////////////////////////////////////////
// Description: compute and return the Mean Squared Deviation between 
//              a pair of MolFrameworks 
//              for now, assumes that both molFrameworks have identical structure
// Returns: Mean Squared Deviation or a negative number if an error occurs
//
// FIXME: DRY - refactor to eliminate the need for polymorphism   
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computePairMSD(MolFramework & molFramework,
                                    std::vector<Vector> & vectorOfVectors,
                                    SiteSelector  & siteSelector) {
  size_t totalSites = molFramework.total_sites; 
  if (totalSites != molFramework.total_sites) {
    // molFramework and molFramework have a different number of sites

    return -1.0; 
  }
  
  if (totalSites < 1) {
    return -1.0;
  }
  
  float sumOfSquaredDeviations = 0.0;
  float secondSiteCoordinates[3];
  size_t totalMatchingSites = 0;
  
  for (SiteID siteID = 1; siteID <= totalSites; siteID++) {
    if (!siteSelector.includeSite(siteID)) {
      continue;
    }
        
    Site_Info firstSiteInfo = molFramework.site_info[siteID];
    float * firstSiteCoordinates = firstSiteInfo.coords;
    
    Vector secondSiteCoordinatesVector = vectorOfVectors.at(siteID);
    
    if (secondSiteCoordinatesVector == NULL_VEC) {
      std::cout << " NULL_VEC " << std::endl;
      continue;
    }
    
    //float *secondSiteCoordinates = new float[3];
    // FIXME - refactor to merge the two Vector classes
    secondSiteCoordinates[0] = secondSiteCoordinatesVector.x;
    secondSiteCoordinates[1] = secondSiteCoordinatesVector.y;
    secondSiteCoordinates[2] = secondSiteCoordinatesVector.z;    
     
    sumOfSquaredDeviations += 
      VectorAlgebra::distanceSquared(firstSiteCoordinates,
                                     secondSiteCoordinates);
    totalMatchingSites ++;
  }
  
  float meanSquaredDeviation = 
    sumOfSquaredDeviations / (float)totalMatchingSites;
  
  return meanSquaredDeviation;
}

////////////////////////////////////////////////////////////////////////////////
// Description: compute and return the Mean Squared Deviation between 
//              a pair of MolFrameworks 
//              for now, assumes that both molFrameworks have identical structure
// Returns: Mean Squared Deviation or a negative number if an error occurs
//
// FIXME: DRY - refactor to eliminate the need for polymorphism   
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computePairMSD(std::vector<Vector> & firstVectorOfVectors,
                                    std::vector<Vector> & secondVectorOfVectors,
                                    SiteSelector  & siteSelector) {
  size_t totalSites = firstVectorOfVectors.size()-1; // FIXME - use sane indices 
  if (totalSites != secondVectorOfVectors.size()-1) {
    // molFramework and molFramework have a different number of sites
    
    return -1.0; 
  }
  
  float sumOfSquaredDeviations = 0.0;
  
  size_t totalMatchingSites = 0;
  
  for (SiteID siteID = 1; siteID <= totalSites; siteID++) {
    if (!siteSelector.includeSite(siteID)) {
      continue;
    }
        
    Vector firstSiteCoordinatesVector = firstVectorOfVectors.at(siteID);
    
    if (firstSiteCoordinatesVector == NULL_VEC) {
      std::cout << " NULL_VEC " << std::endl;
      continue;
    }
    
    Vector secondSiteCoordinatesVector = secondVectorOfVectors.at(siteID);
    
    if (secondSiteCoordinatesVector == NULL_VEC) {
      std::cout << " NULL_VEC " << std::endl;
      continue;
    }
    
    // FIXME - refactor to merge the two Vector classes
    float *firstSiteCoordinates = new float[3];
    float *secondSiteCoordinates = new float[3];
    
    // FIXME - make an operator(Vector) that returns a float* (or a better 
    //         encapsulation)
    firstSiteCoordinates[0] = firstSiteCoordinatesVector.x;
    firstSiteCoordinates[1] = firstSiteCoordinatesVector.y;
    firstSiteCoordinates[2] = firstSiteCoordinatesVector.z;        
    
    secondSiteCoordinates[0] = secondSiteCoordinatesVector.x;
    secondSiteCoordinates[1] = secondSiteCoordinatesVector.y;
    secondSiteCoordinates[2] = secondSiteCoordinatesVector.z;    
    
    totalMatchingSites ++;
 
    sumOfSquaredDeviations += 
      VectorAlgebra::distanceSquared(firstSiteCoordinates,
                                     secondSiteCoordinates);
    delete [] firstSiteCoordinates;
    delete [] secondSiteCoordinates;
  }
    
  float meanSquaredDeviation = 
    sumOfSquaredDeviations / (float)totalMatchingSites;
  
  return meanSquaredDeviation;
}


////////////////////////////////////////////////////////////////////////////////
// Description: compute and return the Root Mean Squared Deviation between 
//              a pair of MolFrameworks 
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computePairRMSD(MolFramework & firstMolFramework,
                                     MolFramework & secondMolFramework,
                                     SiteSelector  & siteSelector) {
  float msd = computePairMSD(firstMolFramework,
                             secondMolFramework,
                             siteSelector);
  
  return sqrt(msd);
}

////////////////////////////////////////////////////////////////////////////////
// Description: compute and return the Root Mean Squared Deviation between 
//              a MolFramework and a std::vector<Vector>
// FIXME: DRY - refactor to eliminate the need for polymorphism 
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computePairRMSD(MolFramework & molFramework,
                                     std::vector<Vector> & vectorOfVectors,
                                     SiteSelector  & siteSelector) {
  float msd = computePairMSD(molFramework,
                             vectorOfVectors,
                             siteSelector);
  
  return sqrt(msd);
}

////////////////////////////////////////////////////////////////////////////////
// Description: compute and return the Root Mean Squared Deviation between 
//              two std::vector<Vector>s
// FIXME: DRY - refactor to eliminate the need for polymorphism 
////////////////////////////////////////////////////////////////////////////////
float MolFramework::computePairRMSD(std::vector<Vector> & firstVectorOfVectors,
				    std::vector<Vector> & secondVectorOfVectors,
				    SiteSelector  & siteSelector) {
  float msd = computePairMSD(firstVectorOfVectors,
                             secondVectorOfVectors,
                             siteSelector);
  
  return sqrt(msd);
}

////////////////////////////////////////////////////////////////////////////////
// Description: populate the mapFromChainIDtoChain 
// mapFromChainIDtoChain[chainID] -> Chain (~ vector of Residues in the Chain)
////////////////////////////////////////////////////////////////////////////////
void MolFramework::populateStructureHierarchy() {

  for (SiteID siteID = 1; siteID <= total_sites; siteID ++) { // TODO - verify that this should range from 1 to total_sites and not from 0
    Site_Info *siteInfo = &site_info[siteID];
    siteInfo->site_number = siteID;
    
    unsigned int chainID = siteInfo->FIRST_chain_ID;
    
    Chain *chain = NULL;
    if (mapFromChainIDtoChain.find(chainID) != mapFromChainIDtoChain.end()){
      chain = mapFromChainIDtoChain[chainID];
    } else {
      chain = new Chain();
      mapFromChainIDtoChain[chainID] = chain;
    }
    
    mapFromChainToChainID[chain] = chainID;
    
    unsigned int residueID = siteInfo->seq_number;
    Residue * residue = NULL;
    
    if (chain->find(residueID) != chain->end()) {
      residue = (*chain)[residueID];
    } else {
      residue = new Residue();
      (*chain)[residueID] = residue;
    }
    
    mapFromResidueToResidueID[residue] = residueID;
    mapFromResidueToChain[residue] = chain;
    
    residue->addSite(siteInfo);
    
    mapFromSiteIDToResidue[siteID] = residue;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
bool MolFramework::isBonded( SiteID atom1, SiteID atom2 ){
  
  vector<SiteID>::iterator result;
  result = find( site_info[atom1].neighbor_list.begin(),
		 site_info[atom1].neighbor_list.end(),
		 atom2 );
  if( result == site_info[atom1].neighbor_list.end() )
    return false;

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Removes the hydrophobic tethers and stacked rings
//   from the neighbor table.  Hydrogen bonds still 
//   remain in the neighbor table, along with the covalent
//   bonds.
////////////////////////////////////////////////////////////////////////////////
void MolFramework::stripNonCovalentNeighbors(){
 
  //first eliminate the phobic tethers from the system.

  //cerr << "(Not counting phobic tethers as bonds)" << endl;
  if( parameters.verbose >= 2 ){
    cout << "  Stripping hydrophobic tethers from covalent bonds..." << endl;
  }
  for (unsigned int hydrophobicTetherNumber =0; hydrophobicTetherNumber < hydrophobic_tethers.size(); hydrophobicTetherNumber++){
    int site_1 = hydrophobic_tethers[hydrophobicTetherNumber].site_1;
    int site_2 = hydrophobic_tethers[hydrophobicTetherNumber].site_2;

    remove_from_site_info_array( site_1, site_2 );
  }

  //if aromatic ring identification is on, then remove these as well
  if ( !parameters.skip_identify_aromatics ){
    if( parameters.verbose >= 2 ){
      cout << "  Stripping stacked aromatics from covalent bonds..." << endl;
    }
    for (unsigned int stackedRingNumber =0; stackedRingNumber < stacked_rings.size(); stackedRingNumber++){
      int site_1 = stacked_rings[stackedRingNumber].site_1;
      int site_2 = stacked_rings[stackedRingNumber].site_2;
  
      remove_from_site_info_array( site_1, site_2 );
    }
  }
}

bool MolFramework::isNucleicAcidAtom( SiteID atom_1 ){
    
    bool nb_atom = is_nucleobase_atom( atom_1 );
    bool ribose_atom = is_ribose_atom(  atom_1 );
    bool nucleic_backbone_atom = is_nucleic_backbone(  atom_1 );
    bool mnb_atom = is_modifiedNucleobase_atom( atom_1 );
    if (nb_atom) return true;
    else if (ribose_atom) return true;
    else if (nucleic_backbone_atom) return true;
    else if (mnb_atom) return true;

    return false;
}


