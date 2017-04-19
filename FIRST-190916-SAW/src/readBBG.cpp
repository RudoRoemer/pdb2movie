#include "global_defs.h"
#include "generalUtils.h"
#include "interact.h"
#include "readBBG.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description: 
//   Analyze the flexibility of a generic body-bar graph.
////////////////////////////////////////////////////////////////////////////////
void readBBG_File( MolFramework &structure ){

  // Read though the file once to detemine the number of site, and the 
  // number of bonds. 
  ////////////////////////////////////////////////////////////////////////////
  int total_sites = scanBBG_File( structure );

  // Initialize site_info array and the pebble game.
  ////////////////////////////////////////////////////////////////////////////////
  structure.initialize( total_sites );

  // Read in the body-bar graph file. The format of this file is three white-space
  // separated columns. The first two columns indicate the site numbers of an 
  // edge in the network. The third column is the number of bars that should be
  // placed between them. If the third column is blank, a default value of five
  // bars is used. 
  ////////////////////////////////////////////////////////////////////////////////
  readBBG_Data( structure );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description: Read the body-bar graph file to determine the total number of 
//   sites. Skip blank lines and lines that begin with the hash character (#).
////////////////////////////////////////////////////////////////////////////////
int scanBBG_File( MolFramework &structure ){

  int site_1 = 0;
  int site_2 = 0;
  int total_sites = 0;
  int total_bonds = 0;

  vector<int> vertices;

  string linebuf;

  // 1. Open the file. 
  //////////////////////////////////////////////////////////////////////
  fstream infile( structure.infile_name.c_str() );
  if( !infile ){
    cout << "Error: Could not open file " << structure.infile_name << "." << endl;
    exit(1);
  }

  // 2. Count the number of lines and the maximum site number. Skip blank 
  //    lines.
  //////////////////////////////////////////////////////////////////////
  while( !infile.eof() ){
    getline(infile, linebuf);

    if( linebuf.find_first_of( WHITESPC_TAB_NEWLINE, 0 ) != string::npos &&
	linebuf.at(0) != '#' ){

      vector<string> fields = tokenize_string( linebuf );
      site_1 = site_2 = 0;

      // Error check for too many numbers on a line.
      //////////////////////////////////////////////////////////////////////
      if( fields.size() > 3 ){
	cout << "Error: Too many fields in input file." << endl;
	cout << "       This line caused the error ->[" << linebuf << "]" << endl;
	exit(1);
      }

      // Check to see if this is an isolated site
      //////////////////////////////////////////////////////////////////////
      site_1 = atoi( fields.at(0).c_str() );
      if( fields.size() > 1 ){
	site_2 = atoi( fields.at(1).c_str() );
	total_bonds++;
      }
	
      if( parameters.keep_missing_bbg_vertices ){
	if( site_1 > total_sites ) total_sites = site_1;
	if( site_2 > total_sites ) total_sites = site_2;
      }
      else{

	if( site_1 > total_sites ) total_sites = site_1;
	if( site_2 > total_sites ) total_sites = site_2;
	vertices.push_back( site_1 );
	if( site_2 ) 
	  vertices.push_back( site_2 );
      }
      
    } 
  }
  infile.close();
  cout << "1. total sites " << total_sites << endl;
  // If we're not including missing site numbers in the total, get the
  // number of unique site numbers from the input file.
  //////////////////////////////////////////////////////////////////////
  if( vertices.size() ){
    stable_sort( vertices.begin(), vertices.end() );
    vector<int>::iterator pos = unique( vertices.begin(), vertices.end() );
    vertices.erase( pos, vertices.end() );
    total_sites = vertices.size();
  }
  cout << "2. total sites " << total_sites << endl;
  structure.total_bonds = total_bonds;
 
  return( total_sites );
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void readBBG_Data( MolFramework &structure ){

  size_t field_start = 0;
  size_t field_end   = 0;

  int site_1 = 0;
  int site_2 = 0;
  int bars   = 0;

  string linebuf;
  string number;

  //////////////////////////////////////////////////////////////////////
  fstream infile( structure.infile_name.c_str() );
  if( !infile ){
    cout << "Error: Could not open file " << structure.infile_name << ". " << endl;
    exit(1);
  }

  //////////////////////////////////////////////////////////////////////
  while( !infile.eof() ){
    getline(infile, linebuf);

    // Skip blank lines and lines that begin with '#'.
    //////////////////////////////////////////////////////////////////////
    if( linebuf.find_first_of(" \n1234567890", 0 ) == string::npos ||
	isComment(linebuf) )
      break;
    
    // Read in site 1 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of( WHITESPC_TAB_NEWLINE, 0);
    field_end   = linebuf.find_first_of( WHITESPC_TAB_NEWLINE, field_start); 
    number = linebuf.substr( field_start, field_end-field_start );
    site_1 = atoi( number.c_str() );
    
    // Read in site 2 of the current edge
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of( WHITESPC_TAB_NEWLINE, field_end);
    field_end   = linebuf.find_first_of( WHITESPC_TAB_NEWLINE, field_start);
    number = linebuf.substr( field_start, field_end-field_start );
    site_2 = atoi( number.c_str() );
    
    // Read the number of bars to place between site_1 and site_2. If it's
    // blank, add the default number of bars (should be 5). 
    //////////////////////////////////////////////////////////////////////
    field_start = linebuf.find_first_not_of( WHITESPC_TAB_NEWLINE, field_end);
    field_end   = linebuf.find_first_of( WHITESPC_TAB_NEWLINE, field_start);
    if( field_start != string::npos ){
      number = linebuf.substr( field_start, field_end-field_start );
      bars = atoi( number.c_str() );
      if( bars <= 0 || bars > 6 ){
	cout << " * Warning: Can not use " << bars << " between sites " << site_1 << " and " << site_2 << ".";
	cout << linebuf << endl;
	cout << " This will be set to the default number of bars (defined in global_defs.h)" << endl;
	bars = DEFAULT_NUMBER_OF_BARS;
      }
    }
    else
      bars = DEFAULT_NUMBER_OF_BARS;
    
    (structure.site_info[site_1].neighbor_list).push_back(site_2);
    (structure.site_info[site_1].number_of_bars).insert( make_pair(site_2,bars) );
    (structure.site_info[site_2].neighbor_list).push_back(site_1);
    (structure.site_info[site_2].number_of_bars).insert( make_pair(site_1,bars) );
    
    (structure.site_info[site_1].cm_label).insert( pair<int,int> (site_2, 0) );
    (structure.site_info[site_2].cm_label).insert( pair<int,int> (site_1, 0) );
    
    if( parameters.bond_dilution )
      structure.dilution_list.push_back( new_bonds(site_1, site_2) );
    
  }

  infile.close();
  
}
////////////////////////////////////////////////////////////////////////////////

