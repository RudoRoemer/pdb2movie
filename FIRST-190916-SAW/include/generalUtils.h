#ifndef _GENERAL_UTILS_
#define _GENERAL_UTILS_

#include "global_defs.h"
#include "Parameters.h"
#include "MolFramework.h"

// the Dot structure is used when finding exposed atoms
//////////////////////////////////////////////////////////////////////
struct Dot{

public:

  Vector position;
  unsigned int nearestAtom;
  double distance2;
  bool going;
  bool gone;
  Dot(){
    position = Vector(0,0,0);
    nearestAtom = 0;
    distance2 = 0.0;
    going = false;
    gone = false;
  }
  ~Dot(){};
};

void analyzeFlexibilityAndOutput( MolFramework &structure );
void buildFramework( MolFramework &structure );
void excludeAltSideChainLocation( MolFramework &structure );
void excludeSites( MolFramework &structure, bool shrink_and_remap_only = false );
void hydrogenBondDilutionEnergyOnly( MolFramework &structure );
void hydrogenBondDilutionAssemblyCodeEnergyOnly( MolFramework &structure );
void printLogo();
void printUsage();
void readCommandLineOptions( int, char **, Parameters & );
void readData( MolFramework &structure, string infile_name );
void readParametersFile( int argc, char **argv, Parameters & );
void remapOrigToFIRSTNumbers( MolFramework &structure );
void setVdwRadii( MolFramework &structure );

string toUpper( string & );
string toLower( string & );
bool sortOnSize( const vector<unsigned int> &a, const vector<unsigned int> &b );
bool sortSetsOnSize( const set<SiteID> &a, const set<SiteID> &b ); // FIXME - replace with a generic sortOnSize method

vector<string> tokenize_string( string &input, string delimiters = " ");
void removeWhiteSpacePadding( string &str );
string hr();
string next_rcd_color_as_name( bool reset = false );
string next_jmol_color_as_name( bool reset = false );
string next_svg_color_as_name( bool reset = false );
bool isComment( string &linebuf );

void  renumber_structures_from_morph_server( MolFramework &structure, MolFramework &target );

void findExposedAtoms( MolFramework &structure );

void checkDistanceBMH( MolFramework &structure );
bool isNumber(string &str);
#endif
