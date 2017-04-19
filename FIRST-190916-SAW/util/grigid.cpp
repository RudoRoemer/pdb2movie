// SAW Jan 4 2007
//using tokenize_string function of BMH
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <bitset>
#include <typeinfo>
#include <algorithm>
#include <new>
#include <sys/types.h>
#include <cctype>

using namespace std;

vector<string> tokenize_string( string &input, string delimiters = " " );
int readGhosts( string filename );
int writeGhosts( string filename );

bool isInGhost( unsigned int atomID, unsigned int ghostID );

struct Ghost{
  public:

  vector< unsigned int > atom;

  Ghost() {
    atom.clear();
  }
};

vector<Ghost> ghostIn, ghostOut;
string inFilename, outFilename;
unsigned int ghostA, ghostB;

int main(int argc, char **argv ) {
  if ( argc < 5 ) {
    cerr << " Usage: grigid inputFilename outputFilename ghostA ghostB" << endl;
    cerr << " Function: Reads a ghost file from inputFilename." << endl;
    cerr << " Rigidifies ghostA and ghostB into a single ghost." << endl;
    cerr << " Writes a new ghost file to outputFilename." << endl;
    exit(1);
  }

  inFilename = argv[1];
  outFilename = argv[2];

  readGhosts( inFilename );

  ghostA = atoi( argv[3] );
  ghostB = atoi( argv[4] );

  if ( ghostB == ghostA ) {
    cerr << "ERROR: ghostA = ghostB! " << ghostA << ", " << ghostB << endl;
    exit(1);
  }

  if ( ghostB < ghostA ) {
    unsigned int temp = ghostA;
    ghostA = ghostB;
    ghostB = temp; //swapped
  }

  if ( ghostB >= ghostIn.size() ) {
    cerr << "ERROR: requested ghost " << ghostB << " greater than max ghost ";
    cerr << ghostIn.size() -1 << endl;
    exit(1);
  }

  if ( ghostA < 1 ) {
    cerr << "ERROR: requested ghost " << ghostA << " less than 1." << endl;
    exit(1);
  }

  cerr << "Rigidifying ghost " << ghostA << " and " << ghostB << endl;

  vector< unsigned int > commonAtom;
  unsigned int tryingAtom;
  for ( unsigned int whichAtom = 0; whichAtom < ghostIn[ghostA].atom.size(); whichAtom++ ) {
    tryingAtom = ghostIn[ghostA].atom.at(whichAtom );
    if ( isInGhost(tryingAtom, ghostB ) ) {
      commonAtom.push_back( tryingAtom );
    }
  }

  cerr << "Found " << ghostIn[ghostA].atom.size() << " atoms in ghost " << ghostA << endl;
  cerr << "Found " << ghostIn[ghostB].atom.size() << " atoms in ghost " << ghostB << endl;
  cerr << "Found " << commonAtom.size() << " shared atoms." << endl;

  if ( commonAtom.size() > 0 ) {
    cerr << "Sharing: ";
    for ( unsigned int whichCommon = 0; whichCommon < commonAtom.size(); whichCommon++ ) {
      cerr << commonAtom.at(whichCommon) << " "; 
    }
    cerr << endl;
  }   

  if ( commonAtom.size() == 0 ) {
    cerr << "WARNING: no common atoms. I hope you mean it." << endl;
  }
  else if ( commonAtom.size() > 2 ) {
    cerr << "WARNING: more than two shared atoms; dubious ghost file." << endl;
  } 

  //now do the recombination 
  ghostOut.clear();
  //all ghosts up to ghostA are pushed over as is
  for ( unsigned int whichG = 0; whichG < ghostA; whichG++ ) {
    ghostOut.push_back( ghostIn.at(whichG) );
  }
  //now, ghost A gets its own atoms, plus ghostB's atoms, not duplicating shared atoms
  ghostOut.push_back( ghostIn.at(ghostA) );

  for ( unsigned int whichAtom = 0; whichAtom < ghostIn.at(ghostB).atom.size(); whichAtom ++ ) {
    tryingAtom = ghostIn.at(ghostB).atom.at(whichAtom);

    bool isShared = false;
    for ( unsigned int whichCommon = 0; whichCommon < commonAtom.size(); whichCommon++ ) {
      if ( commonAtom.at(whichCommon) == tryingAtom ) {
        cerr << "Not duplicating atom " << tryingAtom << endl;
        isShared = true;
        break;
      }
    }
    if ( isShared ) continue; //skip to next atom
    ghostOut.at(ghostA).atom.push_back( tryingAtom );
  }

  //the ghosts between A and B are passed as is
  for ( unsigned int whichG = ghostA +1; whichG < ghostB; whichG++ ) {
    ghostOut.push_back( ghostIn.at(whichG) );
  }
 
  //ghost B is skipped and remaining ghosts are passed
  for ( unsigned int whichG = ghostB +1; whichG < ghostIn.size(); whichG++ ) {
    ghostOut.push_back( ghostIn.at(whichG) );
  }
   
  writeGhosts( outFilename );

}

////////////////////////////////////////////////////////////////////////////////
// Description: Generic string tokenizer. The function takes as arguments the
//   string to be parsed and an optional string containing the delimiters. By 
//   default, the function will assume white-space characters as delimiters. The
//   function returns a vector of strings containing the tokens parsed from the
//   input string. NOTE: If your input string contains newline characters, these
//   will be treated as a delimiter. The string will continue past intervening 
//   newlines until it reaches the end of the string, as determined by string::npos.
////////////////////////////////////////////////////////////////////////////////
vector<string> tokenize_string( string &input, string delimiters ){

  // Add the newline character to the list of tokens, just in case it 
  // wasn't included.
  //////////////////////////////////////////////////////////////////////
  if( delimiters.find("\n") == string::npos )
    delimiters += "\n";

  vector<string> tokens;

  size_t start = input.find_first_not_of(delimiters,0);
  size_t end   = input.find_first_of(delimiters,start);

  if( start == string::npos )
    return( tokens );

  while( end != string::npos ){
    tokens.push_back( input.substr(start, end-start) );
    start = input.find_first_not_of( delimiters, end );
    end   = input.find_first_of( delimiters, start );
  }

  if( start != string::npos )
    tokens.push_back( input.substr(start, end-start) );

  return( tokens );
}
////////////////////////////////////////////////////////////////////////////////

int readGhosts( string filename ){

  ifstream gfile;
  gfile.open( filename.c_str() );

  if (!gfile ) {
    cerr << "Found no input file " << filename << endl;
    exit(1);
  }

  string linein;
  vector< string > splitLine;  

  unsigned int maxGhostID = 0;

  //first pass; get highest ghost id

  while ( !gfile.eof() ) {
    getline ( gfile, linein);
    splitLine = tokenize_string( linein );
    if (splitLine.size() < 2 ) continue;
    if ( splitLine[0] != "GHOST" ) continue; //ignore nonghost lines
   
    unsigned int thisGhostID = atoi( splitLine[1].c_str() ); 
    //cerr << "Found GID " << thisGhostID << endl;
    if ( thisGhostID > maxGhostID ) {
      maxGhostID = thisGhostID;
    } 
  }

  ghostIn.resize(maxGhostID + 1 );

  //rewind the file:
  gfile.clear();
  gfile.seekg( 0, ios::beg );

  //now actually get the data
  while ( !gfile.eof() ) {
    getline ( gfile, linein);
    splitLine = tokenize_string( linein );
    if (splitLine.size() < 4 ) continue; //ignore blank entries
    if ( splitLine[0] != "GHOST" ) continue; //ignore nonghost lines
   
    unsigned int thisGhostID = atoi( splitLine[1].c_str() ); 
    //cerr << "Reading Ghost ID " << thisGhostID << " ";
    unsigned int thisAtomID; 

    for ( unsigned int whichEntry = 3; whichEntry < splitLine.size(); whichEntry++ ) {
      thisAtomID = atoi( splitLine[whichEntry].c_str() );
      //cerr << "Reading atom " << thisAtomID << " ";

      ghostIn[thisGhostID].atom.push_back( thisAtomID );

    }
    //cerr << endl;   
 
  }

  gfile.close();

  for ( unsigned int whichGhost = 1; whichGhost <= maxGhostID; whichGhost++ ) {
    if ( ghostIn[whichGhost].atom.size() == 0 ) {
      cerr << "Warning: ghost " << whichGhost << " has no members!" << endl;
    }
  }

  return(0);
}



////////////////////////////////////////////////////////////////////////////////
// Description: Outputting the ghost member lists
////////////////////////////////////////////////////////////////////////////////
int writeGhosts( string filename ){
  cerr << "Outputting ghost member lists." << endl;

  unsigned int atomsLeft;
  unsigned int atomsThisLine;
  unsigned int currentAtom;
  unsigned int atomsPerLine = 10;

  ofstream gfile;
  gfile.open( filename.c_str() );

  for ( unsigned int whichGhost = 1; whichGhost < ghostOut.size() ; whichGhost++ ) {
    atomsLeft = ghostOut[whichGhost].atom.size();
    bool stillWriting = true;
    currentAtom = 0;

    if (atomsLeft == 0 ) {
      gfile << "GHOST " << whichGhost << " GEN " << endl;
      //empty ghost
      stillWriting = false;
    }

    while (stillWriting) {
      if ( atomsLeft > atomsPerLine ) {
        atomsThisLine = atomsPerLine;
        gfile << "GHOST " << whichGhost << " GEN ";
        for ( unsigned int mycount = 0; mycount < atomsThisLine; mycount++ ) {
          unsigned int thisAtom = ghostOut[whichGhost].atom.at(currentAtom);
          gfile << thisAtom << " ";
          atomsLeft--;
          currentAtom++;
        }
        gfile << endl;
      }
      else if ( atomsLeft == 0 ) {
        stillWriting = false;
      }
      else {
        atomsThisLine = atomsLeft;
        gfile << "GHOST " << whichGhost << " GEN ";
        for ( unsigned int mycount = 0; mycount < atomsThisLine; mycount++ ) {
          unsigned int thisAtom = ghostOut[whichGhost].atom.at(currentAtom);
          gfile << thisAtom << " ";
          atomsLeft--;
          currentAtom++;
        }
        gfile << endl;
        stillWriting = false;
      }
    }
  }

  gfile.close();
  return(0);
}


bool isInGhost( unsigned int atomID, unsigned int ghostID ){
  for ( unsigned int whichAtom = 0; whichAtom < ghostIn.at(ghostID).atom.size(); whichAtom ++ ) {
    if ( atomID == ghostIn.at(ghostID).atom.at(whichAtom) ) {
      return true; //atom is in ghost
    }
  }
  return false;
}

