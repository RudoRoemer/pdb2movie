#include "PDB.h"
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

void PDBAtomLine::parseAtomLine( string line ) {
  stringstream ss;
  string tempString;
  int tempInt;
  double tempDouble;

  atomOrHetatm = ( line[0]=='H' ); // ATOM:0, HETATM:1

  ss << line.substr( 6, 5 );
  ss >> tempString;
  serial = tempString;
  ss.clear();
  ss.str("");
  
  name = line.substr( 12, 4 );

  altLoc = line[16];

  resName = line.substr( 17, 3 );

  chainID = line[21];

  ss << line.substr( 22, 4 );
  ss >> tempInt;
  resSeq = tempInt;
  ss.clear();
  ss.str("");

  iCode = line[26];

  double x;
  double y;
  double z;
  ss << line.substr( 30, 8 );
  ss >> x;
  ss.clear();
  ss.str("");  
  ss << line.substr( 38, 8 );
  ss >> y;
  ss.clear();
  ss.str("");  
  ss << line.substr( 46, 8 );
  ss >> z;
  ss.clear();
  ss.str("");
  position = Vec3(x,y,z);

  if ( line.size() >= 60 ) {  
    ss << line.substr( 54, 6 );
    ss >> tempDouble;
    occupancy = tempDouble;
    ss.clear();
    ss.str("");
  }
  else occupancy = 0.0;
  
  if ( line.size() >= 66 ) {  
    ss << line.substr( 60, 6 );
    ss >> tempDouble;
    tempFactor = tempDouble;
    ss.clear();
    ss.str("");
  }
  else tempFactor = 0.0;

  if ( line.size() >= 76 ) { 
    segID = line.substr( 72, 4 );
    if ( segID.find_first_not_of(' ') == segID.npos ) segID = "";
  }
  else segID = "";

  if ( line.size() >= 78 ) { 
    ss << line.substr( 76, 2 );
    ss >> element;
    ss.clear();
    ss.str("");
    for ( size_t i = 0; i < element.size(); i++ ) {
      if ( !isalpha(element[i]) ) {
        cout << "PDB Error: Element field \"" << element << "\" is non-alphabetic" << endl;
        exit(0);
      }
    }
  }
  else element = "";
  
  if ( line.size() >= 80 ) {   
    charge = line.substr( 78, 2 );
    if ( charge.find_first_not_of(' ') == charge.npos ) charge = "";
  }
  else charge = "";
  
      
}

void ConectRecord::parseLine( string line ) {
  stringstream ss;
  if ( line.size() >= 16 ) {  
    ss << line.substr( 6, 5 );
    ss >> baseAtom;
    ss.clear();
    ss.str("");
  }
  size_t lineSize = line.size();
  size_t startIndex = 11;
  size_t fieldWidth = 5;
  string tempAtomID;
  int count = 0;
  while ( startIndex + fieldWidth <= lineSize && count < 4 ) {
    ss.clear();
    ss.str("");
    ss << line.substr( startIndex, fieldWidth );
    // if this atom id is pure spaces, then don't count it.  The line is done.
    if ( !(ss >> tempAtomID) ) break;
    neighbors.push_back(tempAtomID);
    startIndex += fieldWidth;
    count++;
  }
  //idiom to trim the vector's memory allocation
  //to the minimum necessary
  vector<string>(neighbors).swap(neighbors);

}

ostream& operator<< (ostream& os, const PDBAtomLine& atomLine ) {
  os << showpoint << setiosflags(ios::fixed);
  if (atomLine.atomOrHetatm) os << left << setw(6) << "HETATM"; // 0-5
  else os << left << setw(6) << "ATOM";
  os << right << setw(5) << atomLine.serial << left; // 6-10
  os << " "; // 11
  os << atomLine.name; // 12-15
  os << atomLine.altLoc; // 16
  os << atomLine.resName; // 17-19
  os << " "; // 20
  os << atomLine.chainID; // 21
  os << right << setw(4) << atomLine.resSeq << left; // 22-25
  os << atomLine.iCode; // 26
  os << setw(3) << ""; // 27-29
  os << setw(8) << setprecision(3)<< atomLine.position.x; // 30-37
  os << setw(8) << setprecision(3)<< atomLine.position.y; // 38-45
  os << setw(8) << setprecision(3)<< atomLine.position.z; // 46-53
  os << setw(6) << setprecision(2)<< atomLine.occupancy; // 54-59
  os << setw(6) << setprecision(2)<< atomLine.tempFactor; // 60-65
  os << setw(6) << ""; // 66-71
  os << setw(4) << atomLine.segID; // 72-75
  os << right << setw(2) << atomLine.element << left; // 76-77
  os << setw(2) << atomLine.charge; // 78-79
  return os;
}

ostream& operator<< (ostream& os, const ConectRecord& conectRecord ) {
  os << "CONECT" << right << setw(5) << conectRecord.baseAtom;
  for ( size_t i = 0; i < conectRecord.neighbors.size(); i++ ) {
    os << right << setw(5) << conectRecord.neighbors[i]; 
  }
  return os;
}

PDB::PDB()
{
}

PDB::~PDB()
{
}

void PDB::read( string filename ) {
  //first clear all vectors
  ifstream pdbfile( filename.c_str(), ios::in );
  if (!pdbfile) {
    cout << "Could not open PDB file: " << filename << endl;
    exit(1);
  }
  int nAtoms = 0;
  string currentline;
  while (!pdbfile.eof()) {
    getline( pdbfile, currentline );
    LineInfo lineInfo;
    PDBAtomLine atomLine;
    ConectRecord conectRecord;
    lineInfo.clear();
    if ( (currentline.substr(0,4) == "ATOM" ||
          currentline.substr(0,6) == "HETATM" ) &&
          currentline.size() >= 54 ) {
      lineInfo.isAtomLine = true;
      lineInfo.atomIndex = nAtoms++;
      lookupLineInfo.push_back( lineInfo );
      
      atomLine.parseAtomLine( currentline );
      atomLines.push_back( atomLine );
      
    }
    else if ( currentline.substr(0,3) == "TER" ) {
      lineInfo.isTerLine = true;
      lookupLineInfo.push_back( lineInfo );
    }
    else if ( currentline.substr(0,6) == "CONECT" ) {
      conectRecord.parseLine( currentline );
      conectRecords.push_back( conectRecord );
    }
    else continue;
  }
  pdbfile.close(); 
}

void PDB::write( string filename ) {
  ofstream outfile( filename.c_str(), ios::out );
  for ( size_t lineIndex=0; lineIndex<lookupLineInfo.size(); lineIndex++ ) {
    int atomIndex;
    if (lookupLineInfo[lineIndex].isAtomLine) {
      atomIndex = lookupLineInfo[lineIndex].atomIndex;
      outfile << atomLines[atomIndex] << '\n';
    }
    else if (lookupLineInfo[lineIndex].isTerLine) {
      outfile << "TER\n";
    } 
  }
  for ( size_t i = 0; i < conectRecords.size(); i++ ) {
    outfile << conectRecords[i] << '\n';
  }
  outfile << "END\n";
  outfile.flush();
  outfile.close();
}

