#ifndef PDB_H_
#define PDB_H_

#include <vector>
#include <cstdio>
#include <cstdlib>
#include "Vec3.h"

class LineInfo
{
public:
  bool isAtomLine;
  bool isTerLine;
  int atomIndex;
  LineInfo() : isAtomLine(false), isTerLine(false), atomIndex(0) {}
  void clear() { isAtomLine = false; isTerLine = false; atomIndex = 0; }
};

class PDBAtomLine
{
public:
  int atomOrHetatm;
  std::string serial;
  std::string name;
  char altLoc;
  std::string resName;
  char chainID;
  int resSeq;
  char iCode;
  Vec3 position;
  double occupancy;
  double tempFactor;
  std::string segID;
  std::string element;
  std::string charge;
  friend std::ostream& operator<< (std::ostream& os, const PDBAtomLine& atomLine );
  void parseAtomLine( std::string line );
};

class ConectRecord
{
public:
  std::string baseAtom;
  std::vector<std::string> neighbors;
  friend std::ostream& operator<< (std::ostream& os, const ConectRecord& conectRecord );
  void parseLine( std::string line );
};

class PDB
{
public:
  PDB();
  PDB( std::string filename ) { read(filename); }
  ~PDB();
  void read( std::string filename );
  void write( std::string filename );
  void getPositions( std::vector<Vec3>& positions ) const {
    size_t N = atomLines.size();
    positions.resize( N );
    for ( size_t i = 0; i < N; i++ ) {
      positions[i] = atomLines[i].position;
    }
  }
  void setPositions( const std::vector<Vec3>& positions ) {
    size_t N = atomLines.size();
    if ( N != positions.size() ) {
      std::cout << "PDB Error: size of input positions does not match PDB odbject number of atoms" << std::endl;
      exit(0);
    }
    for ( size_t i = 0; i < N; i++ ) {
      atomLines[i].position = positions[i];
    }
  }
  
  std::vector<LineInfo> lookupLineInfo;
  std::vector<PDBAtomLine> atomLines;
  std::vector<ConectRecord> conectRecords;
};

#endif /*PDB_H_*/
