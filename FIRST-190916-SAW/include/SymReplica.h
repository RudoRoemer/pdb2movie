// Author name:  Craig Jolley
// Created:      04 Jan 2008

#ifndef SYM_REPLICA_H
#define SYM_REPLICA_H

#include <vector>
#include <string>
#include "Vec3.h"

using namespace std;

typedef vector< vector<double> > Matrix;

class SymReplica {
private:
  struct Atom {
    string first30; // first 30 characters of ATOM entry; contains atom identifiers
    string last26;  // last 26 characters of ATOM entry; contains occupancy, etc.
    Vec3 pos;
  };
  string fileName;
  vector<Matrix> matrices;
  vector<Atom> atoms;
  void copyAtom(Atom &oldAtom, Atom &newAtom, Matrix m); 
  // transforms and copies an indivdual atom
  void readBioMT(); 
  // reads transformation matrices from file
  void readPDB(); 
  // reads atom positions from PDB file
  void showBioMT(Matrix m);
  void writePDB();
  // writes a PDB file containing all symmetry replicas
  char uniqueChain(string exclude); 
  // returns a unique valid chain identifier, excluding the ones in string 
public:
  SymReplica(string fname);
  void getMatrices(vector<Matrix> &exportMatrices);
};

#endif
