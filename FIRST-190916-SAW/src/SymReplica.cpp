// Symmetry.cpp -- Routines to create symmetry replicas of a PDB file based on 
//                 BIOMT matrices.  These procedures will only be used as a
//                 front-end before FIRST any analysis.
// Craig Jolley
// January 2008
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include "SymReplica.h"
#include "Parameters.h"
#include "Vec3.h"

extern Parameters parameters;

using namespace std;

SymReplica::SymReplica(string fname) {
  fileName = fname;
  readPDB();
  readBioMT();
  writePDB();
  return;
}

void SymReplica::copyAtom(SymReplica::Atom &oldAtom, SymReplica::Atom &newAtom, Matrix m) {
  newAtom.first30 = oldAtom.first30;
  newAtom.last26 = oldAtom.last26; 
  newAtom.pos.x = m[0][0]*oldAtom.pos.x + m[0][1]*oldAtom.pos.y + m[0][2]*oldAtom.pos.z + m[0][3];
  newAtom.pos.y = m[1][0]*oldAtom.pos.x + m[1][1]*oldAtom.pos.y + m[1][2]*oldAtom.pos.z + m[1][3];
  newAtom.pos.z = m[2][0]*oldAtom.pos.x + m[2][1]*oldAtom.pos.y + m[2][2]*oldAtom.pos.z + m[2][3]; 
}

void SymReplica::readBioMT() {
  ifstream inFile;
  inFile.open(fileName.c_str());
  if (!inFile.is_open()) {
    cerr << "Could not open " << fileName << " to read BioMT matrices!\n";
    exit(EXIT_FAILURE);
  }
  while (inFile.good()) {
    string input;  
    getline(inFile,input);
    if (input.substr(0,19) == "REMARK 350   BIOMT1") {
      vector <double>  nullVec(4);
      Matrix tempMatrix(3,nullVec);
      for (int i = 0; i < 3; i++) {
        tempMatrix[0][i] = atof(input.substr(24+i*10,9).c_str());
      }
      tempMatrix[0][3] = atof(input.substr(61,7).c_str());
      getline(inFile,input); // get second line
      if (input.substr(0,19) != "REMARK 350   BIOMT2") {
        cerr << "ERROR: BIOMT1 should be followed by BIOMT2! Skipping matrix.\n";
      } else {
        for (int i = 0; i < 3; i++) {
          tempMatrix[1][i] = atof(input.substr(24+i*10,9).c_str());
        }
        tempMatrix[1][3] = atof(input.substr(61,7).c_str());
        getline(inFile,input); // get third line
        if (input.substr(0,19) != "REMARK 350   BIOMT3") {
          cerr << "ERROR: BIOMT2 should be followed by BIOMT3! Skipping matrix.\n";
        } else {
          for (int i = 0; i < 3; i++) {
            tempMatrix[2][i] = atof(input.substr(24+i*10,9).c_str());
          }
          tempMatrix[2][3] = atof(input.substr(60,8).c_str());
          // tempMatrix is now complete; save it
          matrices.push_back(tempMatrix);
        }
      }
    }
  }
  if (!inFile.eof() && inFile.fail()) {
    cerr << "Error in reading file " << fileName << "; BIOMT matrices may not be complete.\n";
  }
  inFile.close();
  cout << "Read " << matrices.size() << " matrices from " << fileName << endl;
  // if matrices[0] isn't the identity matrix, there will be trouble
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      if (i == j && matrices[0][i][j] != 1.0) {
        cerr << "ERROR: First BIOMT entry should be the identity matrix!\n";
        exit(EXIT_FAILURE);
      } else if (i != j && matrices[0][i][j] != 0.0) {
        cerr << "ERROR: First BIOMT entry should be the identity matrix!\n";
        exit(EXIT_FAILURE);
      }        
    }    
  }
  return;
}

void SymReplica::readPDB() {
  ifstream inFile;
  inFile.open(fileName.c_str());
  if (!inFile.is_open()) {
    cerr << "Could not open " << fileName << " to read atoms.\n";
    exit(EXIT_FAILURE);
  }
  while (inFile.good()) {
    string input;  
    getline(inFile,input);
    if (input.substr(0,6) == "ATOM  " || input.substr(0,6) == "HETATM") {
      Atom tempAtom;
      tempAtom.first30 = input.substr(0,30);
      tempAtom.last26 = input.substr(54,26);
      tempAtom.pos.x = atof(input.substr(30,8).c_str());
      tempAtom.pos.y = atof(input.substr(38,8).c_str());
      tempAtom.pos.z = atof(input.substr(46,8).c_str());
      atoms.push_back(tempAtom);
    }
  }
  if (!inFile.eof() && inFile.fail()) {
    cerr << "Error in reading file " << fileName << " for symmetry-building; structure may not be complete.\n";
  }
  inFile.close();
  cout << "Read " << atoms.size() << " atoms from " << fileName << endl;
  return;
}

void SymReplica::showBioMT(Matrix m) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      cout << setw(8) << m[i][j];
    }
    cout << endl;
  }
  return;
}

char SymReplica::uniqueChain(string exclude) {
  string validChains = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890`-=[]\\;',./~!@#$%^&*()_+{}|:\"<>?";
  for (unsigned int i = 0; i < validChains.size(); i++) {
    char thisChain = validChains[i];
    if (exclude.find(thisChain,0) == string::npos) { 
      return thisChain;
    }
  }
  return '\n'; // not a valid PDB chain identifier
}

void SymReplica::writePDB() {
  ofstream outFile;
  string outFileName = fileName;
  outFileName.insert(outFileName.find(".pdb",0),"_sym");
  outFile.open(outFileName.c_str());
  if (!outFile.is_open()) {
    cerr << "writePDB() could not open " << outFileName << " for writing!\n";
    exit(EXIT_FAILURE);
  }
  // get list of chains used in PDB file
  string chains;
  for (unsigned int i = 0; i < atoms.size(); i++) {
    char c = atoms[i].first30[21];
    if (c == ' ') {
      c = atoms[i].first30[21] = 'X';
    }
    if (chains.find(c,0) == string::npos) { 
      chains += c;
    }
  }
  // replica 0 will use the original set of chain names; unique names
  // will be found for the chains in subsequent replicas
  // cout << "Chains contained in " << fileName << ": " << chains << endl;
  char thisChain = uniqueChain(chains);
  for (unsigned int j = 0; j < matrices.size(); j++) {
    //cout << "Writing replica " << j << endl;
    for (unsigned int i = 0; i < atoms.size(); i++) {
      Atom newAtom;
      copyAtom(atoms[i],newAtom,matrices[j]);
      // is a new chain ID needed?
      if (j > 0) {
        if (i == 0 || atoms[i].first30[21] != atoms[i-1].first30[21]) {
          outFile << "TER\n";
          thisChain = uniqueChain(chains); // get a new unique chain ID
          chains += thisChain;
          if (thisChain == '\n') {
            cout << "ERROR: Too many chains!  Chain ID's may not be unique.\n";
            thisChain = atoms[i].first30[21];
          }
        } 
        newAtom.first30[21] = thisChain;
      }
      // fix atom number
      int newAtomNum = j*atoms.size() + i + 1;
      outFile << newAtom.first30.substr(0,6);
      outFile.setf(ios_base::right, ios_base::adjustfield);
      if (newAtomNum < 100000) {
        outFile << setw(5) << newAtomNum;
      } else {
        outFile << "*****";
      }
      outFile << newAtom.first30.substr(11,19);
      // output new coordinates
      outFile.setf(ios_base::fixed);
      outFile.precision(3);
      outFile.setf(ios_base::showpoint);
      outFile.setf(ios_base::right, ios_base::adjustfield);
      outFile << setw(8) << newAtom.pos.x;     
      outFile << setw(8) << newAtom.pos.y;
      outFile << setw(8) << newAtom.pos.z;
      outFile << newAtom.last26 << endl;
    }
  }
  outFile << "END\n";
  outFile.close();
  return;
}

void SymReplica::getMatrices(vector<Matrix> &exportMatrices) {
  exportMatrices = matrices;
  return;
}


    
