#include <cmath>
#include "global_defs.h"
#include "SAXS.h"
#include "Froda.h"

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Initializes the electron-density lookup tables contained in EDStructure.
//   Using the lookup tables is generally much faster than calculating these
//   values each time they're needed; this initialization function is called
//   by both versions of the EDStructure constructor.
//   The values used to initialize ai[][] and bi[][] are taken from
//   Rez et al., Acta Cryst. A50, 481-497 (1994) pg. 493
////////////////////////////////////////////////////////////////////////////////
SAXS::LookupTable::LookupTable(double resolutionFactor) {
  binSize = 0.01;
  //cerr << "Number of bins: " << numBins << endl;
  allZero[0] = allZero[1] = allZero[2] = allZero[3] = 0.0;
  // initialize carbon array
  ai[carbon][0] = 1.7401;    bi[carbon][0] = 38.1927;
  ai[carbon][1] = 2.4357;    bi[carbon][1] = 13.7125;
  ai[carbon][2] = 1.3245;    bi[carbon][2] = 0.6520;
  ai[carbon][3] = 0.4960;    bi[carbon][3] = 0.1488;
  // initialize nitrogen array
  ai[nitrogen][0] = 0.4743;    bi[nitrogen][0] = 0.1041;
  ai[nitrogen][1] = 2.9082;    bi[nitrogen][1] = 9.1890;
  ai[nitrogen][2] = 2.2778;    bi[nitrogen][2] = 27.0869;
  ai[nitrogen][3] = 1.3332;    bi[nitrogen][3] = 0.4612; 
  // initialize oxygen array
  ai[oxygen][0] = 3.5531;    bi[oxygen][0] = 6.8702;
  ai[oxygen][1] = 2.6162;    bi[oxygen][1] = 21.0743;
  ai[oxygen][2] = 1.2120;    bi[oxygen][2] = 0.3871;
  ai[oxygen][3] = 0.6107;    bi[oxygen][3] = 0.0960;  
  // initialize sulfur array
  ai[sulfur][0] = 6.5881;    bi[sulfur][0] = 27.6154;
  ai[sulfur][1] = 6.7734;    bi[sulfur][1] = 1.5353;
  ai[sulfur][2] = 1.1051;    bi[sulfur][2] = 0.4580;
  ai[sulfur][3] = 1.4747;    bi[sulfur][3] = 0.0452;
  // initialize phosphorus array
  ai[phosphorus][0] = 5.5289;    bi[phosphorus][0] = 32.2841;
  ai[phosphorus][1] = 6.9083;    bi[phosphorus][1] = 1.8577;
  ai[phosphorus][2] = 1.0702;    bi[phosphorus][2] = 0.4743;
  ai[phosphorus][3] = 1.4395;    bi[phosphorus][3] = 0.0500;
   // initialize iron 2+ array
  ai[iron][0] = 9.2568;    bi[iron][0] = 4.3145;
  ai[iron][1] = 6.2209;    bi[iron][1] = 10.8026;
  ai[iron][2] = 6.8432;    bi[iron][2] = 0.3454;
  ai[iron][3] = 1.6267;    bi[iron][3] = 0.0200;
  // initialize magnesium 2+ array
  ai[magnesium][0] = 1.1002;    bi[magnesium][0] = 16.6083;
  ai[magnesium][1] = 6.6744;    bi[magnesium][1] = 3.2353;
  ai[magnesium][2] = 0.9324;    bi[magnesium][2] = 0.3778;
  ai[magnesium][3] = 1.2525;    bi[magnesium][3] = 0.0686;
  // set up array for sin(x)/x -- max value is for x = 100
  // because sin(x)/x approaches zero
  sinXOverX = new double[10001];
  sinXOverX[0] = 1.0;
  for (int bin = 1; bin < 10001; bin++) {
    double x = bin*binSize;
    sinXOverX[bin] = sin(x) / x;
  }
  // set up array for exp(-x) -- max value is for x = 25, 
  // exp(-x) = 10^(-11); close enough to zero
  eToTheMinusX = new double[2501];
  for (int bin = 0; bin < 2501; bin++) {
    double x = bin*binSize;
    eToTheMinusX[bin] = exp(-x);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for the LookupTable class.
////////////////////////////////////////////////////////////////////////////////
SAXS::LookupTable::~LookupTable(){
  delete [] ai;
  delete [] bi;
  delete [] sinXOverX;
  delete [] eToTheMinusX;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the square of the distance between a point (x2,y2,z2) in Cartesian
//   space and the location of the atom.
////////////////////////////////////////////////////////////////////////////////
double SAXS::Atom::distanceSquared(const Atom &a2) {
  double temp;
  temp = (*x - *(a2.x))*(*x - *(a2.x)) + (*y - *(a2.y))*(*y - *(a2.y)) + 
    (*z - *(a2.z))*(*z - *(a2.z));
  return temp;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for the SAXS class.  This version of the contstructor
//   is to be used when creating an SAXS object that will be associated 
//   with conformers produced by a FRODA simulation.
// Parameters:
//   Froda &froda -- Froda object associated with the EDStructure object being
//           created.  Pointers to atom positions in the EDStructure class are 
//           assigned to data in this Froda object
//   Mol_framework &structInput -- This will contain (among other things) the
//           element names of the Atom objects contained in the 
//           EDStructure::atoms vector
//   double resolutionFactorInput -- resolution factor used by 
//           EDStructure::densityAt() 
///////////////////////////////////////////////////////////////////////////////
SAXS::SAXS(Froda &froda, MolFramework &structureInput, 
  double resolutionFactorInput, string filename) {
  using namespace std;
  resolutionFactor = resolutionFactorInput;
  lookupTable = new LookupTable(resolutionFactor);
  structure = &structureInput;
  for (unsigned int siteNum = 1; siteNum <= structure->total_sites; siteNum++) {
    string newElement = structure->getElementName(siteNum);
    Atom * newAtom = new Atom();
    if (newElement=="C ") {
      newAtom->ai = lookupTable->ai[carbon];
      newAtom->bi = lookupTable->bi[carbon];
    }
    else if (newElement=="N ") {
      newAtom->ai = lookupTable->ai[nitrogen];
      newAtom->bi = lookupTable->bi[nitrogen];
    }
    else if (newElement=="O ") {
      newAtom->ai = lookupTable->ai[oxygen];
      newAtom->bi = lookupTable->bi[oxygen];
    }
    else if (newElement=="P ") {
      newAtom->ai = lookupTable->ai[phosphorus];
      newAtom->bi = lookupTable->bi[phosphorus];
    }
    else if (newElement=="S ") {
      newAtom->ai = lookupTable->ai[sulfur];
      newAtom->bi = lookupTable->bi[sulfur];
    }
    else if (newElement=="MG") { // magnesium
      newAtom->ai = lookupTable->ai[magnesium];
      newAtom->bi = lookupTable->bi[magnesium];
    }
    else if (newElement=="FE") { // iron
      newAtom->ai = lookupTable->ai[iron];
      newAtom->bi = lookupTable->bi[iron];
    }
    else if (newElement =="H ") {
      newAtom->ai = lookupTable->allZero;
      newAtom->bi = lookupTable->allZero;
    }
    else {
      cerr << "Unknown element " << newElement << endl;
      newAtom->ai = lookupTable->allZero;
      newAtom->bi = lookupTable->allZero;
    }
    newAtom->x = &(froda.currentPos.at(siteNum).x);
    newAtom->y = &(froda.currentPos.at(siteNum).y);
    newAtom->z = &(froda.currentPos.at(siteNum).z);
    atoms.push_back(newAtom);
    // the index of the new atom in atoms[] will be siteNum-1;
  }
  // load SAXS profile from file
  if (filename != "") { // if filename empty, don't load anything
    ifstream inFile;
    inFile.open(filename.c_str());
    if (!inFile.is_open()) {
      cerr << "Failed to open file " << filename << " for reading.\n";
      exit(EXIT_FAILURE);
    }
    while (inFile.good()) {
      SAXSPoint tempPoint;
      inFile >> tempPoint.k;
      //cout << "Read '" << tempPoint.k << "' from file.\n";
      inFile >> tempPoint.intensity;
      //cout << "Read '" << tempPoint.intensity << "' from file.\n";
      //CJ TODO: My I/O checking here could be more elegant
      if (tempPoint.intensity != 0) {
        saxsTarget.push_back(tempPoint);
        tempPoint.intensity = 0;
        saxsCurrent.push_back(tempPoint); // set k entries in current array
      }
    }
    inFile.close();
  }
  writeSAXS(structure->base_name + "_initial.saxs"); 
  for (unsigned int i = 0; i < saxsCurrent.size(); i++) {
    //cout << "saxsCurrent[ " << i << "] = " << saxsCurrent[i].k << ',' << saxsCurrent[i].intensity;
    cout << "saxsTarget[ " << i << "] = " << saxsTarget[i].k << ',' << saxsTarget[i].intensity << endl;
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for the SAXS class.  Frees up the memory allocated by
//   constructor.
///////////////////////////////////////////////////////////////////////////////
SAXS::~SAXS() {
  for (unsigned int i = 0; i < atoms.size(); i++) {
    delete atoms[i];
  }
  delete lookupTable;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the SAXS scattering intensity returned by the structure at
//   wavenumber k. 
//   See Craig's lab notebook for 10/31/06, 12/08/06 for details on math.
///////////////////////////////////////////////////////////////////////////////
double SAXS::saxs(double k) {
  double result = 0;
  int n = atoms.size();
  for (int i = 0; i < n; i++) { // first sum off-diagonal
    for (int j = 0; j < i; j++) {
      double kr = k * sqrt(atoms[i]->distanceSquared(*atoms[j]));
      if (kr < 100) { // maximum value at which lookup table is defined
        int bin2 = (int) (kr / lookupTable->binSize);
        for (int a = 0; a < 4; a++) {
          for (int b = 0; b < 4; b++) {
            /*result += 2 * atoms[i]->ai[a] * atoms[j]->ai[b] 
              * exp(-k*k*(4*PI*PI*resolutionFactor + atoms[i]->bi[a] + atoms[j]->bi[b])/(4*PI*PI))  
	      * lookupTable->sinXOverX[bin2];*/
            double x = k*k*(4*PI*PI*resolutionFactor + atoms[i]->bi[a] + atoms[j]->bi[b])/(4*PI*PI);
            if (x < 25) {
              int bin1 = (int) (x / lookupTable->binSize); 
              result += 2 * atoms[i]->ai[a] * atoms[j]->ai[b] 
                * lookupTable->eToTheMinusX[bin1] * lookupTable->sinXOverX[bin2];
            }
          }
        }
      }
    }
  }
  for (int i = 0; i < n; i++) { // now sum diagonal; i=j
    for (int a = 0; a < 4; a++) {
      for (int b = 0; b < 4; b++) {
        /*result += atoms[i]->ai[a] * atoms[i]->ai[b] 
	 * exp(-k*k*(4*PI*PI*resolutionFactor + atoms[i]->bi[a] + atoms[i]->bi[b])/(4*PI*PI));*/
        double x = k*k*(4*PI*PI*resolutionFactor + atoms[i]->bi[a] + atoms[i]->bi[b])/(4*PI*PI);
        if (x < 25) {
          int bin = (int) (x / lookupTable->binSize);
          result += atoms[i]->ai[a] * atoms[i]->ai[b] 
            * lookupTable->eToTheMinusX[bin];
        }
      }
    }
  }
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Similar to SAXS::saxs, except that only the contribution to the 
//   scattering from pairs of atoms with a specified degree of separation in 
//   the bond network is included.  Returns only the diagonal term if 
//   degree = 0.
////////////////////////////////////////////////////////////////////////////////
double SAXS::restrictedSAXS(double k, int degree) {
  double result = 0;
  int n = atoms.size();
  if (degree == 0) { // sum diagonal only
    for (int i = 0; i < n; i++) {
      for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
          double x = k*k*(4*PI*PI*resolutionFactor + atoms[i]->bi[a] + atoms[i]->bi[b])/(4*PI*PI);
          if (x < 25) {
            int bin = (int) (x / lookupTable->binSize);
            result += atoms[i]->ai[a] * atoms[i]->ai[b] 
              * lookupTable->eToTheMinusX[bin];
          }
        }
      }
    }
  } else { // looking for off-diagonal terms
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        if (structure->NthNearestNeighbor(i+1,j+1,degree) && 
            !structure->NthNearestNeighbor(i+1,j+1,degree-1)) {
      // element 0 in the vector corresponds to atom 1 in FIRST's internal 
      // numbering; hence the +1 in the arguments to NthNearestNeighbor
          double kr = k * sqrt(atoms[i]->distanceSquared(*atoms[j]));
          if (kr < 100) { // maximum value at which lookup table is defined
            int bin2 = (int) (kr / lookupTable->binSize);
            for (int a = 0; a < 4; a++) {
              for (int b = 0; b < 4; b++) {
                double x = k*k*(4*PI*PI*resolutionFactor + atoms[i]->bi[a] + atoms[j]->bi[b])/(4*PI*PI);
                if (x < 25) {
                  int bin1 = (int) (x / lookupTable->binSize); 
                  result += 2 * atoms[i]->ai[a] * atoms[j]->ai[b] 
                    * lookupTable->eToTheMinusX[bin1] * lookupTable->sinXOverX[bin2];
                }
              }
            }
          }
        }
      }
    }
  }
  return result;
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Output the SAXS profile stored in currentSAXS to the specified file
////////////////////////////////////////////////////////////////////////////////
void SAXS::writeSAXS(string filename) {
  ofstream outFile;
  outFile.open(filename.c_str());
  if (!outFile.is_open()) {
    cerr << "Failed to open file " << filename << " for writing.\n";
    exit(EXIT_FAILURE);
  }
  cout << "Writing SAXS profile to " << filename << "...    0%";
  //for (unsigned int i = 0; i < saxsCurrent.size(); i++) {
  //  outFile << saxsCurrent[i].k << ' ' << saxsCurrent[i].intensity << endl;
  //}
  double k = 0;
  double maxk = 0.5;
  double step = 0.01;
  double saxsZero = saxs(0);
  outFile << "0 1.0\n";
  k += step;
  while (k <= maxk) {
    outFile << k << ' ' << saxs(k)/saxsZero << endl;
    k += step;
    cout << "\b\b\b\b\b" << setw(4) << setprecision(3) << k*100/maxk << '%' << flush;
  }
  cout << endl;
  outFile.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Uses the current conformation contained in *structure to update 
//   currentSAXS
////////////////////////////////////////////////////////////////////////////////
void SAXS::updateProfile() {
  double saxsZero = saxs(0);
  for (unsigned int i = 0; i < saxsCurrent.size(); i++) {
    saxsCurrent[i].intensity = saxs(saxsCurrent[i].k) / saxsZero;
    //cout << saxsCurrent[i].k << ' ' << saxsCurrent[i].intensity << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//    Returns the correlation between saxsCurrent[] and saxsTarget[].  The 
//    correlation function can be altered as needed to give better performance.
////////////////////////////////////////////////////////////////////////////////
double SAXS::correlate() {
  updateProfile();
  if (saxsCurrent.size() == 1) {
    return 1.0 / abs(saxsCurrent[0].intensity - saxsTarget[0].intensity);
  }
  double cNumerator = 0.0;     // numerator of C
  double cDenominator1 = 0.0;  // first sum in denominator of C
  double cDenominator2 = 0.0;  // second sum in denominator of C
  if (saxsCurrent.size() != saxsTarget.size()) {
    cerr << "ERROR: saxsCurrent has " << saxsCurrent.size()  << "elements,";
    cerr << " saxsTarget has " << saxsTarget.size() << "!\n";
    exit(EXIT_FAILURE);
  }
  for (unsigned int i = 0; i < saxsCurrent.size(); i++) {
    if (saxsCurrent[i].k != saxsTarget[i].k) {
      cerr << "ERROR: Discrepancy in k-column of SAXS profiles at position ";
      cerr << i << endl;
      exit(EXIT_FAILURE);
    }
    //double k4 = saxsCurrent[i].k*saxsCurrent[i].k*saxsCurrent[i].k*saxsCurrent[i].k;
    cNumerator += saxsCurrent[i].intensity * saxsTarget[i].intensity;
    cDenominator1 +=  saxsCurrent[i].intensity * saxsCurrent[i].intensity;
    cDenominator2 +=  saxsTarget[i].intensity * saxsTarget[i].intensity;
  }  
  //cout << "components = " << cNumerator << ',' << cDenominator1 << ',' << cDenominator2 << endl;
  double result = cNumerator/sqrt(cDenominator1*cDenominator2);
  //cout << "result = " << result << endl;
  return result;
}
