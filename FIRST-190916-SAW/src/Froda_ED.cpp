// Author name: Craig Jolley
// Created:     June 2006

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <cmath>
#include <vector>
#include <ctime>
#include "GenericMap.h"
#include "EDStructure.h"
#include "MolFramework.h"
#include "Froda.h"
#include "Grid.h"
#include "global_defs.h"
#include "mt19937ar.h"
#include "Parameters.h"

extern const Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Returns the square of the distance between a point (x2,y2,z2) in Cartesian
//   space and the location of the atom.
////////////////////////////////////////////////////////////////////////////////
double EDStructure::Atom::distanceSquared(double x2, double y2, double z2) {
  double temp;
  temp = (*x - x2)*(*x - x2) + (*y - y2)*(*y - y2) + (*z - z2)*(*z - z2);
  return temp;
}


////////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the electron density contributed by an atom at a distance r 
// Parameters:
//   double r2 -- the square of the distance from the atom
//   double resolutionFactor -- the resolution factor by which the atomic
//                              Gaussians should be broadened
// Return Value List:
//   Returns the calculated contribution to the electron density at r.
////////////////////////////////////////////////////////////////////////////////
double EDStructure::Atom::density(double r2, double binSize) {
  double rho;
  if (lookup==NULL) {
    rho = 0.0;
  } else {
    int bin = (int) (r2 / binSize);
    rho = lookup[bin];
  }
  return rho;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Initializes the electron-density lookup tables contained in EDStructure.
//   Using the lookup tables is generally much faster than calculating these
//   values each time they're needed; this initialization function is called
//   by both versions of the EDStructure constructor.
//   The values used to initialize ai[][] and bi[][] are taken from
//   Rez et al., Acta Cryst. A50, 481-497 (1994) pg. 493
////////////////////////////////////////////////////////////////////////////////
EDStructure::LookupTable::LookupTable(double resolutionFactor, double &cutoff) {
  binSize = 0.01;
  cutoff = 3.6 * sqrt(resolutionFactor) + 3.0;
  //cerr << "cutoff = " << cutoff << endl;
  int numBins = 2 + (int) (cutoff*cutoff/binSize);
  //cerr << "Number of bins: " << numBins << endl;
  cDensity = new double[numBins];
  nDensity = new double[numBins];
  oDensity = new double[numBins];
  sDensity = new double[numBins];
  pDensity = new double[numBins];
  feDensity = new double[numBins];
  mgDensity = new double[numBins];
  double a[4];
  double b[4];
  allZero[0] = allZero[1] = allZero[2] = allZero[3] = 0.0;
  // initialize carbon array
  ai[carbon][0] = 1.7401;    bi[carbon][0] = 38.1927;
  ai[carbon][1] = 2.4357;    bi[carbon][1] = 13.7125;
  ai[carbon][2] = 1.3245;    bi[carbon][2] = 0.6520;
  ai[carbon][3] = 0.4960;    bi[carbon][3] = 0.1488;
  for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[carbon][j]*pow(PI/(bi[carbon][j]+8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[carbon][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    cDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      cDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  }
  // initialize nitrogen array
  ai[nitrogen][0] = 0.4743;    bi[nitrogen][0] = 0.1041;
  ai[nitrogen][1] = 2.9082;    bi[nitrogen][1] = 9.1890;
  ai[nitrogen][2] = 2.2778;    bi[nitrogen][2] = 27.0869;
  ai[nitrogen][3] = 1.3332;    bi[nitrogen][3] = 0.4612; 
    for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[nitrogen][j]*pow(PI/(bi[nitrogen][j] + 8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[nitrogen][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    nDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      nDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  }
  // initialize oxygen array
  ai[oxygen][0] = 3.5531;    bi[oxygen][0] = 6.8702;
  ai[oxygen][1] = 2.6162;    bi[oxygen][1] = 21.0743;
  ai[oxygen][2] = 1.2120;    bi[oxygen][2] = 0.3871;
  ai[oxygen][3] = 0.6107;    bi[oxygen][3] = 0.0960;  
  for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[oxygen][j]*pow(PI/(bi[oxygen][j] + 8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[oxygen][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    oDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      oDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  }
  // initialize sulfur array
  ai[sulfur][0] = 6.5881;    bi[sulfur][0] = 27.6154;
  ai[sulfur][1] = 6.7734;    bi[sulfur][1] = 1.5353;
  ai[sulfur][2] = 1.1051;    bi[sulfur][2] = 0.4580;
  ai[sulfur][3] = 1.4747;    bi[sulfur][3] = 0.0452;
  for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[sulfur][j]*pow(PI/(bi[sulfur][j] + 8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[sulfur][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    sDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      sDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  }

  // initialize phosphorus array
  ai[phosphorus][0] = 5.5289;    bi[phosphorus][0] = 32.2841;
  ai[phosphorus][1] = 6.9083;    bi[phosphorus][1] = 1.8577;
  ai[phosphorus][2] = 1.0702;    bi[phosphorus][2] = 0.4743;
  ai[phosphorus][3] = 1.4395;    bi[phosphorus][3] = 0.0500;
  for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[phosphorus][j]*pow(PI/(bi[phosphorus][j] + 8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[phosphorus][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    pDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      pDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  } 
   // initialize iron 2+ array
  ai[iron][0] = 9.2568;    bi[iron][0] = 4.3145;
  ai[iron][1] = 6.2209;    bi[iron][1] = 10.8026;
  ai[iron][2] = 6.8432;    bi[iron][2] = 0.3454;
  ai[iron][3] = 1.6267;    bi[iron][3] = 0.0200;
  for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[iron][j]*pow(PI/(bi[iron][j] + 8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[iron][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    pDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      pDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  }
  // initialize magnesium 2+ array
  ai[magnesium][0] = 1.1002;    bi[magnesium][0] = 16.6083;
  ai[magnesium][1] = 6.6744;    bi[magnesium][1] = 3.2353;
  ai[magnesium][2] = 0.9324;    bi[magnesium][2] = 0.3778;
  ai[magnesium][3] = 1.2525;    bi[magnesium][3] = 0.0686;
  for (int j = 0; j < 4; j++) {
    a[j] = 8*ai[magnesium][j]*pow(PI/(bi[magnesium][j] + 8*PI*PI*resolutionFactor),1.5);
    b[j] = bi[magnesium][j] + 8*PI*PI*resolutionFactor;
  }
  for (int bin = 0; bin < numBins; bin++) {
    double r2 = bin*binSize;
    pDensity[bin] = 0.0;
    for (int i = 0; i < 4; i++) {
      pDensity[bin] += a[i]*exp(-4*PI*PI*r2/b[i]);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for the LookupTable class.
////////////////////////////////////////////////////////////////////////////////
EDStructure::LookupTable::~LookupTable(){
  delete [] cDensity;
  delete [] nDensity;
  delete [] oDensity;
  delete [] sDensity;
  delete [] pDensity;
  delete [] feDensity;
  delete [] mgDensity;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for the EDStructure class.  Loads in a PDB file, paying 
//   attention only to the element names and coordinates contained in the ATOM 
//   and HETATM records.  At the moment, this constructor is used by the 
//   constructor of the TheoMap class, which constructs a theoretical electron
//   density map based on the data contained in the PDB file.  In this version,
//   I'm not bothering with setting up the coarse grid, since an EDStructure
//   created by this constructor only calculates its electron density once;
//   so it's not so tragic if it scales badly. 
// Parameters:
//   string fileName -- name of the PDB file to be loaded
//   double resolutionFactorInput -- value of the resolution factor used in
//                                   Atom::density()
////////////////////////////////////////////////////////////////////////////////
EDStructure::EDStructure(string fileName, double resolutionFactorInput) {
  using namespace std;
  resolutionFactor = resolutionFactorInput;
  lookupTable = new LookupTable(resolutionFactor,cutoff);
  // loading from file, so this object won't be connected to any MolFramework 
  // object; set structure to null pointer
  structure = NULL;

  if( parameters.verbose ){
    cout << "Loading " << fileName << "...\n";
  }
  ifstream inFile;
  inFile.open(fileName.c_str());
  if (!inFile.is_open()) {
    cerr << "Failed to open file " << fileName << " for reading.\n";
    exit(EXIT_FAILURE);
  }
  string nextLine;
  string knownElements = "HCNOPS"; // one-letter element names
  while(inFile.good()) {
    getline(inFile,nextLine); // read the next line
    if ((nextLine.substr(0,4)=="ATOM")||(nextLine.substr(0,6)=="HETATM")) {
      string newElement;
      // positions given correspond to column assignments in the 
      // official PDB format
      // check two-letter element names first
      if (nextLine.substr(12,2) == "MG") { newElement = "MG"; }
      else if (nextLine.substr(12,2) == "FE") { newElement = "FE"; }
      else if ((knownElements.find(nextLine[13],0) == string::npos)&&(nextLine.length() > 77)) {
        newElement = nextLine[77];
      }
      else {
        newElement = nextLine[13];
      }
      if (newElement.length() == 1) {
        newElement += " ";
      }
      Atom * newAtom = new Atom();
      if (newElement=="C ") {
        newAtom->lookup = lookupTable->cDensity;
        newAtom->ai = lookupTable->ai[carbon];
        newAtom->bi = lookupTable->bi[carbon];
      }
      else if (newElement=="N ") {
        newAtom->lookup = lookupTable->nDensity;
        newAtom->ai = lookupTable->ai[nitrogen];
        newAtom->bi = lookupTable->bi[nitrogen];
      }
      else if (newElement=="O ") {
        newAtom->lookup = lookupTable->oDensity;
        newAtom->ai = lookupTable->ai[oxygen];
        newAtom->bi = lookupTable->bi[oxygen];
      }
      else if (newElement=="P ") {
        newAtom->lookup = lookupTable->pDensity;
        newAtom->ai = lookupTable->ai[phosphorus];
        newAtom->bi = lookupTable->bi[phosphorus];
      }
      else if (newElement=="S ") {
        newAtom->lookup = lookupTable->sDensity;
        newAtom->ai = lookupTable->ai[sulfur];
        newAtom->bi = lookupTable->bi[sulfur];
      }
      else if (newElement=="MG") { // magnesium
        newAtom->lookup = lookupTable->mgDensity;
        newAtom->ai = lookupTable->ai[magnesium];
        newAtom->bi = lookupTable->bi[magnesium];
      }
      else if (newElement=="FE") { // iron
        newAtom->lookup = lookupTable->feDensity;
        newAtom->ai = lookupTable->ai[iron];
        newAtom->bi = lookupTable->bi[iron];
      }
      else if (newElement =="H ") {
        newAtom->lookup = NULL;
        newAtom->ai = lookupTable->allZero;
        newAtom->bi = lookupTable->allZero;
      }
      else {
        cerr << "Unknown element " << newElement << endl;
        newAtom->ai = lookupTable->allZero;
        newAtom->bi = lookupTable->allZero;
        newAtom->lookup = NULL;
      }
      if (newAtom->lookup != NULL) { // don't bother if atom contributes no density
        newAtom->x = new double (atof(nextLine.substr(30,37).c_str()));
        newAtom->y = new double (atof(nextLine.substr(38,45).c_str()));
        newAtom->z = new double (atof(nextLine.substr(46,53).c_str()));
        atoms.push_back(newAtom);  // add to vector
        if (atoms.size()==1) {
        // if this is the first atom loaded into memory, use its
        // coordinates as initial guesses for the max/min values
          max.x = *(newAtom->x);
          max.y = *(newAtom->y);
          max.z = *(newAtom->z);
          min.x = *(newAtom->x);
          min.y = *(newAtom->y);
          min.z = *(newAtom->z);
        } 
        else {
          // if not, see if the current atom can be used to improve
          // the max/min values
          if (*(newAtom->x) > max.x) {
            max.x = *(newAtom->x);
          }
          else if (*(newAtom->x) < min.x) {
            min.x = *(newAtom->x);
          }
          if (*(newAtom->y) > max.y) {
            max.y = *(newAtom->y);
          }
          else if (*(newAtom->y) < min.y) {
            min.y = *(newAtom->y);		
          }
          if (*(newAtom->z) > max.z) {
            max.z = *(newAtom->z);
          }
          else if (*(newAtom->z) < min.z) {
            min.z = *(newAtom->z);
          }
        }
      }	 
    }
  }
  cout << "Read " << atoms.size() << " atoms.\n";
  inFile.close();
  // not using coarseGrid in this version
  coarseGrid = NULL;
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for the EDStructure class.  This version of the contstructor
//   is to be used when creating an EDStructure object that will be associated 
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
EDStructure::EDStructure(Froda &froda, MolFramework &structureInput, 
                         GenericMap &gmap, double resolutionFactorInput) {
  using namespace std;
  // Set default value for cutoff;
  resolutionFactor = resolutionFactorInput;
  lookupTable = new LookupTable(resolutionFactor,cutoff);
  // assign pointers
  structure = &structureInput;
  max.x = structure->max_coords[0];
  max.y = structure->max_coords[1];
  max.z = structure->max_coords[2];
  min.x = structure->min_coords[0];
  min.y = structure->min_coords[1];
  min.z = structure->min_coords[2];
  for (unsigned int siteNum = 1; siteNum <= structure->total_sites; siteNum++) {
    string newElement = structure->getElementName(siteNum);
    Atom * newAtom = new Atom();
    if (newElement=="C ") {
      newAtom->lookup = lookupTable->cDensity;
      newAtom->ai = lookupTable->ai[carbon];
      newAtom->bi = lookupTable->bi[carbon];
    }
    else if (newElement=="N ") {
      newAtom->lookup = lookupTable->nDensity;
      newAtom->ai = lookupTable->ai[nitrogen];
      newAtom->bi = lookupTable->bi[nitrogen];
    }
    else if (newElement=="O ") {
      newAtom->lookup = lookupTable->oDensity;
      newAtom->ai = lookupTable->ai[oxygen];
      newAtom->bi = lookupTable->bi[oxygen];
    }
    else if (newElement=="P ") {
      newAtom->lookup = lookupTable->pDensity;
      newAtom->ai = lookupTable->ai[phosphorus];
      newAtom->bi = lookupTable->bi[phosphorus];
    }
    else if (newElement=="S ") {
      newAtom->lookup = lookupTable->sDensity;
      newAtom->ai = lookupTable->ai[sulfur];
      newAtom->bi = lookupTable->bi[sulfur];
    }
    else if (newElement=="MG") { // magnesium
      newAtom->lookup = lookupTable->sDensity;
      newAtom->ai = lookupTable->ai[magnesium];
      newAtom->bi = lookupTable->bi[magnesium];
    }
    else if (newElement=="FE") { // iron
      newAtom->lookup = lookupTable->sDensity;
      newAtom->ai = lookupTable->ai[iron];
      newAtom->bi = lookupTable->bi[iron];
    }
    else if (newElement =="H ") {
      newAtom->lookup = NULL;
      newAtom->ai = lookupTable->allZero;
      newAtom->bi = lookupTable->allZero;
    }
    else {
      cerr << "Unknown element " << newElement << endl;
      newAtom->lookup = NULL;
      newAtom->ai = lookupTable->allZero;
      newAtom->bi = lookupTable->allZero;
    }
    newAtom->x = &(froda.currentPos.at(siteNum).x);
    newAtom->y = &(froda.currentPos.at(siteNum).y);
    newAtom->z = &(froda.currentPos.at(siteNum).z);
    atoms.push_back(newAtom);
    // the index of the new atom in atoms[] will be siteNum-1;
  }
  // if GenericMap map is larger than the structure, extend the boundaries
  gridPoint origin;
  origin.i = 0;
  origin.j = 0;
  origin.k = 0;
  if (gmap.gridPointToVector(origin).x < structure->min_coords[X]) {
    structure->min_coords[X] = gmap.gridPointToVector(origin).x; 
  }
  if (gmap.gridPointToVector(origin).y < structure->min_coords[Y]) {
    structure->min_coords[Y] = gmap.gridPointToVector(origin).y; 
  }
  if (gmap.gridPointToVector(origin).z < structure->min_coords[Z]) {
    structure->min_coords[Z] = gmap.gridPointToVector(origin).z;
  }
  if (gmap.gridPointToVector(gmap.getExtent()).x > structure->max_coords[X]) {
    structure->max_coords[X] = gmap.gridPointToVector(gmap.getExtent()).x;
  }
  if (gmap.gridPointToVector(gmap.getExtent()).y > structure->max_coords[Y]) {
    structure->max_coords[Y] = gmap.gridPointToVector(gmap.getExtent()).y;
  }
  if (gmap.gridPointToVector(gmap.getExtent()).z > structure->max_coords[Z]) {
    structure->max_coords[Z] = gmap.gridPointToVector(gmap.getExtent()).z;
  }
  // now set up the coarse grid
  coarseGrid = new Grid(*structure, cutoff/sqrt(2.0));
  // FWIW, I'm a little confused by the fact that I need to use
  // cutoff/sqrt(2) -- from the way the grid routines are written
  // it looks like cutoff/2 should do the trick, but this seems
  // to leave out some atoms.
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Destructor for the EDStructure class.  Frees up the memory allocated by
//   one of the two constructor functions.
///////////////////////////////////////////////////////////////////////////////
EDStructure::~EDStructure() {

  for (unsigned int i = 0; i < atoms.size(); i++) {
    // if current_atom.structure points to NULL, then delete x, y, z in the 
    // current atom, because it was created using new double operators in the 
    // file-based constructor
    if (structure == NULL) { 
      delete atoms[i]->x;
      delete atoms[i]->y;
      delete atoms[i]->z;
    }
    delete atoms[i];
  }
  delete coarseGrid;
  delete lookupTable;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the electron density at a given point in space resulting from 
//   the atoms contained in an EDStructure object.
// Parameters:
//   Vector v -- the location in Cartesian space where the density is to be
//               calculated
// Return Value List:
//   Returns the calculated electron density
///////////////////////////////////////////////////////////////////////////////
double EDStructure::densityAt(Vector v) {
  double totalDensity = 0.0;
  if (coarseGrid==NULL) { // if coarse grid isn't being used, check every atom
    for (size_t i = 0; i < atoms.size(); i++) {
      double r2 = atoms[i]->distanceSquared(v.x,v.y,v.z);
      if (r2 < cutoff*cutoff) {
        totalDensity += atoms[i]->density(r2,lookupTable->binSize);
      } 
    }
  }
  else { // use the coarse grid to identify atoms that are near v
    std::vector<SiteID> nearbyAtoms;
    coarseGrid->getAtomsWithin2LengthsOfPoint(v,nearbyAtoms);
    // Next, identify the atoms that are within the cutoff radius.
    for (size_t i = 0; i < nearbyAtoms.size(); i++) {
      double r2 = atoms[nearbyAtoms[i]-1]->distanceSquared(v.x,v.y,v.z);
      // element 0 in the vector corresponds to atom 1 in FIRST's internal 
      // numbering; hence the -1 in the index of atoms[]
      if (r2 < cutoff*cutoff) {
        totalDensity += atoms[nearbyAtoms[i]-1]->density(r2,lookupTable->binSize);
      }  
    }
  }
  //cerr << totalDensity << endl;
  return totalDensity;
}


////////////////////////////////////////////////////////////////////////////////
// Description:  
//   Constructor for GenericMap class.  Currently doesn't do anything.
////////////////////////////////////////////////////////////////////////////////
GenericMap::GenericMap(){
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
GenericMap::~GenericMap(){
  delete [] gridData;
}


///////////////////////////////////////////////////////////////////////////////
// Description:
//   Given crystallographic grid indices i,j,k, the location of the 
//   corresponding point in the 1-D array GenericMap::gridData is given by:
//   n = i + extent.i*j + extent.i*extent.j*k
//   This function does the reverse transformation
// Parameters:
//   int n -- a position in the 1-D array GenericMap::gridData
// Return Value List:
//   Returns a gridPoint structure containing the i,j,k indices corresponding 
//   to the 1-D index n
///////////////////////////////////////////////////////////////////////////////
gridPoint GenericMap::gridIndex(int n) {
  gridPoint temp;
  temp.i = n % extent.i;
  temp.j = (n % (extent.i*extent.j) - temp.i)/extent.i;
  temp.k = (n - temp.j * extent.i - temp.i)/(extent.i*extent.j);
  return temp;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Given a grid point i,j,k, finds its location in Cartesian coordinates.
// Parameters:
//   gridPoint gp -- the grid point whose position is sought
// Return Value list:
//   Returns a Vector object containing the location of gp in Cartesian 
//   coordinates.
///////////////////////////////////////////////////////////////////////////////
Vector GenericMap::gridPointToVector(gridPoint gp) {
  Vector temp(0.0, 0.0, 0.0);
  temp.x += gp.i * i.x / extent.i;
  temp.x += gp.j * j.x / extent.j;
  temp.x += gp.k * k.x / extent.k;
  temp.x += offset.x;
  temp.y += gp.i * i.y / extent.i;
  temp.y += gp.j * j.y / extent.j;
  temp.y += gp.k * k.y / extent.k;
  temp.y += offset.y;
  temp.z += gp.i * i.z / extent.i;
  temp.z += gp.j * j.z / extent.j;
  temp.z += gp.k * k.z / extent.k;
  temp.z += offset.z;
  return temp;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculates the correlation function C between the an electron density map
//   and an EDStructure object.  Because the calculation of this function is 
//   slow, evaluate it only at grid points within EDStructure::cutoff of an 
//   atom; those are the only points whose contribution to C will be non-zero.
// Parameters:
//   EDStructure &pdb -- the EDStructure object whose correlation with the ED
//                       map is to be computed
// Return Value List:
//   The correlation is a real number between 0.0 and 1.0, with a score of 1.0 
//   indicating perfect correlation and 0.0 indicating no overlap at all 
//   between the two electron densities.  This function returns ln(1-C), which
//   shows more variation in the region of interest.
///////////////////////////////////////////////////////////////////////////////
double GenericMap::correlate(EDStructure &pdb) {
  using namespace std;
  double cNumerator = 0.0;     // numerator of C
  double cDenominator1 = 0.0;  // first sum in denominator of C
  double cDenominator2 = 0.0;  // second sum in denominator of C
  pdb.coarseGrid->update();
  for (long n = 0; n < dataSize; n++) {
    double localDensity = pdb.densityAt(gridPointToVector(gridIndex(n)));
    cNumerator += gridData[n]*localDensity;
    cDenominator1 += gridData[n]*gridData[n];
    cDenominator2 += localDensity*localDensity;
  }  
  double result = cNumerator/sqrt(cDenominator1*cDenominator2);
  return result;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   After the map has been initialized, call this to output the map in EZD
//   format.
// Parameters:
//   string fileName -- the filename to output
///////////////////////////////////////////////////////////////////////////////
void GenericMap::writeEZD(string fileName) {
  ofstream outFile;
  outFile.open(fileName.c_str());
  if (!outFile.is_open()) {
    cerr << "Failed to output ED map to " << fileName << endl;
    exit(EXIT_FAILURE);
  }
  cout << "Writing ED map to " << fileName << endl;
  outFile << "EZD_MAP\n";
  outFile << "! Generated from theoretical electron density map.\n";
  outFile << "CELL " << a << ' ' << b << ' ' << c << ' ';
  outFile << alpha << ' ' << beta << ' ' << gamma << endl;
  outFile << "ORIGIN " << origin.i << ' ' << origin.j << ' ' << origin.k << endl;
  outFile << "EXTENT " << extent.i << ' ' << extent.j << ' ' << extent.k << endl;
  outFile << "GRID " << extent.i << ' ' << extent.j << ' ' << extent.k << endl;
  // find maximum value of electron density
  int gridSize = extent.i*extent.j*extent.k;
  double maxDensity = 0;
  for (long n = 0; n < gridSize; n++) {
    if (gridData[n] > maxDensity) {
      maxDensity = gridData[n];
    }
  }
  double scale = 9999.0 / maxDensity;
  outFile << "SCALE ";
  outFile << scientific << setprecision(8) << uppercase << scale << endl;
  outFile << "MAP\n";
  for (long n = 0; n < gridSize; n++) {
    outFile << fixed << setprecision(1) << gridData[n] * scale;
    if (((n+1) % 7 == 0)||n==gridSize-1) {  // seven numbers per line
      outFile << endl;
    }
    else {
      outFile << ' ';
    }
  }
  outFile << "END\n";
  outFile.close();
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Trims a map down by elminating some data elements.  Used when an 
//   experimental ED map has more data than is needed for efficient fitting.
// Parameters:
//   int scalingFactor    -- The factor by which the amount of data should be 
//                           scaled down.  For example, a scalingFactor of 3
//                           means that every 3rd data point is kept.
////////////////////////////////////////////////////////////////////////////////
void GenericMap::trimMap(int scalingFactor) {
  // no need to change a,b,c,alpha,beta,gamma,i,j,k
  cout << "Keeping every " << scalingFactor << " data points..." << endl;
  gridPoint newExtent;
  newExtent.i = extent.i / scalingFactor;
  newExtent.j = extent.j / scalingFactor;
  newExtent.k = extent.k / scalingFactor;
  int newDataSize = newExtent.i * newExtent.j * newExtent.k;  
  double * newGridData = new double[newDataSize];
  int n_new = 0;
  for (int n = 0; n < dataSize; n++) {
    gridPoint gp = gridIndex(n);
    if (((gp.i+1) % scalingFactor == 0) && ((gp.j+1) % scalingFactor == 0)
        && ((gp.k+1) % scalingFactor == 0)) {
      newGridData[n_new] = gridData[n];
      n_new++;
    }
  }
  // the offset is the position of (0,0,0) in the new map; this should 
  // correspond to i=j=k=scalingFactor-1 in the old map
  // CJ TODO: The alignment still isn't quite perfect if scalingFactor isn't 
  //          a factor of inputMap.extent -- not sure if this is fixable
  gridPoint newZero;
  newZero.i = newZero.j = newZero.k = scalingFactor - 1;
  Vector newOffset = gridPointToVector(newZero);
  gridPoint newOrigin;
  newOrigin.i = int (newOffset.x * newExtent.i / a);
  newOrigin.j = int (newOffset.y * newExtent.j / b);
  newOrigin.k = int (newOffset.z * newExtent.k / c);
  // now change old to new values
  extent = newExtent;
  dataSize = newDataSize;
  offset = newOffset;
  origin = newOrigin;
  delete [] gridData;
  gridData = newGridData;
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for EZDMap class.  Loads in the information given in the file 
//   header (unit cell parameters, origin, extent, number of grid points, and 
//   scaling factor), loads the ED information from the map, calculates the 
//   unit cell vectors and the real-space offset from the origin.
// Parameters:
//   char *fileName -- name of the EZD-format file containing the ED map
///////////////////////////////////////////////////////////////////////////////
EZDMap::EZDMap(string fileName) {
  using namespace std;
  char err1[36] = "Error in file header: problem with ";
  char err2[22] = " line.\n";
  ifstream inFile;
  inFile.open(fileName.c_str());
  // First, make sure the file is kosher
  if (!inFile.is_open()) {
    cerr << "Failed to open file " << fileName << " for reading.\n";
    exit(EXIT_FAILURE);
  }
  string test;
  inFile >> test;  // first line should be "EZD_MAP"
  if (test!="EZD_MAP") {
    cerr << "Error in file header: not an EZD map!\n";
    inFile.close();
    exit(EXIT_FAILURE);
  }
  // Extract information from the file header
  cout << "Loading header information from " << fileName << " ...\n";
  inFile.get();                // throw away \n at end of first line
  inFile.getline(comment,80);  // second line is a comment
  inFile >> test;              // third line should begin with "CELL"
  if (test!="CELL") {
    cerr << err1 << "CELL" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> a;                // get unit cell parameters
  inFile >> b;
  inFile >> c;
  inFile >> alpha;
  inFile >> beta;
  inFile >> gamma;
  inFile >> test;              // fourth line should begin with "ORIGIN"
  if (test != "ORIGIN")	{
    cerr << err1 << "ORIGIN" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> origin.i;         // get origin
  inFile >> origin.j;
  inFile >> origin.k;
  inFile >> test;             // fifth line should begin with "EXTENT"
  if (test != "EXTENT") {
    cerr << err1 << "EXTENT" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> extent.i;        // get grid extent
  inFile >> extent.j;
  inFile >> extent.k;
  inFile >> test;             // sixth line should begin with "GRID"
  if (test != "GRID")	{
    cerr << err1 << "GRID" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> num.i;        // get number of grid points in each direction
  inFile >> num.j;
  inFile >> num.k;
  inFile >> test;             // seventh line should begin with "SCALE"
  if (test != "SCALE") {
    cerr << err1 << "SCALE" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  inFile >> scale;      // get scaling factor
  // Load in grid data from file
  cout << "Loading grid data from " << fileName << ":   0%";
  dataSize = extent.i*extent.j*extent.k;  // # of points on grid
  gridData = new double[dataSize];
  inFile >> test;
  if (test!="MAP") {  // eighth line should begin with "MAP"
    cerr << err1 << "MAP" << err2;
    inFile.close();
    exit(EXIT_FAILURE);
  }
  int lastPC = 0; // previous percentage complete; see below
  for (long n = 0; n < dataSize; n++) {
    double inputNumber;
    inFile >> inputNumber;
    gridData[n] = inputNumber / scale;
    // Progress report for loading big files
    if (n % 1000 == 0) {
      int percentComplete = int (100.0 * n / dataSize);
      if ((percentComplete % 5 == 0)&&(percentComplete > lastPC)) {
        cout << "\b\b\b\b" << setw(3) << percentComplete << '%' << flush;
	lastPC = percentComplete;
      }
    }
  }
  cout << endl;
  // Calculate unit-cell vectors i, j, and k
  double twoPi360 = 0.0174532925199;
  double radAlpha, radBeta, radGamma;
  radAlpha = twoPi360 * alpha;
  radBeta = twoPi360 * beta;
  radGamma = twoPi360 * gamma;
  double iScale = a * extent.i / num.i;
  double jScale = b * extent.j / num.j;
  double kScale = c * extent.k / num.k;
  i = Vector(iScale, 0.0, 0.0);
  j = Vector(jScale*cos(radGamma), jScale*sin(radGamma), 0.0);
  k.x = kScale*cos(radBeta);
  k.y = kScale*(cos(radAlpha) - cos(radBeta)*cos(radGamma))/sin(radGamma);
  k.z = sqrt(pow(kScale,2) - pow(k.x,2) - pow(k.y,2));
  // Now calculate the real-space offset
  offset.x = (i.x * (origin.i) + j.x * (origin.j) + k.x * origin.k)/extent.i;
  offset.y = (i.y * (origin.i) + j.y * (origin.j) + k.y * origin.k)/extent.j;
  offset.z = (i.z * (origin.i) + j.z * (origin.j) + k.z * origin.k)/extent.k;
  displayHeader();

  inFile.close();
  return;
}
///////////////////////////////////////////////////////////////////////////////
// Description:
//   Displays the information contained in the EZD file header
///////////////////////////////////////////////////////////////////////////////
void EZDMap::displayHeader() {
  using namespace std;

  if( parameters.verbose >= 2 ){
    cout << "  Information from EZD header file:" << endl << endl;
    cout << comment << endl;
    cout << "  a = " << a << ", b = " << b << ", c = " << c << endl;
    cout << "  alpha = " << alpha << ", beta = " << beta << ", gamma = " << gamma << endl;
    cout << "  Origin: " << origin.i << ", " << origin.j << ", " << origin.k << endl;
    cout << "  Grid extent is " << extent.i << " x " << extent.j << " x " << extent.k << endl;
    cout << "  Num of grid points: " << num.i << " x " << num.j << " x " << num.k << endl;
    cout << "  Scaling factor = " << scale << "\n";
  }
  return;
}

///////////////////////////////////////////////////////////////////////////////
// Description:
//   Constructor for the TheoMap class.  Creates a theoretical ED map starting 
//   from a PDB structure.  In the function prototype, resolutionFactorInput has 
//   a default value of 0.0 and gridDensity has a default value of 1.0
///////////////////////////////////////////////////////////////////////////////
TheoMap::TheoMap(string fileName, double resolutionFactorInput, double gridDensity, double noiseLevel) {
  using namespace std;
  // First, load the PDB file into memory
  EDStructure pdb(fileName,resolutionFactorInput);
  // Now, construct the unit cell.  Choice of unit cell is arbitrary; 
  // choose a cubic cell that extends 5 Angstroms beyond the largest 
  // dimensions of the protein in each dimension with gridDensity grid
  // points per Angstrom
  alpha = 90.0;
  beta = 90.0;
  gamma = 90.0;
  a = int (pdb.max.x - pdb.min.x + 10.0);
  b = int (pdb.max.y - pdb.min.y + 10.0);
  c = int (pdb.max.z - pdb.min.z + 10.0);
  i.x = a;   i.y = 0.0; i.z = 0.0;
  j.x = 0.0; j.y = b;   j.z = 0.0;
  k.x = 0.0; k.y = 0.0; k.z = c;
  extent.i = int (ceil(a*gridDensity));
  extent.j = int (ceil(b*gridDensity));
  extent.k = int (ceil(c*gridDensity));
  cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
  cout << "extent = " << extent.i << ',' << extent.j << ',' << extent.k << endl;
  cout << "i = " << i.x << ',' << i.y << ',' << i.z << endl;
  cout << "j = " << j.x << ',' << j.y << ',' << j.z << endl;
  cout << "k = " << k.x << ',' << k.y << ',' << k.z << endl;
  // To get the offset, remember that the i=j=k=0 point
  // will have the lowest values of x,y, and z; in this
  // case they'll be 5 Angstroms below the minimum in the 
  // PDB file.  Rounding to an integer value is needed
  // to get the EZD output to align correctly.
  offset.x = round(pdb.min.x - 5.0);
  offset.y = round(pdb.min.y - 5.0);
  offset.z = round(pdb.min.z - 5.0);
  origin.i = int (offset.x * extent.i / a);
  origin.j = int (offset.y * extent.j / b);
  origin.k = int (offset.z * extent.k / c);
  dataSize = extent.i*extent.j*extent.k;  // # of points on grid
  if (dataSize == 1) {
    cerr << "ERROR: Only one point in ED grid!" << endl;
    cerr << "       Theoretical ED map will contain no data." << endl;
  }
  gridData = new double[dataSize];
  gridPoint gp;
  cout << "Constructing theoretical ED map...     ";
  // The math here is the inverse of generic_map::grid_index()
  long pointsDone = 0;
  double totalDensity = 0.0;
  double totalNoise = 0.0;
  for (int i = 0; i < extent.i; i++) {
    gp.i = i;
    for (int j = 0; j < extent.j; j++) {
      gp.j = j;
      for (int k = 0; k < extent.k; k++) {
	gp.k = k;
        double localDensity = pdb.densityAt(gridPointToVector(gp));
        if (noiseLevel > 0.0) { // add Gaussian noise to data
          double r = genrand_real3();
          double phi = genrand_real3();
          // The width of the noise Gaussian is proportional to the square
          // root of the local density with coefficient noiseLevel.
          // Use a Box-Muller transformation to turn the uniform distribution
          // from the random number generator into a Gaussian distribution.
          double noise = noiseLevel*sqrt(localDensity)*cos(2*PI*phi)*sqrt(-2*log(r));
          totalDensity += localDensity;
          totalNoise += abs(noise);
          localDensity += noise;
        }
	gridData[i + extent.i*j + extent.i*extent.j*k] = localDensity;
        pointsDone++;
        if (dataSize / 100 > 0) {
          if (pointsDone % (dataSize / 100) == 0) {
            cout << "\b\b\b\b" << setw(3) << 100 * pointsDone / dataSize << '%' << flush;
          }
        }
      }
    }
  }
  cout << "\b\b\b\b100%" << endl;
  if (noiseLevel > 0.0) {
    cout << "Signal/Noise ratio: " << totalDensity / totalNoise << endl;
  }
  return;
}


////////////////////////////////////////////////////////////////////////////////
// Setup routine, taken from runFRODA, for theoMap
////////////////////////////////////////////////////////////////////////////////
void Froda::setupTheoEDObjects( MolFramework &structure, EDStructure* &myPdb, GenericMap* &targetMap ){
    double gridDensity = 10;  // Default 0.1 Angstrom grid
    /*if (edResFac >= 2.46954 && edResFac < 3.22552) {
      gridDensity = 1/3.0; // use 3 Angstrom grid for 7-8 A resolution
    }
    else if (edResFac >= 1.25997 && edResFac < 2.46954) {
      gridDensity = 0.5; // use 2 Angstrom grid for 5-7 A resolution
    }
    else if (edResFac < 1.25997) {
      gridDensity = 1; // use 1 Angstrom grid for < 5 Angstrom resolution
      }*/
    targetMap = new TheoMap(theoMapFile,edResFac,gridDensity,edNoise);
    targetMap->writeEZD(theoMapFile.substr(0,theoMapFile.find(".pdb",1))+".ezd");
    myPdb = new EDStructure(*this,structure,*targetMap,edResFac);
    myPdb->cachedCorrelation = targetMap->correlate(*myPdb);
    cout << "Initial correlation with ED map: " << myPdb->cachedCorrelation << endl;
  return;
}


////////////////////////////////////////////////////////////////////////////////
// Setup routine, taken from runFRODA, for EZDMap
////////////////////////////////////////////////////////////////////////////////
void Froda::setupRealEDObjects( MolFramework &structure, EDStructure* &myPdb, GenericMap* &targetMap ){
    targetMap = new EZDMap(ezdMapFile);
    myPdb = new EDStructure(*this,structure,*targetMap,edResFac);
    myPdb->cachedCorrelation = targetMap->correlate(*myPdb);
    cout << "Initial correlation with ED map: " << myPdb->cachedCorrelation << endl;
  return;
}

////////////////////////////////////////////////////////////////////////////////
// Setup routine, taken from runFRODA, for trimMap
////////////////////////////////////////////////////////////////////////////////
void Froda::setupEDTrimMap( MolFramework &structure, EDStructure* &myPdb, GenericMap* &targetMap ){
    if (targetMap == NULL) {
      cerr << "ERROR: Can't use -trimmap without loading an ED map!" << endl;
      exit(EXIT_FAILURE);
    } 
    targetMap->trimMap(trimMapFactor);
    targetMap->writeEZD("trimmed.ezd");
    myPdb->cachedCorrelation = targetMap->correlate(*myPdb);
    cout << "Correlation with trimmed ED map: " << myPdb->cachedCorrelation << endl;
  return;
}



////////////////////////////////////////////////////////////////////////////////
// Destructor
////////////////////////////////////////////////////////////////////////////////
TheoMap::~TheoMap(){

  delete [] gridData;
  cout << "Destroying TheoMap" << endl;
}


   

  
  
  
