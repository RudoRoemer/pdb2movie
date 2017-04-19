// SymmetricPerturber.cpp
// Perturbs atoms consistent with symmetry specified by
// BioMT matrices
// Craig Jolley, January 2008
////////////////////////////////////////////////////////////////////////////////

#include "RandomCenterPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include "FitPerturber.h"
#include "SymmetricPerturber.h"
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"
#include "SymmetryMatrices.h"

using namespace std;

SymmetricPerturber::SymmetricPerturber(RigidUnitSystem *rigidUnitSystem_,
                                       SymmetryMatrices *symMat_,double size_ ) {
  rigidUnitSystem = rigidUnitSystem_;
  symMat = symMat_;
  size = size_;
  fitPert = new FitPerturber(rigidUnitSystem);
}

SymmetricPerturber::~SymmetricPerturber() {
  delete fitPert;
}

const Vec3 SymmetricPerturber::transform(Vec3 v, Matrix m) {
  Vec3 newV;
  newV.x = m[0][3];
  newV.y = m[1][3];
  newV.z = m[2][3]; 
  for (int j = 0; j < 3; j++) {
    newV.x += v.x*m[0][j];
    newV.y += v.y*m[1][j];
    newV.z += v.z*m[2][j];
  }
  return newV;
}

void SymmetricPerturber::perturb() {
  vector <Vec3> newPositions = rigidUnitSystem->meanPositions();
  int nMatrices = symMat->size();
  size_t monomerAtoms = newPositions.size() / nMatrices;
  for (size_t atomNum = 0; atomNum < monomerAtoms; atomNum++) {
    // randomly perturb atom in monomer 0
    Vec3 atomPerturbation;
    generateRandomUnitVector(atomPerturbation);
    atomPerturbation *= genrand_real2()*size;
    newPositions[atomNum] += atomPerturbation;
    // now perturb the symmetry-related atoms
    for (int m = 1; m < nMatrices; m++) {
      int n = m * monomerAtoms + atomNum; 
      Vec3 symPerturbation = transform(atomPerturbation,symMat->getMatrix(m));
      newPositions[n] += symPerturbation;
    }
  }
  // atoms have been moved; use FitPerturber to fit rigid units
  fitPert->setMeanPointsTarget(&newPositions);
  fitPert->perturb();
  return;
}
