// SymmetricPerturber.h
// Perturbs atoms consistent with symmetry specified by
// BioMT matrices
// Craig Jolley, January 2008
////////////////////////////////////////////////////////////////////////////////

#ifndef SYMMETRICPERTURBER_H_
#define SYMMETRICPERTURBER_H_

#include "SymReplica.h"
#include "RigidUnitSystem.h"
#include "FitPerturber.h"
#include "Vec3.h"
#include "SymmetryMatrices.h"

class RigidUnitSystem;

class SymmetricPerturber {
public:
  SymmetricPerturber(RigidUnitSystem *rigidUnitSystem_, 
                     SymmetryMatrices *symMat_, double size_);
  virtual ~SymmetricPerturber();
  void perturb();
private:
  RigidUnitSystem *rigidUnitSystem;
  SymmetryMatrices *symMat;
  double size;
  FitPerturber *fitPert;
  const Vec3 transform(Vec3 v, Matrix m);
};

#endif /* SYMMETRICPERTURBER_H_ */
