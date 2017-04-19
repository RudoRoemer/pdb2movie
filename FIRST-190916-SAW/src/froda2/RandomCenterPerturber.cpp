#include "RandomCenterPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"

using namespace std;

RandomCenterPerturber::RandomCenterPerturber( RigidUnitSystem *rigidUnitSystem_,
                                              double size_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  size(size_) {}

void RandomCenterPerturber::perturb() {
  size_t nRU = rigidUnitSystem->nRigidUnits();
  
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    Vec3 centerPerturbation;
    generateRandomUnitVector( centerPerturbation );
    centerPerturbation *= genrand_real2()*size;
    rigidUnitSystem->addToCenter( ru, centerPerturbation );
  }
  
  rigidUnitSystem->update();
}
