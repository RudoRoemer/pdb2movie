#include "RandomRotorPerturber.h"
#include "RandomVector.h"
#include "RigidUnitSystem.h"
#include <algorithm>
#include "Vec3.h"
#include <cmath>
#include "mt19937ar.h"

using namespace std;

RandomRotorPerturber::RandomRotorPerturber( RigidUnitSystem *rigidUnitSystem_,
                                            double size_ ) :
  rigidUnitSystem( rigidUnitSystem_ ),
  size(size_) {}

void RandomRotorPerturber::perturb() {
  size_t nRU = rigidUnitSystem->nRigidUnits();
  
  for ( size_t ru = 0; ru < nRU; ru++ ) {
    if ( rigidUnitSystem->hasZeroRadius(ru) ) continue;

    Vec3 rotorPerturbation;
    double maxRotationAngle = min( 3.141592653589793238462643383279, size/rigidUnitSystem->radius(ru) );
    double randomRotationAngle = genrand_real1()*maxRotationAngle;
    double rotorMagnitude = 2.0 * sin(randomRotationAngle/2.0);
    generateRandomUnitVector( rotorPerturbation );
    rotorPerturbation *= rotorMagnitude;
    
    rigidUnitSystem->addToRotor( ru, rotorPerturbation );
  } 
  rigidUnitSystem->collapseRotors();
  rigidUnitSystem->update();
}
