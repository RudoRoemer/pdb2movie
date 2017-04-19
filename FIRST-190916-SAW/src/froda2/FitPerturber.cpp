#include "FitPerturber.h"
#include "RigidUnitData.h"
#include "RigidUnitLookup.h"
#include "RigidUnitSystem.h"
#include "RigidUnitFitter.h"

using namespace std;

FitPerturber::FitPerturber(
                RigidUnitSystem *rigidUnitSystem_ ) : 
  rigidUnitSystem( rigidUnitSystem_ ),
  fitter( rigidUnitSystem )
{
}

FitPerturber::~FitPerturber()
{
}

void FitPerturber::setMeanPointsTarget( const vector<Vec3>* meanPointsTarget_ ) {
  meanPointsTarget = meanPointsTarget_;
}

void FitPerturber::perturb() {
  rigidUnitSystem->collapseRotors();
  fitter.calcFitToMeanPoints( *meanPointsTarget );
  size_t nRU = rigidUnitSystem->nRigidUnits();
  Vec3 rotor;
  Vec3 deltaCenter;
  for ( size_t ru = 0; ru < nRU; ru++ ) {

    rotor = fitter.getFitRotation( ru );
    deltaCenter = fitter.getFitTranslation( ru );
    double radius = rigidUnitSystem->radius(ru);

    //Sometimes, rigid units having just three points
    //are in a near-straight-line configuration.  Very easily,
    //the perturbation fit can require a large rotation (near 180
    //degrees).  After these perturbations, the subsequent
    //relaxation often needs to find a way to undo the rotation,
    //but this is difficult, and it dramatically affects the 
    //relaxation.  So, here, we check the magnitude of the fit
    //rotation.  If a rotation is large, then we do not
    //perform the rotation as part of the perturbation.
    if ( radius*sqrt( rotor.norm2() ) < 0.5 ) {
      rigidUnitSystem->setRotor( ru, rotor );
      rigidUnitSystem->collapseRotor( ru );
    }
    rigidUnitSystem->addToCenter( ru, deltaCenter );

  }
  rigidUnitSystem->update();

}
