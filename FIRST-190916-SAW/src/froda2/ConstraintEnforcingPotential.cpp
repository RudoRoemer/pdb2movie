#include "ConstraintEnforcingPotential.h"
#include "VL.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"

ConstraintEnforcingPotential::ConstraintEnforcingPotential(
      VL *verletList__,
      Energy *energy__,
      Gradient *gradient__,
      Mismatch *mismatch__ ) :
  verletList_( verletList__ ),
  energy_( energy__ ),
  gradient_( gradient__ ),
  mismatch_( mismatch__ )
{
}

ConstraintEnforcingPotential::~ConstraintEnforcingPotential()
{
  delete verletList_;
  delete energy_;
  delete gradient_;
  delete mismatch_;
}
