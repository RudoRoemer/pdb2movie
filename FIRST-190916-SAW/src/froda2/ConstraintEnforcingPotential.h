#ifndef CONSTRAINTENFORCINGPOTENTIAL_H_
#define CONSTRAINTENFORCINGPOTENTIAL_H_

class VL;
class Energy;
class Gradient;
class Mismatch;
#include "VL.h"
#include "Energy.h"
#include "Gradient.h"
#include "Mismatch.h"
#include "GeneralizedCoords.h"

class ConstraintEnforcingPotential
{
public:
  ConstraintEnforcingPotential(
    VL *verletList__,
    Energy *energy__,
    Gradient *gradient__,
    Mismatch *mismatch__  
  );
  
  virtual ~ConstraintEnforcingPotential();
  
  double energy() { return energy_->calc(); }
  const GeneralizedCoords& gradient() { return gradient_->calc(); }
  const GeneralizedCoords& d2V_dQ2_diagonal() { return gradient_->calc_d2V_dQ2_diagonal(); }
  double maxMismatch() { return mismatch_->calc(); }
  
  Mismatch *mismatch() { return mismatch_; }
  
private:  
  VL *verletList_;
  Energy *energy_;
  Gradient *gradient_;
  Mismatch *mismatch_;   
};

#endif /*CONSTRAINTENFORCINGPOTENTIAL_H_*/
