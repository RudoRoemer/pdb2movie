#ifndef GRADIENT_H_
#define GRADIENT_H_

#include <list>
#include <vector>
#include "GeneralizedCoords.h"
#include "Observable.h"
#include "Vec3.h"
class RigidUnitSystem;

struct SecondDerivative {
  double d2V_dx2;
  double d2V_dy2;
  double d2V_dz2;
  double d2V_dxdy;
  double d2V_dydz;
  double d2V_dzdx;
};

class Gradient;
class GradientTerm
{
public:
  GradientTerm() {}
  virtual ~GradientTerm() {}  
  virtual void addToGradient( std::vector<Vec3> &dV_dr_rigidUnitPoint,
      std::vector<SecondDerivative> &secondDerivative_rigidUnitPoint ) = 0;
};

class Gradient : public Observer {
public:
  Gradient( RigidUnitSystem *rigidUnitSystem_ );
  virtual ~Gradient();
  void addTerm( GradientTerm* term ) { gradientTerms.push_back(term); }
  const GeneralizedCoords& calc() {
    if ( !isUpdated ) update();
    return gradientComponents;
  }
  const GeneralizedCoords& operator()() { return calc(); }
  const GeneralizedCoords& calc_d2V_dQ2_diagonal() {
    if ( !isUpdated ) update();
    return d2V_dQ2_diagonal;
  }
  void update();
  void receiveNotification( Observable *observable );

private:
  RigidUnitSystem *rigidUnitSystem;
  std::list<GradientTerm*> gradientTerms;
  
  //first derivative
  std::vector<Vec3> dV_dr_rigidUnitPoint;
  GeneralizedCoords gradientComponents;
  
  //second derivative
  std::vector<SecondDerivative> secondDerivative_rigidUnitPoint; 
  GeneralizedCoords d2V_dQ2_diagonal;
  bool isUpdated;
  double rotormagcutoff2;  
  void applyChainRule();
  
};


#endif /*GRADIENT_H_*/
