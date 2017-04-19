#ifndef FIT_H_
#define FIT_H_

#include <vector>
#include "Vec3.h"
class MultiVarFunction_Adapter_Fit;
#include "ConjugateGradientMinimizer.h"
#include "Rotator.h"

class Fit
{
public:
	Fit();
	virtual ~Fit();
  
  void setSourceBasePoints( const std::vector<Vec3> &base_, double radiusGyration_ );
  void setSourceAbsolutePoints( const std::vector<Vec3> &absolute_ );
  void setTargetAbsolutePoints( const std::vector<Vec3> &absolute_ );
  void setRotor( const Vec3& rotor_ );
  double energy();
  const Vec3& gradient();
  const Vec3& secondDerivativeDiagonal();
  const Vec3& getRotor() const { return rotor; }
  double getRadiusGyration() const { return radiusGyration; }
  void simpleFit();
  //void doFit();
  const Vec3& getFitRotor() const { return savedRotation.rotor(); }
  const Vec3& getFitCenter() const { return targetAbsoluteCenter; }
  const vector<Vec3>& getRotatedPoints() const { return rotated; }
  
private:
  std::vector<Vec3> base;
  std::vector<Vec3> target;
  std::vector<Vec3> rotated;
  Vec3 rotor;
  Vec3 targetAbsoluteCenter;
  double radiusGyration;
  Vec3 d2V_dB2;
  Vec3 gradrotor;
  double E;
  bool isUpdatedE;
  bool isUpdatedGrad;
  bool isUpdated2ndDeriv;

  MultiVarFunction_Adapter_Fit *func;
  Rotator currentRotation;
  Rotator savedRotation;
  int numberOfSavedRotations;
  ConjugateGradientMinimizer minimizer;
  double rotormagcutoff2;
  bool rotorIsZero;
  bool hasZeroRadius;
  void collapseRotor();
};

#endif /*FIT_H_*/
