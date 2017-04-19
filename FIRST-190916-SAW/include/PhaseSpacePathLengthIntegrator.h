#ifndef _PHASESPACEPATHLENGTHINTEGRATOR_H_
#define _PHASESPACEPATHLENGTHINTEGRATOR_H_
#include <vector>
#include <cmath>

#include "Vec3.h"

////////////////////////////////////////////////////////////////////////////////
//
//   Scott Menor
//   Arizona State University
//   Department of Physics and Astronomy
//   Biophysics Theory Group
//   PO Box 871504
//   Tempe, AZ 85287
//   Scott.Menor@asu.edu
//
//   Copyright (c) 2005-2007, Arizona State University, All Rights Reserved.
////////////////////////////////////////////////////////////////////////////////
class PhaseSpacePathLengthIntegrator {
public:
  PhaseSpacePathLengthIntegrator( const std::vector<Vec3> *coordinates );
  ~PhaseSpacePathLengthIntegrator();
  
  void activateMassWeighting(const std::vector<double>& masses_ );
  void activateMasking( const std::vector<char>& boolMask_ );

  void updatePath();
  
  double getIntegratedPathLength() const { return integratedPathLength; }
  double getSegmentLength() const { return length; }
  
private:
  const std::vector<Vec3> *coordinates; 
  std::vector<Vec3> previousCoordinates;
  std::vector<double> masses;
  std::vector<char> boolMask;
  double length;
  double integratedPathLength;
  bool weighting;
  bool masking;

  double computeDistance( 
      const std::vector<Vec3> &initialCoordinates,
      const std::vector<Vec3> &finalCoordinates );
};

#endif // _PHASESPACEPATHLENGTHINTEGRATOR_H_
