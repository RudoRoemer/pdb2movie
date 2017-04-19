////////////////////////////////////////////////////////////////////////////////
// Description:
//   VectorAlgebra.cpp - class to hadle common vector algebra calculations
//
//   Scott Menor
//   Arizona State University
//   Department of Physics and Astronomy
//   Biophysics Theory Group
//   PO Box 871504
//   Tempe, AZ 85287
//   Scott.Menor@asu.edu
//
//   Copyright (c) 2005, Arizona State University, All Rights Reserved.

#include "VectorAlgebra.h"

// TODO - ideally this should use STL templates so we can have vector operations on arbitrary stl::vector<>s
// TODO - for now, considering float arrays of length 3 only (TODO - generalize to n-dimensions)

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float VectorAlgebra::dot(float *x, float *y, float *period) {

  float dot = 0;
	
  for (int dimension=0;dimension<3;dimension++) {
    dot += x[dimension] * y[dimension];
  }
  
  return dot;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float VectorAlgebra::lengthSquared( float *x, float *period) {
  
  float lengthSquared = VectorAlgebra::dot(x, x, period);
  return lengthSquared;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float VectorAlgebra::length(float *x, float * period) {
  
  float lengthSquared = VectorAlgebra::lengthSquared(x, period);
  float length = sqrt(lengthSquared);
  
  return length;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float *VectorAlgebra::subtract(float *x, float *y, float * period) {

  float *difference = new float[3];
  
  for (int dimension=0;dimension<3;dimension++) {
    float rawDifference = x[dimension] - y[dimension];
    
    float actualDifference = rawDifference;
    
    if (period != NULL) {
      float dimensionPeriod = period[dimension];
      
      if (dimensionPeriod > 0.0) {
        float xValue = x[dimension];
        float yValue = y[dimension];

        float xInRange = xValue - dimensionPeriod * floor(xValue / dimensionPeriod);
        float yInRange = yValue - dimensionPeriod * floor(yValue / dimensionPeriod);
 
        float xyDifference = xInRange - yInRange;//x[dimension] - y[dimension];
        float yxDifference = dimensionPeriod + yInRange - xInRange;//y[dimension] - x[dimension];

        if (fabs(xyDifference) < fabs(yxDifference)) {
          actualDifference = xyDifference;

        } else {
          actualDifference = yxDifference;

        }
      }
    }
    
    difference[dimension] = actualDifference;
  } 
  
  return difference;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float VectorAlgebra::distanceSquared(float *x, float *y, float * period) {

  float *difference = VectorAlgebra::subtract(x, y, period);
  float distanceSquared = VectorAlgebra::lengthSquared( difference, period);
  delete [] difference;
  return distanceSquared;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
float VectorAlgebra::distance(float *x, float *y, float * period) {

  float * difference = VectorAlgebra::subtract(x, y, period);
  float distance = VectorAlgebra::length(difference, period);
  delete [] difference;
  return distance;
}

