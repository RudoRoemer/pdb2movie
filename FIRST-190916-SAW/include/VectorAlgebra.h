////////////////////////////////////////////////////////////////////////////////
// Description:
//   VectorAlgebra.h - class to hadle common vector algebra calculations
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

#ifndef _VECTOR_ALGEBRA_H_
#define _VECTOR_ALGEBRA_H_

#include <math.h>
#include <cstdlib>
//#include "global_defs.h" 

class VectorAlgebra {
public:
	static float dot(float *, float *, float *period=NULL);
	static float length(float *, float *period=NULL);	
	static float lengthSquared(float *, float *period=NULL);
	static float *subtract(float *, float*, float *period=NULL);
	static float distance(float *, float*, float *period=NULL);
        static float distanceSquared(float *, float*, float *period=NULL);
};

#endif // _VECTOR_ALGEBRA_H_
