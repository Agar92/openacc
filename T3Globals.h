#pragma once
#ifndef T3GLOBALS_H
#define T3GLOBALS_H

#include "T3RNG.h"

#ifdef OPENACC
#include <accelmath.h>
#include <openacc.h>
#ifdef CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#else
#include <omp.h>
  struct float3{
    float x, y, z;
  };
  struct float3{
    float x, y, z, w;
  };
#endif

  using FloatingType=float;
  using RandGen=t3::RNDGenerator;

//variables:
  const long int G=2000;
  const unsigned int GL=200000;
  const int N=0;
  const long int K=GL;
  const unsigned int Nbin=1;
  const int BLt=GL/Nbin;

  const float ag=1.0f;
  const float lambda=1.0f;
  const float da=ag*1e-4;
  const float rho0=1.0f;
  const float dtheta0=0.01f;
  const float Dan=0.1f;


namespace t3{
  
  typedef int PDG_t;
  typedef int MatID_t;

}//end of namespace t3.

#endif//T3GLOBALS_H
