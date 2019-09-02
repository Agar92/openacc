#include <chrono>
#include <iostream>
#include "unistd.h"

#include "T3DataHolder.h"

//////////////////////////////////////////////////////////////////////////////
// Program main
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
	auto begin=std::chrono::steady_clock::now();
#ifdef OPENACC
  #ifdef CUDA
    std::cout<<"OPENACC IS DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#else
  #ifdef CUDA
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS DEFINED"<<std::endl;
  #else
    std::cout<<"OPENACC IS NOT DEFINED, CUDA IS NOT DEFINED"<<std::endl;
  #endif
#endif


  DataHolder d;
    
  
  LIFE=0;
  MAX_ELEMENT=0;
  int step=0;
#ifdef OPENACC
//////\\\\\\//////#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,extrapdg,extram,outPDG1,outPDG2,outP1,outP2,csBorderDataCS,csBorderDataFS,csMultipleScattering,aParticleTable,aMaterialTable/*multiplescatteringProcess*/,d) copy(ag,lambda,dtheta0,da,rho0,Dan)
/////\\\\\/////#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,extrapdg,extram,outPDG1,outPDG2,outP1,outP2,csBorderDataCS,csBorderDataFS,csMultipleScattering,d) copy(ag,lambda,dtheta0,da,rho0,Dan)
////\\\\////#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,extrapdg,extram,outPDG1,outPDG2,outP1,outP2,csBorderDataCS,csBorderDataFS) copy(d) copy(ag,lambda,dtheta0,da,rho0,Dan)
///\\\///#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,extrapdg,extram,outPDG1,outPDG2,outP1,outP2) copy(d) copy(ag,lambda,dtheta0,da,rho0,Dan)
//\\//#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3,extrapdg,extram) copy(d) copy(ag,lambda,dtheta0,da,rho0,Dan)
#pragma acc data create(extrapdg,extram) copy(d) copy(ag,lambda,dtheta0,da,rho0,Dan)
  {
#endif
    d.InitParticle();
    while(LIFE>0)
    {
      
      d.injector();
      d.propagator();
      d.compressor();
      d.reactor();
      ++step;
    }
#ifdef OPENACC
  }
#endif
  const double pw=std::pow(N+1,2.0/3.0);
  const double E=G/2;
  const float multiplier=0.1*pw/(1.0+0.093*pw);
  const float res=(0.1*pw/(1.0+0.093*pw)*std::pow(E,0.5357));
  const unsigned int NBIN_OPT=(int)res;
  std::cout<<"NBIN_OPT="<<NBIN_OPT<<" res="<<res
           <<" m="<<multiplier<<" pw="<<pw<<std::endl;
  
	auto end=std::chrono::steady_clock::now();
	auto elapsed_ms=std::chrono::duration_cast<std::chrono::milliseconds>(end-begin);
	std::cout<<"time="<<elapsed_ms.count()<<" ms, G="<<G<<", K="<<K<<", Ntop="
           <<Ntop<<", SumDG="<<SumDGam<<std::endl;
}
