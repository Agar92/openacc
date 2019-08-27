
////#define ACCEL
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <assert.h>
#include <math.h>
#include "unistd.h"

#include <string>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <vector>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <cstring>
#include <cstdint>
#include <climits>
#include <array>
#include <stdio.h>
#include <cmath>
#include <sys/time.h>

#include "T3RNG.h"
#include "T3Globals.h"
#include "T3Particle.h"

using namespace t3;

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

//e - blue p - red photon - green n - goluboi

// build with:

#define MIN(a, b) ((a<b)?a:b)

Particle<FloatingType> particles[GL] __attribute__((aligned(64)));
int ind01[Nbin][BLt] __attribute__((aligned(64)));
int ind23[Nbin][BLt] __attribute__((aligned(64)));
Particle<FloatingType> arr1[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr2[GL] __attribute__((aligned(64)));
Particle<FloatingType> arr3[GL] __attribute__((aligned(64)));
long int MAX_ELEMENT;
unsigned int SHIFT;
unsigned int POSITION3;
unsigned int POSITION2;
unsigned int POSITION1;
unsigned int POSITION0;
unsigned int POSITION23;
int LIFE=0;

ParticleTable aParticleTable;

unsigned int sizep=sizeof(Particle<FloatingType>);
int push=0;
int over=0;
unsigned int Ntop=0;
double SumDGam=0.;
double GAMMA=0.;

int mini[Nbin];
int count01[Nbin];
int count23[Nbin];
int count0[Nbin];
int count1[Nbin];
int count2[Nbin];
int pos2[Nbin];
int count3[Nbin];
int ii1[Nbin];
int ii3[Nbin];
int ii23[Nbin];

int init[Nbin];
int fin[Nbin];

int pointer1[Nbin];
int pointer2[Nbin];
int pointer3[Nbin];
int GL1=BLt-1;
int dL;
int DL;
int n;
int numbin;

float3 scalevec(float3 vector, float scalar)
{
  float3 rt = vector;
  rt.x *= scalar;
  rt.y *= scalar;
  rt.z *= scalar;
  return rt;
}
float dot(float3 v0, float3 v1)
{
  return v0.x*v1.x+v0.y*v1.y+v0.z*v1.z;
}

float3 cross(float3 v0, float3 v1)
{
  float3 rt;
  rt.x = v0.y*v1.z-v0.z*v1.y;
  rt.y = v0.z*v1.x-v0.x*v1.z;
  rt.z = v0.x*v1.y-v0.y*v1.x;
  return rt;
}

float3 normalize(float3 vector)
{
  float dist = sqrtf(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
  float3 rt;
  rt.x=vector.x/dist;
  rt.y=vector.y/dist;
  rt.z=vector.z/dist;
  return rt;
}

unsigned int Rand32(unsigned int xn)
{
  u_quad_t a=0x5DEECE66D;
  u_quad_t c=0xB;
  return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
}

double rndv(unsigned int xn)
{
  return (double) xn / (double) 0x100000000LL;
}

void propagator()
{
#ifdef OPENACC
#pragma acc parallel loop copy(LIFE) present(particles)
#else
#pragma omp parallel for
#endif
  for(int i=0; i<LIFE; ++i)
  {
    float3 r0={particles[i].r.x, particles[i].r.y, particles[i].r.z};
    float   P= particles[i].GetEtot();
    float3 v = {particles[i].p.x/P, particles[i].p.y/P, particles[i].p.z/P};
    //l1 l2 l3
    float l1x;
    float l1y;
    float l1z;
    if(v.x>=0.0f) l1x=((particles[i].ix()+1)*ag-r0.x)/v.x;
    else          l1x=(particles[i].ix()*ag-r0.x)/v.x;
    if(v.y>=0.0f) l1y=((particles[i].jy()+1)*ag-r0.y)/v.y;
    else          l1y=(particles[i].jy()*ag-r0.y)/v.y;
    if(v.z>=0.0f) l1z=((particles[i].kz()+1)*ag-r0.z)/v.z;
    else          l1z=(particles[i].kz()*ag-r0.z)/v.z;
    float l1=MIN(MIN(l1x, l1y), MIN(l1y, l1z));
    double f2=particles[i].GenerateCanonical();
    float l2=fabs(lambda*log(f2));
    double f3=particles[i].GenerateCanonical();
    float l3=fabs(lambda*log(f3));
    int irc;
    float l23;
    if(l3<l2)
    {
      irc=3;
      l23=l3;
    }
    else
    {
      irc=2;
      l23=l2;
    }
    float l=l1;
    if(l23<l1) l=l23;
    else       irc=1;
    //v e1 e2
    float3 e1;
    if(v.x<v.y)
    {
      if(v.x<v.z) e1={0., v.z, -v.y};
      else        e1={v.y, -v.x, 0.};
    }
    else
    {
      if(v.y<v.z) e1={v.z, 0., -v.x};
      else        e1={v.y, -v.x, 0.};
    }
    float3 e2=cross(v,e1);
    e1=normalize(e1);
    e2=normalize(e2);
    double phi=particles[i].GenerateCanonical()*2*M_PI;
    //угол многократного рассеяния
    float p1=float(1.-2*particles[i].GenerateCanonical());
    float p12=p1*p1;
    float lips=sqrtf(1.f-p1*p1)-1.E-5;
    float p2=lips*float(1.-2*particles[i].GenerateCanonical());
    float p=p12+p2*p2;
    if(p>=1.) p=.9999999;
    float theta=fabs(dtheta0*p1*sqrtf(-2*fabs(l)*logf(p)/p));
    if(theta>1.) theta=1.;
    float cost=cosf(theta);
    float sint=0.;
    if(cost<1.f&&cost>-1.f) sint=sqrtf(1.f-cost*cost);
    float3 t1=scalevec(v,cost);
    float3 t2=scalevec(e1,sint*sinf(phi));
    float3 t3=scalevec(e2,sint*cosf(phi));
    float3 new_vel;
    new_vel.x=t1.x+t2.x+t3.x;
    new_vel.y=t1.y+t2.y+t3.y;
    new_vel.z=t1.z+t2.z+t3.z;
    new_vel=normalize(new_vel);
    float loss=0.;
    float dl=fabs(l);
    if(particles[i].pdg!=PDG_t(22) && particles[i].pdg!=PDG_t(2112)) loss=rho0*dl;
    if(loss>=particles[i].GetEtot() || particles[i].GetEtot()<=1.)
    {
      particles[i].de=particles[i].GetEtot();
      particles[i].p.w=0.0;
      particles[i].SetEtot(0.0, aParticleTable);
      particles[i].ir=0;
    }
    else
    {
      particles[i].p.w=particles[i].p.w-loss;
      particles[i].de=loss;
      particles[i].ir=irc;
    }
    P=particles[i].p.w;
    particles[i].p.x=new_vel.x*P;
    particles[i].p.y=new_vel.y*P;
    particles[i].p.z=new_vel.z*P;
    particles[i].r.x+=particles[i].p.x/P*dl+da;
    particles[i].r.y+=particles[i].p.y/P*dl+da;
    particles[i].r.z+=particles[i].p.z/P*dl+da;
  }
}//End of Propagator

void reactor()
{
#ifdef OPENACC
#pragma acc parallel loop copy(SHIFT) present(particles)
#else
#pragma omp parallel for
#endif
  for(int i=0; i<POSITION23; ++i)
  {
    float part=particles[i].GenerateCanonical();
	  float gamma1;
    float gamma2;    
    float gamma=particles[i].p.w;
    float g=gamma-2.0f;
    if(g<0.0f)
    {
      gamma1=g/2 + 1;
      gamma2=gamma1;
    }
    else
    {
	    gamma1=g*part+1.0f;
      gamma2=g*(1.0f-part)+1.0f;
    }
	  float3 e1;
    float3 e2;
    float xy=particles[i].p.x*particles[i].p.y;
    float xz=particles[i].p.x*particles[i].p.z;
    float yz=particles[i].p.y*particles[i].p.z;
    float x2=particles[i].p.x*particles[i].p.x;
    float y2=particles[i].p.y*particles[i].p.y;
    float z2=particles[i].p.z*particles[i].p.z;
	  if(particles[i].p.x < particles[i].p.y)
    {
      if(particles[i].p.x < particles[i].p.z)
      {
        e1={0., particles[i].p.z, -particles[i].p.y};
        e2={y2+z2, -xy, -xz};
      }
      else
      {
        e1={particles[i].p.y, -particles[i].p.x, 0.};
        e2={-xz, -yz, y2+x2};
      }
    }
    else
    {
	    if(particles[i].p.y < particles[i].p.z)
      {
        e1={particles[i].p.z, 0., -particles[i].p.x};
        e2={xy, -x2-z2, yz};
      }
	    else
      {
        e1={particles[i].p.y, -particles[i].p.x, 0.};
        e2={-xz, -yz, y2+x2};
      }
    }
	  e1=normalize(e1);
	  e2=normalize(e2);
	  float phi=particles[i].GenerateCanonical()*2*M_PI;
    // first particle
	  float Theta1=1.;
    if(gamma>Dan) Theta1=Dan/gamma;
    float cost1=cosf(Theta1);
    float sint1=0.;
    if(cost1<1.f&&cost1>-1.f) sint1=sqrtf(1.f-cost1*cost1);
	  float3 v1;
    float sf=sinf(phi);
    float cf=cosf(phi);
    float ss=sint1*sf;
    float sc=sint1*cf;
    float P=particles[i].p.w;
	  v1.x=particles[i].p.x/P*cost1+e1.x*ss+e2.x*sc;
	  v1.y=particles[i].p.y/P*cost1+e1.y*ss+e2.y*sc;
	  v1.z=particles[i].p.z/P*cost1+e1.z*ss+e2.z*sc;
	  v1=normalize(v1);
	  particles[i].p.w=gamma1;
	  float Theta2=1.;
    if(gamma>Dan) Theta2=Dan/gamma;
    float cost2=cosf(Theta2);
    float sint2=0.;
    if(cost2<1.f&&cost2>-1.f) sint2=sqrtf(1.f-cost2*cost2);
	  float3 v2;
    ss=sint2*sf;
    sc=sint2*cf;
	  v2.x=particles[i].p.x/P*cost2-e1.x*ss-e2.x*sc;
	  v2.y=particles[i].p.y/P*cost2-e1.y*ss-e2.y*sc;
	  v2.z=particles[i].p.z/P*cost2-e1.z*ss-e2.z*sc;
	  v2=normalize(v2);
    particles[i].p.x=v1.x*particles[i].p.w;
	  particles[i].p.y=v1.y*particles[i].p.w;
	  particles[i].p.z=v1.z*particles[i].p.w;
    particles[i+SHIFT].p.w=gamma2;
    particles[i+SHIFT].p.x=v2.x*particles[i+SHIFT].p.w;
    particles[i+SHIFT].p.y=v2.y*particles[i+SHIFT].p.w;
    particles[i+SHIFT].p.z=v2.z*particles[i+SHIFT].p.w;
  }
}//End of Reactor

void InitParticle()
{
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) present(particles[0:GL]) copy(LIFE,MAX_ELEMENT)
  {
#endif
    particles[LIFE].id=MAX_ELEMENT;
    particles[LIFE].rs=RandGen(MAX_ELEMENT+27);
    particles[LIFE].r.x=0.0f;
    particles[LIFE].r.y=0.0f;
    particles[LIFE].r.z=0.0f;
    particles[LIFE].r.w=0.0f;
    particles[LIFE].p.w=G;
    particles[LIFE].de=0.0f;
    particles[LIFE].p.x=1.0f/sqrtf(3.0)*particles[LIFE].p.w;
    particles[LIFE].p.y=1.0f/sqrtf(2.0)*particles[LIFE].p.w;
    particles[LIFE].p.z=1.0f/sqrtf(6.0)*particles[LIFE].p.w;
    particles[LIFE].ir=-1;
    particles[LIFE].pdg=PDG_t(22);
#ifdef OPENACC
#pragma acc atomic update
#endif
    ++LIFE;
#ifdef OPENACC
#pragma acc atomic update
#endif
    ++MAX_ELEMENT;
#ifdef OPENACC
  }
#endif
}//End of InitParticle

void compressor()
{
  double sdg=SumDGam;
#ifdef OPENACC
#pragma acc parallel loop vector reduction(+:sdg) present(particles)
#else
#pragma omp parallel for reduction(+:sdg)
#endif
  for(int i=0; i<LIFE; ++i) sdg += particles[i].de;
  SumDGam=sdg;
  dL=LIFE/Nbin;
  DL=dL+1;
  if(LIFE%Nbin==0) n=Nbin;
  else     n=Nbin*DL-LIFE;
    
  POSITION0=POSITION1=POSITION2=POSITION3=POSITION23=0;
#ifdef OPENACC
  #pragma acc parallel loop copy(GL1,dL,DL,n,init,fin)
  {
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      count01[b]=GL1;
      count23[b]=0;
      count0[b]=GL1;
      count1[b]=0;
      count2[b]=GL1;
      count3[b]=0;
      if(b<n)
      {
        init[b]=b*dL;
        fin[b]=(b+1)*dL;
      }
      else if(b==n)
      {
        init[b]=n*dL;
        fin[b]=n*dL+DL;
      }
      else if(b>n)
      {
        init[b]=n*dL+DL*(b-n);
        fin[b]=n*dL+DL*(b-n+1);
      }
    }
#ifdef OPENACC
  }
#endif

#pragma acc parallel loop copy(count01,count23,init,fin) present(particles,ind23)//num_gangs(1) vector_length(1) copy(count01,count23,cnt01,cnt23)
  {
    for(int b=0; b<Nbin; ++b)    
    {
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if(particles[i].ir<2) ind23[b][count01[b]--]=i;
        else                  ind23[b][count23[b]++]=i;
      }
    }
  }

#ifdef OPENACC
#pragma acc parallel loop copy(count0,count1,count2,count3,count01,count23,mini,ii23,init,fin) present(particles,ind01,ind23)//num_gangs(1) vector_length(1)
#else
#pragma omp parallel for
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      ii23[b]=count23[b]-1;
      mini[b]=GL1-count01[b];
      if(count23[b]<mini[b]) mini[b]=count23[b];
      int js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma simd
#endif
      for(int j=0; j<mini[b]; ++j) js+=int(ind23[b][ii23[b]-j]>ind23[b][GL1-j]);
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii23[b]-j]],particles[ind23[b][GL1-j]]);
      
      for(int i=init[b]; i<fin[b]; ++i)
      {
        if     (particles[i].ir==0) ind01[b][count0[b]--]=i;
        else if(particles[i].ir==1) ind01[b][count1[b]++]=i;
        else if(particles[i].ir==2) ind23[b][count2[b]--]=i;
        else                        ind23[b][count3[b]++]=i;
      }
    }
#ifdef OPENACC
#pragma acc parallel loop copy(count0,count1,count2,count3,mini,ii1,ii3,GL1) present(particles)//num_gangs(1) vector_length(1)
#endif
    for(int b=0; b<Nbin; ++b)    
    {
      ii1[b]=count1[b]-1;
      mini[b]=GL1-count0[b];
      if(count1[b]<mini[b]) mini[b]=count1[b];
      int js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma simd
#endif
      for(int j=0; j<mini[b]; ++j) js+=int(ind01[b][ii1[b]-j]>ind01[b][GL1-j]);
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind01[b][ii1[b]-j]],particles[ind01[b][GL1-j]]);
      
      ii3[b]=count3[b]-1;
      mini[b]=GL1-count2[b];
      if(count3[b]<mini[b]) mini[b]=count3[b];
      js=0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma simd
#endif
      for(int j=0; j<mini[b]; ++j) js+=int(ind23[b][ii3[b]-j]>ind23[b][GL1-j]);
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma simd
#endif
      for(int j=0; j<js; ++j) std::swap(particles[ind23[b][ii3[b]-j]],particles[ind23[b][GL1-j]]);
    }
  
  // Reorder the pointers limits
#ifdef OPENACC
#pragma acc parallel loop reduction(+:POSITION0,POSITION1,POSITION2,POSITION3,POSITION23)
#else
#pragma omp parallel for
#endif
  for(int b=0; b<Nbin; ++b)
  {
    count0[b]=GL1-count0[b];
    count2[b]=GL1-count2[b];
    POSITION0+=count0[b];
    POSITION1+=count1[b];
    POSITION2+=count2[b];
    POSITION3+=count3[b];
    POSITION23+=count23[b];
  }
  LIFE-=POSITION0;
  SHIFT=LIFE;
  LIFE+=POSITION23;
  
//здесь должно идти слияние мини ящиков в ящики для 0, 1, 2, 3, удаление 0, и перекладывание 3, 2, 1 в исходный ящик
  pointer1[0]=pointer2[0]=pointer3[0]=0;
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(pointer1[0:Nbin],pointer2[0:Nbin],pointer3[0:Nbin])
#endif
  for(int b=0; b<Nbin-1; ++b)
  {
    pointer1[b+1]=pointer1[b]+count1[b];
    pointer2[b+1]=pointer2[b]+count2[b];
    pointer3[b+1]=pointer3[b]+count3[b];
  }
  
  for(int b=0; b<Nbin; ++b)
  {
#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep,cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep,cudaMemcpyDeviceToDevice);
    }
#else
    memcpy(&arr1[pointer1[b]],&particles[init[b]+count3[b]+count2[b]],count1[b]*sizep);
    memcpy(&arr2[pointer2[b]],&particles[init[b]+count3[b]],count2[b]*sizep);
    memcpy(&arr3[pointer3[b]],&particles[init[b]],count3[b]*sizep);
#endif
  }
  
  // слияние ящиков для 1, 2, 3 в массив particles
#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
  {
    cudaMemcpy(&particles[0],&arr3[0],POSITION3*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep,cudaMemcpyDeviceToDevice);
    cudaMemcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep,cudaMemcpyDeviceToDevice);
  }
#else
  memcpy(&particles[0],&arr3[0],POSITION3*sizep);
  memcpy(&particles[POSITION3],&arr2[0],POSITION2*sizep);
  memcpy(&particles[POSITION23],&arr1[0],POSITION1*sizep);
#endif
  
#ifdef OPENACC
#pragma acc parallel loop copy(SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for
#endif
  for(int i=0; i<POSITION3; ++i)
  {
    int i2=i+i;
    particles[i].id=MAX_ELEMENT+i2;
    particles[i].rs=RandGen(particles[i].id);
    particles[i].ir=-1;
    particles[i+SHIFT].r.x=particles[i].r.x;
    particles[i+SHIFT].r.y=particles[i].r.y;
    particles[i+SHIFT].r.z=particles[i].r.z;
    particles[i+SHIFT].id=MAX_ELEMENT+i2+1;
    particles[i+SHIFT].rs=RandGen(particles[i+SHIFT].id);
    particles[i+SHIFT].ir=-1;
    if(particles[i].pdg==PDG_t(22))
    {
      particles[i].pdg=PDG_t(11);
      particles[i+SHIFT].pdg=PDG_t(-11);
    }
    else if(particles[i].pdg==PDG_t(11))
    {
      particles[i].pdg=PDG_t(22);
      particles[i+SHIFT].pdg=PDG_t(11);
    }
    else if(particles[i].pdg==PDG_t(-11))
    {
      particles[i].pdg=PDG_t(22);
      particles[i+SHIFT].pdg=PDG_t(22);
    }
    else if(particles[i].pdg==PDG_t(2112))
    {
      particles[i].pdg=PDG_t(22);
      particles[i+SHIFT].pdg=PDG_t(22);
    }
  }    
  MAX_ELEMENT+=POSITION3+POSITION3;
#ifdef OPENACC
#pragma acc parallel loop copy(POSITION3,POSITION23,SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for
#endif
  for(int i=POSITION3; i<POSITION23; ++i)
  {
    particles[i].ir=-1;
    particles[i+SHIFT].r.x=particles[i].r.x;
    particles[i+SHIFT].r.y=particles[i].r.y;
    particles[i+SHIFT].r.z=particles[i].r.z;
    particles[i+SHIFT].id=MAX_ELEMENT+i;
    particles[i+SHIFT].rs=RandGen(MAX_ELEMENT+i);
    particles[i+SHIFT].ir=-1;
    if(particles[i].pdg==PDG_t(22)) particles[i+SHIFT].pdg=PDG_t(11);
    else                            particles[i+SHIFT].pdg=PDG_t(22);
  }
  MAX_ELEMENT+=POSITION23-POSITION3;
}//End of Compressor

void injector()
{
  double sg=0.;
#ifdef OPENACC
#pragma acc data copy(sg)//// without this line the program fails with error 9
#pragma acc parallel loop reduction(+:sg) present(particles)
#else
#pragma omp parallel for reduction(+:sg)
#endif
  for(int j=0; j<LIFE; ++j) sg+=particles[j].p.w;
  GAMMA=sg;
  if(LIFE>Ntop) Ntop=LIFE;
  if(push<N)
  {
    if(LIFE<K) {InitParticle(); ++push;}
    std::cout<<" push particle n="<<push<<", L="<<LIFE
             <<",G="<<GAMMA<<",DG="<<SumDGam<<",Tot="<<GAMMA+SumDGam<<std::endl;
  }
  else
  {
    std::cout<<"fin o="<<over++<<", G="<<GAMMA<<", SDG="<<SumDGam<<", L="
                <<LIFE<<", Tot="<<GAMMA+SumDGam<<std::endl;
  }
}//End of Injector


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
  
  LIFE=0;
  MAX_ELEMENT=0;
  int step=0;
#ifdef OPENACC
#pragma acc data create(particles,ind01,ind23,arr1,arr2,arr3) copy(ag,lambda,dtheta0,da,rho0,Dan,aParticleTable)
  {
#endif
    InitParticle();
    while(LIFE>0)
    {
      injector();
      propagator(); 
      compressor();
      reactor();
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
