#pragma once
#ifndef T3DATAHOLDER_H
#define T3DATAHOLDER_H

#include <chrono>
#include <iostream>
#include <cmath>
#include <cstring>
#include "unistd.h"

#include "T3Globals.h"
#include "T3RNG.h"
#include "T3Particle.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3ThreeVector.h"
#include "T3Process.h"
#include "T3ProcessImplFromCSandFS.h"
#include "T3MultipleScatteringCSImpl.h"
#include "T3MultipleScatteringFSImpl.h"

using namespace t3;

//e - blue p - red photon - green n - goluboi

// build with: cmake . -DCMAKE_C_COMPILER=pgcc -DCMAKE_CXX_COMPILER=pgc++
// -DCMAKE_CXX_FLAGS="-acc -mcmodel=medium -ta=tesla:cc30 -Mnollvm -Minline
// -Mcuda=cuda10.1" -DCMAKE_CXX_STANDARD=17 -DACC=ON -DCUDA=ON -DPROBLEM=OFF/ON

//build from the command line:
//pgc++ -acc -mcmodel=medium -ta=tesla:cc30 -Mnollvm -Minline -Mcuda=cuda10.1 -std=c++17 -g -DOPENACC -DCUDA -DPROBLEM *.cpp -o Test
//pgc++ -std=c++17 -g -DPROBLEM *.cpp -o Test

#define MIN(a, b) ((a<b)?a:b)

//=======================================================================//
using Multiple_scattering_t =
  t3::Process<t3::ProcessImplFromCSandFS<t3::MultipleScatteringCS<FloatingType>,
                                         t3::MultipleScatteringFS<true, FloatingType>>>;

// typename is necessaru for the compiler
// Multiple_scattering_t::Base_t::CS_t()
// and Multiple_scattering_t::Base_t::FS_t()
// are the anonymous temporary objects (r-values).
Multiple_scattering_t multiplescatteringProcess =
      Multiple_scattering_t(typename Multiple_scattering_t::Base_t::CS_t(),
                   typename Multiple_scattering_t::Base_t::FS_t());

//=======================================================================//


t3::PDG_t extrapdg[GL];
FloatingType extram[GL];

//********************************************//
//Host variables and arrays:
//********************************************//
long int MAX_ELEMENT;
unsigned int SHIFT;
unsigned int POSITION3;
unsigned int POSITION2;
unsigned int POSITION1;
unsigned int POSITION0;
unsigned int POSITION23;
int LIFE=0;


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
//********************************************//
//End of host variables and arrays.
//********************************************//

class DataHolder
{
public:
  Particle<FloatingType> particles[GL] __attribute__((aligned(64)));
  int ind01[Nbin][BLt] __attribute__((aligned(64)));
  int ind23[Nbin][BLt] __attribute__((aligned(64)));
  Particle<FloatingType> arr1[GL] __attribute__((aligned(64)));
  Particle<FloatingType> arr2[GL] __attribute__((aligned(64)));
  Particle<FloatingType> arr3[GL] __attribute__((aligned(64)));  
  ParticleTable aParticleTable;
  MaterialTable aMaterialTable;
  FloatingType csMultipleScattering[GL];
  t3::PDG_t outPDG1[GL];
  t3::PDG_t outPDG2[GL];
  t3::LorentzVector<FloatingType> outP1[GL];
  t3::LorentzVector<FloatingType> outP2[GL];
//array for T3MultipleScatteringFSImpl.h for storing
//csBorder arrays for each particle. In T3MultipleScatteringFSImpl.h
//the arrays are allocated dynamically, and here, i think it will
//be allocated statically.
  FloatingType csBorderDataCS[GL];//this is for GetCS() in Propagator().
  FloatingType csBorderDataFS[GL];  //this is for GetFS() in Reactor().


  
  void propagator();
  void reactor();
  void InitParticle();
  void compressor();
  void injector();
private:

};

void DataHolder::propagator()
{
//CALCULATE CROSS SECTIONS:
  auto const material=t3::MatID_t(2u);
  multiplescatteringProcess.GetCS(particles, material, LIFE,
                                  /*for storing cross sections on all isotopes of the material*/
                                  csMultipleScattering, csBorderDataCS,
                                  aParticleTable, aMaterialTable);

#ifdef OPENACC
#pragma acc parallel loop copy(LIFE) present(particles,csMultipleScattering)
#else
#pragma omp parallel for
#endif
  for(int i=0; i<LIFE; ++i)
  {
    ThreeVector<FloatingType> r0(particles[i].r.x(), particles[i].r.y(), particles[i].r.z());
    float   P= particles[i].GetEtot();
    ThreeVector<FloatingType> v(particles[i].p.x()/P, particles[i].p.y()/P, particles[i].p.z()/P);
    //l1 l2 l3
    float l1x;
    float l1y;
    float l1z;
    if(v.x()>=0.0f) l1x=((particles[i].ix()+1)*ag-r0.x())/v.x();
    else            l1x=(particles[i].ix()*ag-r0.x())/v.x();
    if(v.y()>=0.0f) l1y=((particles[i].jy()+1)*ag-r0.y())/v.y();
    else            l1y=(particles[i].jy()*ag-r0.y())/v.y();
    if(v.z()>=0.0f) l1z=((particles[i].kz()+1)*ag-r0.z())/v.z();
    else            l1z=(particles[i].kz()*ag-r0.z())/v.z();
    float l1=MIN(MIN(l1x, l1y), MIN(l1y, l1z));
    double f2=particles[i].GenerateCanonical();
    float l2=fabs(lambda/csMultipleScattering[i]*log(f2));
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
    ThreeVector<FloatingType> e1;
    if(v.x()<v.y())
    {
      if(v.x()<v.z()) e1={0., v.z(), -v.y()};
      else            e1={v.y(), -v.x(), 0.};
    }
    else
    {
      if(v.y()<v.z()) e1={v.z(), 0., -v.x()};
      else            e1={v.y(), -v.x(), 0.};
    }
    ThreeVector<FloatingType> e2=cross(v,e1);
    e1.Unit();
    e2.Unit();
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
    ThreeVector<FloatingType> t1=v*cost;
    ThreeVector<FloatingType> t2=e1*sint*sinf(phi);
    ThreeVector<FloatingType> t3=e2*sint*cosf(phi);
    ThreeVector<FloatingType> new_vel;
    new_vel.SetX(t1.x()+t2.x()+t3.x());
    new_vel.SetY(t1.y()+t2.y()+t3.y());
    new_vel.SetZ(t1.z()+t2.z()+t3.z());
    new_vel.Unit();
    float loss=0.;
    float dl=fabs(l);
    if(particles[i].pdg!=PDG_t(22) && particles[i].pdg!=PDG_t(2112)) loss=rho0*dl;
    if(loss>=particles[i].GetEtot() || particles[i].GetEtot()<=1.)
    {
      particles[i].de=particles[i].GetEtot();
      //ERROR:
      ///\\\///particles[i].p.SetE(0.0);
#ifdef PROBLEM1
      particles[i].SetEtot(0.0, aParticleTable, extrapdg, extram, i);
#endif
      particles[i].ir=0;
    }
    else
    {
      //ERROR:
      ///\\\///particles[i].p.SetE(particles[i].p.E()-loss);
#ifdef PROBLEM1
      particles[i].SetEtot(particles[i].p.E()-loss, aParticleTable, extrapdg, extram, i);
#endif      
      particles[i].de=loss;
      particles[i].ir=irc;
    }
    P=particles[i].p.E();
    particles[i].p.SetX(new_vel.x()*P);
    particles[i].p.SetY(new_vel.y()*P);
    particles[i].p.SetZ(new_vel.z()*P);
    particles[i].r.SetX(particles[i].p.x()/P*dl+da);
    particles[i].r.SetY(particles[i].p.y()/P*dl+da);
    particles[i].r.SetZ(particles[i].p.z()/P*dl+da);
  }
}//End of Propagator

void DataHolder::reactor()
{
  const auto material=t3::MatID_t(2);
  multiplescatteringProcess.GetFS(particles, material, POSITION23,
                                  outPDG1, outP1, outPDG2, outP2, csBorderDataFS,
                                  aParticleTable, aMaterialTable);

  #pragma acc parallel loop
  for (unsigned int i = 0; i < POSITION23; ++i) {
    particles[i].p = outP1[i];
    particles[i+SHIFT].p = outP2[i];
  }  
}//End of Reactor

void DataHolder::InitParticle()
{
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) present(particles[0:GL]) copy(LIFE,MAX_ELEMENT)
  {
#endif
    particles[LIFE].id=MAX_ELEMENT;
    particles[LIFE].rs=RandGen(MAX_ELEMENT+27);
    particles[LIFE].r.SetX(0.0f);
    particles[LIFE].r.SetY(0.0f);
    particles[LIFE].r.SetZ(0.0f);
    particles[LIFE].r.SetE(0.0f);
    particles[LIFE].p.SetE(G);
    particles[LIFE].de=0.0f;
    particles[LIFE].p.SetX(1.0f/sqrtf(3.0)*particles[LIFE].p.E());
    particles[LIFE].p.SetY(1.0f/sqrtf(2.0)*particles[LIFE].p.E());
    particles[LIFE].p.SetZ(1.0f/sqrtf(6.0)*particles[LIFE].p.E());
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

void DataHolder::compressor()
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
    particles[i+SHIFT].r.SetX(particles[i].r.x());
    particles[i+SHIFT].r.SetY(particles[i].r.y());
    particles[i+SHIFT].r.SetZ(particles[i].r.z());    
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
    particles[i+SHIFT].r.SetX(particles[i].r.x());
    particles[i+SHIFT].r.SetY(particles[i].r.y());
    particles[i+SHIFT].r.SetZ(particles[i].r.z());
    particles[i+SHIFT].id=MAX_ELEMENT+i;
    particles[i+SHIFT].rs=RandGen(MAX_ELEMENT+i);
    particles[i+SHIFT].ir=-1;
    if(particles[i].pdg==PDG_t(22)) particles[i+SHIFT].pdg=PDG_t(11);
    else                            particles[i+SHIFT].pdg=PDG_t(22);
  }
  MAX_ELEMENT+=POSITION23-POSITION3;
}//End of Compressor

void DataHolder::injector()
{
  double sg=0.;
#ifdef OPENACC
#pragma acc data copy(sg)//// without this line the program fails with error 9
#pragma acc parallel loop reduction(+:sg) present(particles)
#else
#pragma omp parallel for reduction(+:sg)
#endif
  for(int j=0; j<LIFE; ++j) sg+=particles[j].p.E();
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


#endif//T3DATAHOLDER_H
