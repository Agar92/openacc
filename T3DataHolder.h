#pragma once
#ifndef T3DATAHOLDER_H
#define T3DATAHOLDER_H

#include <algorithm>
#include <array>
#include <assert.h>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <cmath>
#include <omp.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include <iomanip>
#include <cfloat>

#include "T3ArrayView.h"
#include "T3Process.h"
#include "T3ProcessImplFromCSandFS.h"
#include "T3WrappedProcess.h"
#include "T3MultipleScatteringCSImpl.h"
#include "T3MultipleScatteringFSImpl.h"

#include "T3Defs.h"
#include "T3ParticleTable.h"

#include "T3ThreeVector.h"
#include "T3LorentzVector.h"

#ifdef OPENACC
#include <accelmath.h>
#include <openacc.h>
#include <cuda.h>
#include <cuda_runtime.h>
#endif

//----------------------------------------------//
//Random number generator:
//----------------------------------------------//
class RNDGenerator
{
private:
  unsigned int fseed;
public:
  void seed(unsigned int seed){fseed=seed;}
  using result_type=unsigned int;
  RNDGenerator():fseed(1){}
  RNDGenerator(unsigned int seed):fseed(seed){}
  result_type min() const {return 0u;}
  result_type max() const {return 0xFFFFFFFF;}
  result_type Rand32(unsigned int xn)
  {
    u_quad_t a=0x5DEECE66D;
    u_quad_t c=0xB;
    return (unsigned int)((a*xn+c) & 0xFFFFFFFF);
  }
  result_type operator()()
  {
    fseed=Rand32(fseed);
    return fseed;
  }
};

using FloatingType = double;
using RandGen = RNDGenerator;

using namespace t3::units;

namespace dataholder {

  constexpr double TLS=10.0*MeV;
  constexpr bool report = true;
  constexpr bool reportTimes = false;
  constexpr bool histogram = true;
  constexpr long int G = 27;
  constexpr int N = 99999;
  constexpr int Np = N+1;
  constexpr int INJ = 1000;
  constexpr long int K = N+1;
  constexpr unsigned int max_loop = 10;
  constexpr int cuba = 16;
  constexpr unsigned int Nbin = 8;
  constexpr int DN = (N+1) / Nbin +1;
  constexpr unsigned int BLt = 200000000 / Nbin;
  constexpr unsigned int GL1 = BLt - 1;
  constexpr double cgam = 5.0;
  constexpr auto dcgam = cgam + cgam;
  constexpr int cubn = cuba + cuba;
  constexpr int cub2 = cubn * cubn;
  constexpr int cub3 = cub2 * cubn;
  constexpr double fuse = .25;
//---------------------------//
  unsigned int SHIFT=0;
  unsigned int POSITION1=0;
  unsigned int POSITION2=0;
  unsigned int POSITION23=0;
  unsigned int POSITION3=0;
  int INJECTED_PARTICLES=0;
  long int MAX_ELEMENT=0;
  unsigned int Ntop=0;
  unsigned int LIFE=0;
  double SumDGam=0.0;
  double GAMMA=0.0;
  double NoNew=0.0;
  unsigned int mini[Nbin];
  unsigned int count01[Nbin];
  unsigned int count23[Nbin];
  unsigned int count0[Nbin];
  unsigned int count1[Nbin];
  unsigned int count2[Nbin];
  unsigned int count3[Nbin];
  unsigned int ii1[Nbin];
  unsigned int ii3[Nbin];
  unsigned int ii23[Nbin];
  unsigned int init[Nbin];
  unsigned int fin[Nbin];
  unsigned int pointer1[Nbin];
  unsigned int pointer2[Nbin];
  unsigned int pointer3[Nbin];
  decltype(INJECTED_PARTICLES) GetNumOfInjectedParticles(){ return INJECTED_PARTICLES; }
  decltype(LIFE) GetNumOfAliveParticles(){ return LIFE; }
  decltype(NoNew) GetNoNew(){return NoNew;}
  decltype(SumDGam) GetSumDGam(){return SumDGam;}
  decltype(Ntop) GetNtop(){return Ntop;}
  int GetNp(){return Np;}
  //bins in histogram:
  const int BinNumber1=100;
  const int BinNumber2=200;
  constexpr FloatingType InitParticlex0 = 0.5;
  constexpr FloatingType InitParticley0 = 0.5;
  
} // end of namespace dataholder

constexpr FloatingType ag = 1.0e-5 * cm;//the width of the cell in cm

template <typename Floating> struct Particle {
  Particle()
      : Particle<Floating>(t3::LorentzVector<Floating>(0.0,0.0,0.0,0.0),
                           t3::LorentzVector<Floating>(1.0/sqrt(3.0),1.0/sqrt(2.0),1.0/sqrt(6.0),dataholder::G),
                           0.0, t3::PDG_t(2112), 1.0, 0, 1, 1) {}

  Particle(t3::LorentzVector<Floating> _r, t3::LorentzVector<Floating> _p, Floating _de, t3::PDG_t _pdg,
           Floating _wt, RandGen::result_type _rs, int _ir, int _id);
  t3::LorentzVector<Floating> r;// r = {rx, ry, rz, t}
  t3::LorentzVector<Floating> p;// p = {px, py, pz, en}
  Floating de;
  Floating wt;
  RandGen rs;
  t3::PDG_t pdg;
  int ir;
  void SetEtot(Floating Etot, t3::ParticleTable aParticleTable)
  {
    const auto m = aParticleTable.GetMass(pdg);
    const auto pnew = sqrt(Etot * Etot - m*m);
    p.SetPxPyPzE(pnew*vx(), pnew*vy(), pnew*vz(), Etot);
  }
  Floating GetdE() const { return de; }
  Floating GetEtot() const { return p.E(); }
  Floating vx() const { return p.x() / p.R(); }
  Floating vy() const { return p.y() / p.R(); }
  Floating vz() const { return p.z() / p.R(); }
  int ix() const { return static_cast<int>(std::floor(r.x() / ag)); }
  int jy() const { return static_cast<int>(std::floor(r.y() / ag)); }
  int kz() const { return static_cast<int>(std::floor(r.z() / ag)); }
  int id;
  auto GenerateCanonical() {
     return t3::GenerateSubCanonical<Floating>(rs);
  }
};

template <typename Floating>
Particle<Floating>::Particle(t3::LorentzVector<Floating> _r, t3::LorentzVector<Floating> _p, Floating _de,
                             t3::PDG_t _pdg, Floating _wt, RandGen::result_type _rs, int _ir, int _id)
    : r{_r}, p{_p}, de(_de), pdg(_pdg), wt(_wt), rs(_rs), ir(_ir), id(_id){}

namespace dataholder {

constexpr unsigned int sizep = sizeof(Particle<FloatingType>);
  
template <typename Floating> class DataHolder {
public:
  Particle<Floating> particles[Np+DN];
private:
  unsigned int ind01[Nbin][BLt];
  unsigned int ind23[Nbin][BLt];
  Particle<Floating> arr1[Np];
  Particle<Floating> arr2[Np];
  Particle<Floating> arr3[Np];
  t3::PDG_t outPDG1[Np];
  t3::LorentzVector<Floating> outP1[Np];
  t3::PDG_t outPDG2[Np];
  t3::LorentzVector<Floating> outP2[Np];
  Floating csIsotropic[Np];
  Floating csDuplicate[Np];
  Floating csMultipleScattering[Np];

  t3::ParticleTable aParticleTable;
  t3::MaterialTable aMaterialTable;
  //for histogram:
//-------------------------------------------------
  int HTheta[cubn][BinNumber1];
  Floating OxTheta[BinNumber1];
  const double HThetaMin=1.0e-4, HThetaMax=5.0e-2;
  double lnHThetaMin, lnHThetaMax;
  double deltaTheta;
//-------------------------------------------------
  int HThetax[cubn][BinNumber2];
  Floating OxThetax[BinNumber2];
  const double HThetaxMin=-0.05, HThetaxMax=0.05;
  double deltaThetax;
//-------------------------------------------------
  int HThetay[cubn][BinNumber2];
  Floating OxThetay[BinNumber2];
  const double HThetayMin=-0.05, HThetayMax=0.05;
  double deltaThetay;
  //end of for histogram
  
  using Multiple_scattering_t =
      t3::ProcessFromCSandFS<t3::MultipleScatteringCS<Floating>,
                         t3::MultipleScatteringFS<true, Floating>>;

  Multiple_scattering_t multiplescatteringProcess =
      Multiple_scattering_t(typename Multiple_scattering_t::Base_t::CS_t(),
                   typename Multiple_scattering_t::Base_t::FS_t());
public:
  DataHolder()
    : aParticleTable(), aMaterialTable(), particles{}, ind01{}, ind23{}, arr1{}, arr2{}, arr3{},
      outPDG1{}, outPDG2{}, outP1{}, outP2{}, csIsotropic{}, csDuplicate{},csMultipleScattering{},
      HTheta{}, HThetax{}, HThetay{} {

      //intialize histograms:
 //-------------------------------------------------------
      lnHThetaMin=std::log(HThetaMin);
      lnHThetaMax=std::log(HThetaMax);
      deltaTheta = (lnHThetaMax - lnHThetaMin)/BinNumber1;
      for(int m=0; m<BinNumber1; ++m)
        OxTheta[m]=std::exp(lnHThetaMin+deltaTheta*(m+0.5));
//-------------------------------------------------------
      deltaThetax = (HThetaxMax - HThetaxMin)/BinNumber2;
      for(int m=0; m<BinNumber2; ++m)
        OxThetax[m]=HThetaxMin+deltaThetax*(m+0.5);
      deltaThetay = (HThetayMax - HThetayMin)/BinNumber2;
      for(int m=0; m<BinNumber2; ++m)
        OxThetay[m]=HThetayMin+deltaThetay*(m+0.5);
//-------------------------------------------------------      
  }
  DataHolder(DataHolder const &) = delete;
  ~DataHolder()
  {
  }
  void Propagate();
  void Compress();
  void React();
  void Inject();
  void InitParticle();
  void Histogram_theta();
};

template <typename Floating> void DataHolder<Floating>::Histogram_theta() {
  std::ofstream foutne_theta1;
  foutne_theta1.open("ne_theta1.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_theta1<<m<<"   ";
    for(int j=0; j<cuba; ++j)
    {
      foutne_theta1<<std::setw(3)<<HTheta[j][m]<<" ";
    }
    foutne_theta1<<std::endl;
  }
  foutne_theta1.close();
  std::ofstream foutne_theta2;
  foutne_theta2.open("ne_theta2.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_theta2<<m<<"   ";
    for(int j=cuba; j<cubn; ++j)
    {
      foutne_theta2<<std::setw(3)<<HTheta[j][m]<<" ";
    }
    foutne_theta2<<std::endl;
  }
  foutne_theta2.close();
  std::ofstream foutne_thetai;
  foutne_thetai.open("ne_thetai.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_thetai<<m<<"   ";
    foutne_thetai<<std::setw(10)<<OxTheta[m]<<" ";
    foutne_thetai<<std::endl;
  }
  foutne_thetai.close();
}

template <typename Floating> void DataHolder<Floating>::InitParticle() {
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(LIFE,MAX_ELEMENT,INJECTED_PARTICLES) present(particles)
{
#endif
  t3::PDG_t initPDG = aParticleTable.makePDGfromZandA(1, 2);
  auto const m = aParticleTable.GetMass(initPDG);
  auto const Tkinls = TLS;
  auto const InitEls= m+Tkinls;
  const auto pls=std::sqrt(InitEls*InitEls-m*m);
  particles[LIFE] = Particle<Floating>(
     t3::LorentzVector<Floating>(InitParticlex0*ag,
       InitParticley0*ag,-cuba*ag,0.0),
       t3::LorentzVector<Floating>(0.0,0.0,pls,InitEls),
       0.0, initPDG, 1., MAX_ELEMENT, 1, MAX_ELEMENT-1);
  ++LIFE;
  ++MAX_ELEMENT;
  ++INJECTED_PARTICLES;
#ifdef OPENACC
}
#endif
} // End of make_particle

template <typename Floating> void DataHolder<Floating>::Compress() {
  auto SumDGam_ = SumDGam;
#ifdef OPENACC
#pragma acc parallel loop reduction(+ : SumDGam_) present(particles)
#else
#pragma omp parallel for simd reduction(+ : SumDGam_)
#pragma distribute_point
#endif
  for (size_t i = 0; i < LIFE; ++i) {
    // FIXME de might be uninitialized
    SumDGam_ += static_cast<decltype(SumDGam_)>(particles[i].de) *
                static_cast<decltype(SumDGam_)>(particles[i].wt);
  }
  SumDGam = SumDGam_;
  auto const dL = LIFE / Nbin;
  auto const DL = dL + 1;
  unsigned int const n = Nbin - LIFE % Nbin;

  POSITION1=POSITION2=POSITION3=POSITION23=0;

#ifdef OPENACC
#pragma acc parallel loop copy(GL1,dL,DL,n,count01,count23,count0,count1,count2,count3,init,fin)
#else
#pragma omp parallel for simd
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    count01[b] = GL1;
    count23[b] = 0;
    count0[b] = GL1;
    count1[b] = 0;
    count2[b] = GL1;
    count3[b] = 0;
    if (b < n) {
      init[b] = b * dL;
      fin[b] = (b + 1) * dL;
    } else if (b == n) {
      init[b] = n * dL;
      fin[b] = n * dL + DL;
    } else if (b > n) {
      init[b] = n * dL + DL * (b - n);
      fin[b] = n * dL + DL * (b - n + 1);
    }
  }

#ifdef OPENACC
#pragma acc parallel loop copy(init,fin,count01,count23) present(particles,ind23)
#else
#pragma omp parallel for simd
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    for (unsigned int i = init[b]; i < fin[b]; ++i) {
      if (particles[i].ir < 2)
        ind23[b][count01[b]--] = i; // FIXME count01.at(b) maybe 0?
      else
        ind23[b][count23[b]++] = i;
    }
  }
#ifdef OPENACC
#pragma acc parallel loop copy(count0,count1,count2,count3,count01,count23,mini,ii23,init,fin) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    ii23[b] = count23[b] - 1;
    mini[b] = GL1 - count01[b];
    if (count23[b] < mini[b])
      mini[b] = count23[b];
    unsigned int js = 0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+ : js)
#endif
    for (unsigned int j = 0; j < mini[b]; ++j)
      if (ind23[b][ii23[b] - j] > ind23[b][GL1 - j])
        ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
    for (unsigned int j = 0; j < js; ++j)
      std::swap(particles[ind23[b][ii23[b] - j]],
                particles[ind23[b][GL1 - j]]);
    for (unsigned int i = init[b]; i < fin[b]; ++i) {
      if (particles[i].ir == 0)
        ind01[b][count0[b]--] = i;
      else if (particles[i].ir == 1)
        ind01[b][count1[b]++] = i;
      else if (particles[i].ir == 2)
        ind23[b][count2[b]--] = i;
      else
        ind23[b][count3[b]++] = i;
    }
  }

#ifdef OPENACC
#pragma acc parallel loop copy(GL1,count0,count1,count2,count3,mini,ii1,ii3) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    ii1[b] = count1[b] - 1;
    mini[b] = GL1 - count0[b];
    if (count1[b] < mini[b])
      mini[b] = count1[b];
    unsigned int js = 0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+ : js)
#endif
    for (unsigned int j = 0; j < mini[b]; ++j)
      if (ind01[b][ii1[b] - j] > ind01[b][GL1 - j])
        ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
    for (unsigned int j = 0; j < js; ++j)
      std::swap(particles[ind01[b][ii1[b] - j]], particles[ind01[b][GL1 - j]]);
    ii3[b] = count3[b] - 1;
    mini[b] = GL1 - count2[b];
    if (count3[b] < mini[b])
      mini[b] = count3[b];
    js = 0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+ : js)
#endif
    for (unsigned int j = 0; j < mini[b]; ++j)
      if (ind23[b][ii3[b] - j] > ind23[b][GL1 - j])
        ++js;

#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
    for (unsigned int j = 0; j < js; ++j)
      std::swap(particles[ind23[b][ii3[b] - j]], particles[ind23[b][GL1 - j]]);
  }
  auto POSITION0 = 0u;
  auto POSITION1_ = 0u;
  auto POSITION2_ = 0u;
  auto POSITION3_ = 0u;
  auto POSITION23_ = 0u;
#ifdef OPENACC
#pragma acc parallel loop reduction(+:POSITION0, POSITION1, POSITION2, POSITION3, POSITION23)
#else
#pragma omp parallel for simd reduction(+:POSITION0,POSITION1_,POSITION2_,POSITION3_,POSITION23_)
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    count0[b] = GL1 - count0[b];
    count2[b] = GL1 - count2[b];
    POSITION0 += count0[b];
    POSITION1_ += count1[b];
    POSITION2_ += count2[b];
    POSITION3_ += count3[b];
    POSITION23_ += count23[b];
  }
  POSITION1 = POSITION1_;
  POSITION2 = POSITION2_;
  POSITION3 = POSITION3_;
  POSITION23 = POSITION23_;

  auto prevLIFE = LIFE;
  SHIFT = LIFE - POSITION0;

  if(histogram) LIFE = LIFE - POSITION0;
  else          LIFE = LIFE + POSITION23 - POSITION0;

#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(pointer1,pointer2,pointer3)
{
#endif
  pointer1[0] = pointer2[0] = pointer3[0] = 0;
  for (unsigned int b = 0; b < Nbin - 1; ++b) {
    pointer1[b + 1] = pointer1[b] + count1[b];
    pointer2[b + 1] = pointer2[b] + count2[b];
    pointer3[b + 1] = pointer3[b] + count3[b];
  }
#ifdef OPENACC
}
#endif

  // DO NOT parallelize or vectorize - undefined behavior
  for (unsigned int b = 0; b < Nbin; ++b) {
#ifdef OPENACC
    #pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&arr1[pointer1[b]], &particles[init[b] + count3[b] + count2[b]], count1[b] * sizep, cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2[pointer2[b]],
                   &particles[init[b] + count3[b]],
                   count2[b] * sizep, cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3[pointer3[b]], &particles[init[b]],
                   count3[b] * sizep, cudaMemcpyDeviceToDevice);
    }
#else
    memcpy(&arr1[pointer1[b]], &particles[init[b] + count3[b] + count2[b]], count1[b] * sizep);
    memcpy(&arr2[pointer2[b]],
                 &particles[init[b] + count3[b]],
                 count2[b] * sizep);
    memcpy(&arr3[pointer3[b]], &particles[init[b]],
                 count3[b] * sizep);
#endif
  }

#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
 {
   cudaMemcpy(&particles[0], &arr3[0], POSITION3 * sizep, cudaMemcpyDeviceToDevice);
   cudaMemcpy(&particles[POSITION3], &arr2[0], POSITION2 * sizep, cudaMemcpyDeviceToDevice);
   cudaMemcpy(&particles[POSITION23], &arr1[0], POSITION1 * sizep, cudaMemcpyDeviceToDevice);
 }
#else
 memcpy(&particles[0], &arr3[0], POSITION3 * sizep);
 memcpy(&particles[POSITION3], &arr2[0], POSITION2 * sizep);
 memcpy(&particles[POSITION23], &arr1[0], POSITION1 * sizep);
#endif

 if(!histogram)
 {
#ifdef OPENACC
#pragma acc parallel loop copy(SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (unsigned int i = 0; i < POSITION3; ++i) {
    unsigned int i2 = i + i;
    particles[i+SHIFT] = particles[i];
    particles[i].rs.seed(static_cast<RandGen::result_type>(MAX_ELEMENT + i2));
    particles[i+SHIFT].rs.seed(
        static_cast<RandGen::result_type>(MAX_ELEMENT + i2 + 1));
    if (particles[i].pdg == aParticleTable.makePDGfromZandA(1,2))
    {
      particles[i].pdg = aParticleTable.makePDGfromZandA(1,2);
      particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(1,2);
    }
    else if (particles[i].pdg == aParticleTable.makePDGfromZandA(22,48))
    {
      particles[i].pdg = aParticleTable.makePDGfromZandA(22,48);
      particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(22,48);
    }
    else if (particles[i].pdg == t3::PDG_t(22))
    {
      particles[i].pdg = t3::PDG_t(11);
      particles[i+SHIFT].pdg = t3::PDG_t(-11);
    } else if (particles[i].pdg == t3::PDG_t(11))
    {
      particles[i].pdg = t3::PDG_t(22);
      particles[i+SHIFT].pdg = t3::PDG_t(11);
    } else if (particles[i].pdg == t3::PDG_t(-11))
    {
      particles[i].pdg = t3::PDG_t(22);
      particles[i+SHIFT].pdg = t3::PDG_t(22);
    } else if (particles[i].pdg == t3::PDG_t(2112))
    {
      particles[i].pdg = t3::PDG_t(2112);
      particles[i+SHIFT].pdg = t3::PDG_t(2112);
    }
  }
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(MAX_ELEMENT,POSITION3) present(particles)
{
#endif
  MAX_ELEMENT += POSITION3 + POSITION3;
#ifdef OPENACC
}
#endif

#ifdef OPENACC
#pragma acc parallel loop copy(POSITION3,POSITION23,SHIFT,MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (unsigned int i = POSITION3; i < POSITION23; ++i) {
    particles[i+SHIFT] = particles[i];
    particles[i+SHIFT].rs.seed(static_cast<RandGen::result_type>(MAX_ELEMENT + i));
    if (particles[i].pdg == aParticleTable.makePDGfromZandA(1,2))
      particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(1,2);
    else if (particles[i].pdg == aParticleTable.makePDGfromZandA(22,48))
      particles[i+SHIFT].pdg = aParticleTable.makePDGfromZandA(22,48);
    else if (particles[i].pdg == t3::PDG_t(22))
      particles[i+SHIFT].pdg = t3::PDG_t(11);
    else
      particles[i+SHIFT].pdg = t3::PDG_t(22);
  }
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(MAX_ELEMENT,POSITION23,POSITION3) present(particles) 
{
#endif
  MAX_ELEMENT += POSITION23 - POSITION3;
#ifdef OPENACC
}
#endif
 }//end of (!histogram)
}

template <typename Floating> void DataHolder<Floating>::Propagate() {  
  auto inputlskinE = t3::makeAOSMemberView(
      particles,
      [this](Particle<Floating> const &particle)  {//[this] is necessary to access aParticleTable.
        return particle.GetEtot() - aParticleTable.GetMass(particle.pdg);
      });
  auto inputPDG = t3::makeAOSMemberView(
      particles,
      [](Particle<Floating> const &particle)  {
        return particle.pdg;
      });

  multiplescatteringProcess.GetCS(inputlskinE, inputPDG, t3::MatID_t(2u), LIFE,
                                   csMultipleScattering);

  constexpr Floating da = ag * 1e-10;
  constexpr Floating rho0 = 1.0;
  constexpr auto lMAX = std::numeric_limits<Floating>::max();
  Floating const ntid2 = aMaterialTable.GetConcentrations(t3::MatID_t(2u));
  Floating const rho_tid2 = aMaterialTable.GetDensity(t3::MatID_t(2u));
  Floating DE=0.0;

#ifdef OPENACC
#pragma acc parallel loop copy(da,rho0,lMAX,ntid2,rho_tid2) present(particles)
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (size_t i = 0; i < LIFE; ++i)
  {
    auto En = particles[i].GetEtot();

    if (particles[i].ir > 0)
    {
      auto const csDuplicateI = csDuplicate[i];
      auto const csMultipleScatteringi = csMultipleScattering[i];

      Floating  l1x =
          (particles[i].vx() == 0.)
              ? lMAX
              : ((particles[i].vx() > 0.)
                     ? ((particles[i].ix() + 1) * ag - particles[i].r.x())
                     : (particles[i].ix() * ag - particles[i].r.x())) /
                        particles[i].vx() +
                    da;
      Floating  l1y =
          (particles[i].vy() == 0.)
              ? lMAX
              : ((particles[i].vy() > 0.)
                     ? ((particles[i].jy() + 1) * ag - particles[i].r.y())
                     : (particles[i].jy() * ag - particles[i].r.y())) /
                        particles[i].vy() +
                    da;
      Floating  l1z =
          (particles[i].vz() == 0.)
              ? lMAX
              : ((particles[i].vz() > 0.)
                     ? ((particles[i].kz() + 1) * ag - particles[i].r.z())
                     : (particles[i].kz() * ag - particles[i].r.z())) /
                        particles[i].vz() +
                    da;

      Floating const dEdx = 2.0 * MeV * cm * cm / gr;
      Floating const dEdxfull = dEdx*rho_tid2;
      Floating const range = (particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg))/dEdxfull;
      Floating l0 = range;
      l0=1.0*cm;
      Floating const csMScatteringi = csMultipleScatteringi;
      Floating const lambdar = 1.0/ntid2/csMScatteringi;
      auto const R = particles[i].GenerateCanonical();
      Floating l2 = std::abs(lambdar * std::log(R));
      Floating const l1 = std::min({l1x, l1y, l1z});
      int const irc = (l0 < l2 && l0 < l1) ? 0 : ((l2 < l1) ? 2 : 1);
      Floating const l = std::min({l1, l2, l0});
      Floating const dl = std::abs(l);
      auto const indexz = particles[i].kz();

      if(irc == 1 && std::abs(l-l1z)<da)
      {
          auto const modparticlei = particles[i].p.R();
          auto const costheta = particles[i].p.z()/modparticlei;
          auto const theta = std::acos(costheta);
          const int id = particles[i].id;
          auto const thetax=std::atan( particles[i].p.x()/particles[i].p.z() );//=px/pz.
          auto const thetay=std::atan( particles[i].p.y()/particles[i].p.z() );//=py/pz.
          auto const x=particles[i].r.x() - InitParticlex0 * ag;
          auto const y=particles[i].r.y() - InitParticley0 * ag;
          auto const r=std::sqrt(x*x+y*y);
          const double logtheta=std::log(theta);
          const int ThetaBin=(logtheta-lnHThetaMin)/deltaTheta;
//*          
          if(ThetaBin>=0 && ThetaBin<BinNumber1)
          {
#ifdef OPENACC
            #pragma acc atomic update
#else
            #pragma omp atomic update
#endif
              ++HTheta[indexz+cuba][ThetaBin];
          }
//-------------------------------------------
          const int ThetaxBin=(thetax-HThetaxMin)/deltaThetax;
          if(ThetaxBin>=0 && ThetaxBin<BinNumber2)
#ifdef OPENACC
            #pragma acc atomic update
#else
            #pragma omp atomic update
#endif
              ++HThetax[indexz+cuba][ThetaxBin];          
//-------------------------------------------
          const int ThetayBin=(thetay-HThetayMin)/deltaThetay;
          if(ThetayBin>=0 && ThetayBin<BinNumber2)
#ifdef OPENACC
            #pragma acc atomic update
#else
            #pragma omp atomic update
#endif
              ++HThetay[indexz+cuba][ThetayBin];                   
//*/              
      }      
      particles[i].r.SetPxPyPzE (particles[i].r.x()+particles[i].vx() * (dl + da),
                              particles[i].r.y()+particles[i].vy() * (dl + da),
                              particles[i].r.z()+particles[i].vz() * (dl + da), 0.0);//t=0.0

      bool const out = (particles[i].ix() >= cuba || particles[i].jy() >= cuba ||
                        particles[i].kz() >= cuba || particles[i].ix() < -cuba ||
                        particles[i].jy() < -cuba || particles[i].kz() < -cuba);
      Floating loss = 0.;   
      if (aParticleTable.IsNucleus(particles[i].pdg))
      {
        loss = dEdxfull * dl;
        if (loss >= particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg))
        {
          particles[i].de += particles[i].GetEtot() - aParticleTable.GetMass(particles[i].pdg);
          particles[i].SetEtot(aParticleTable.GetMass(particles[i].pdg), aParticleTable);
          particles[i].ir = 0;
        }
        else
        {
          Floating const old_energy = particles[i].GetEtot();
          Floating const new_energy = old_energy - loss;
          particles[i].SetEtot(new_energy, aParticleTable);
          bool check=false;
          bool cond=indexz*ag<particles[i].r.z() && particles[i].r.z()<=(indexz+1)*ag;
          particles[i].de += loss;
          particles[i].ir = irc;
        }
      }  
      if (out)
      {
        particles[i].ir = 0;
      }
    }//End if ir>0
  }//End of i-particle loop
}

template <typename Floating> void DataHolder<Floating>::React() {  
  auto returnZeroMat = t3::makeConstantView(t3::MatID_t(2));
  auto inputP = t3::makeAOSMemberView(
      particles, [](Particle<Floating> const &particle) {
        return particle.p;
      });
  auto inputPDG = t3::makeAOSMemberView(
      particles, [](Particle<Floating> const &particle) {
        return particle.pdg;
      });
  auto inputRNG = t3::makeAOSMemberView(
      particles,
      [](Particle<Floating> &particle) -> decltype(particle.rs) & {
        return particle.rs;
      });
  auto const outPDGoutput = std::make_tuple(outPDG1, outPDG2);
  auto const outPoutput = std::make_tuple(outP1, outP2);

  multiplescatteringProcess.GetFS(inputP, inputPDG, returnZeroMat, inputRNG,
                         POSITION23, outPDGoutput, outPoutput);

#ifdef OPENACC
#pragma acc parallel loop
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (unsigned int i = 0; i < POSITION23; ++i) {
    particles[i].p = outP1[i];
    particles[i].ir = 1;
    particles[i].pdg = outPDG1[i];
  }
}

template <typename Floating> void DataHolder<Floating>::Inject() {
  double sg = 0.;
  double wt = 0.;
  double sgCompensation = 0.;
  double wtCompensation = 0.;
#ifdef OPENACC
#pragma acc parallel loop reduction(+ : sg, wt)
#else
#pragma omp parallel for simd reduction(+ : sg, wt)
#endif
  for (unsigned int j = 0; j < LIFE; ++j)
  {
    auto sgi = particles[j].GetEtot() * static_cast<double>(particles[j].wt)-sgCompensation;
    auto sgtemp = sg + sgi;
    sgCompensation = (sgtemp - sg) - sgi;
    sg = sgtemp;
    auto wti = static_cast<double>(particles[j].wt) - wtCompensation;
    auto wttemp = wt + wti;
    wtCompensation = (wttemp - wt) - wti;
    wt = wttemp;
  }
  GAMMA = sg;
  NoNew = wt;
  static int push = 0;
  if (LIFE > Ntop)
    Ntop = LIFE;
  if (push < N) {
    if (LIFE < K) {
      for(int i=0; i<INJ; ++i)
      {
        if(push+i<N) InitParticle();
      }
      push+=INJ;
    }
  }
 }
} // end of namespace dataholder
#include "T3MaterialTable.h"
#endif // T3DATAHOLDER_H
