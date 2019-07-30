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

  constexpr double TLS=10.0*MeV;//0.96*MeV;//energy of the initial particle in InitParticle().
  constexpr bool report = true;
  constexpr bool reportTimes = false;
  constexpr bool histogram = true;
  constexpr long int G = 27;
  constexpr int N = 99999;//14999999;//999;//999999;//999999;//999999;//>INJ
  constexpr int Np = N+1;//total number of injected particles.
  constexpr int INJ = 1000;
  constexpr long int K = N+1;//=G;
  constexpr unsigned int max_loop = 10;//30000;//2700;
  constexpr int cuba = 16; // halfWidth of the Critical Cube in cell lengths(ag)
  constexpr unsigned int Nbin = 8;
  constexpr int DN = (N+1) / Nbin +1;
  constexpr unsigned int BLt = 200000000 / Nbin;
  constexpr unsigned int GL1 = BLt - 1;
  constexpr double cgam = 5.0; // energy of secondary fission neutrons
  constexpr auto dcgam =
      cgam + cgam; // doubled energy of secondary fission neutrons
  constexpr int cubn = cuba + cuba; // dimension for 3D buffers
  constexpr int cub2 =
      cubn * cubn; // squared dimension for 3D reposes (safety factor)
  constexpr int cub3 =
      cub2 * cubn;             // squared dimension for 3D reposes (safety factor)
  constexpr double fuse = .25; // Neutron fusion parameter
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
  std::array<unsigned int, Nbin> mini{0};
  std::array<unsigned int, Nbin> count01{0};
  std::array<unsigned int, Nbin> count23{0};
  std::array<unsigned int, Nbin> count0{0};
  std::array<unsigned int, Nbin> count1{0};
  std::array<unsigned int, Nbin> count2{0};
  std::array<unsigned int, Nbin> count3{0};
  std::array<unsigned int, Nbin> ii1{0};
  std::array<unsigned int, Nbin> ii3{0};
  std::array<unsigned int, Nbin> ii23{0};
  std::array<unsigned int, Nbin> init{0};
  std::array<unsigned int, Nbin> fin{0};
  std::array<unsigned int, Nbin> pointer1{0};
  std::array<unsigned int, Nbin> pointer2{0};
  std::array<unsigned int, Nbin> pointer3{0};
  decltype(INJECTED_PARTICLES) GetNumOfInjectedParticles(){ return INJECTED_PARTICLES; }
  decltype(LIFE) GetNumOfAliveParticles(){ return LIFE; }
  decltype(NoNew) GetNoNew(){return NoNew;}
  decltype(SumDGam) GetSumDGam(){return SumDGam;}
  decltype(Ntop) GetNtop(){return Ntop;}
} // end of namespace dataholder

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
  //direction unit vector (vx(), vy(), vz()).
  Floating vx() const { return p.x() / p.R(); }
  Floating vy() const { return p.y() / p.R(); }
  Floating vz() const { return p.z() / p.R(); }
  int ix() const { return static_cast<int>(std::floor(r.x() / ag)); }
  int jy() const { return static_cast<int>(std::floor(r.y() / ag)); }
  int kz() const { return static_cast<int>(std::floor(r.z() / ag)); }
  int id;
  auto GenerateCanonical() {
    return std::generate_canonical<Floating,
                                   std::numeric_limits<Floating>::digits>(rs);
  }
  static constexpr Floating ag = 1.0e-5 * cm;//the width of the cell in cm
  static constexpr Floating InitParticlex0 = 0.5;
  static constexpr Floating InitParticley0 = 0.5;
};

template <typename Floating>
//replaced//Particle<Floating>::Particle(Floating _rx, Floating _ry, Floating _rz,
//                             Floating _ct, Floating _vx, Floating _vy,
//                             Floating _vz, Floating _en, Floating _de, t3::PDG_t _pdg,
//                             Floating _wt, RandGen::result_type _rs, int _ir)
//    : r{_rx, _ry, _rz, _ct}, p{_vx, _vy, _vz, _en},
//      de(_de), pdg(_pdg), wt(_wt), rs(_rs), ir(_ir){};
Particle<Floating>::Particle(t3::LorentzVector<Floating> _r, t3::LorentzVector<Floating> _p, Floating _de,
                             t3::PDG_t _pdg, Floating _wt, RandGen::result_type _rs, int _ir, int _id)
    : r{_r}, p{_p}, de(_de), pdg(_pdg), wt(_wt), rs(_rs), ir(_ir), id(_id){}

namespace dataholder {
template <typename Floating> class DataHolder;

template <typename Floating> class ParticleVector {
public:
  auto &at(size_t ind) { return fParticles.at(ind); }
  auto const &at(size_t ind) const { return fParticles.at(ind); }
  auto& back() { return fParticles.back(); }
  auto const& back() const { return fParticles.back(); }
  auto data() { return fParticles.data(); }
  auto size() const { return fParticles.size(); }
  //added
  auto begin() { return fParticles.begin(); }
  auto end() { return fParticles.end(); }
private:
  std::array<Particle<Floating>, Np+DN> fParticles;
};

template <typename Floating> class DataHolder {
  constexpr static unsigned int sizep = sizeof(Particle<Floating>);
public:
  ParticleVector<Floating> particles;
private:
  std::array<std::array<unsigned int, BLt>, Nbin> ind01;
  std::array<std::array<unsigned int, BLt>, Nbin> ind23;
  std::array<Particle<Floating>, Np> arr1;
  std::array<Particle<Floating>, Np> arr2;
  std::array<Particle<Floating>, Np> arr3;
  std::array<t3::PDG_t, Np> outPDG1;
  std::array<t3::LorentzVector<Floating>, Np> outP1;
  std::array<t3::PDG_t, Np> outPDG2;
  std::array<t3::LorentzVector<Floating>, Np> outP2;

  t3::ParticleTable aParticleTable;
  t3::MaterialTable aMaterialTable;
public:
  ParticleVector<Floating> & GetParticles(){return particles;}
  static int GetNp(){return Np;}
private:
  //for histpogram:
//-------------------------------------------------
  static const auto BinNumber1=100;
  std::array<std::array<int, BinNumber1>, cubn> HTheta;
  std::array<Floating, BinNumber1> OxTheta;
  const double HThetaMin=1.0e-4, HThetaMax=/*5 **/ 5.0e-2;
  double lnHThetaMin, lnHThetaMax;
  double deltaTheta;
//-------------------------------------------------
  static const auto BinNumber2=200;
  std::array<std::array<int, BinNumber2>, cubn> HThetax;
  std::array<Floating, BinNumber2> OxThetax;
  const double HThetaxMin=-/*5 **/ 0.05, HThetaxMax=/*5 **/ 0.05;
  double deltaThetax;
//-------------------------------------------------
  std::array<std::array<int, BinNumber2>, cubn> HThetay;
  std::array<Floating, BinNumber2> OxThetay;
  const double HThetayMin=-/*5 **/ 0.05, HThetayMax=/*5 **/ 0.05;//2x-?
  double deltaThetay;
  //end of for histogram

  using Multiple_scattering_t =
      t3::ProcessFromCSandFS<t3::MultipleScatteringCS<Floating>,
                         t3::MultipleScatteringFS<true, Floating>>;

  Multiple_scattering_t multiplescatteringProcess =
      Multiple_scattering_t(typename Multiple_scattering_t::Base_t::CS_t(),
                   typename Multiple_scattering_t::Base_t::FS_t());

  static constexpr Floating ag = Particle<Floating>::ag;

public:
  DataHolder()
    : aParticleTable(), aMaterialTable(), particles{}, ind01{}, ind23{}, arr1{}, arr2{}, arr3{},
      outPDG1{}, outPDG2{}, outP1{}, outP2{}, HTheta{}, HThetax{}, HThetay{} {

      //intialize histograms:
 //-------------------------------------------------------
      lnHThetaMin=std::log(HThetaMin);
      lnHThetaMax=std::log(HThetaMax);
      deltaTheta = (lnHThetaMax - lnHThetaMin)/BinNumber1;
      for(int m=0; m<BinNumber1; ++m)
        OxTheta.at(m)=std::exp(lnHThetaMin+deltaTheta*(m+0.5));
//-------------------------------------------------------
      deltaThetax = (HThetaxMax - HThetaxMin)/BinNumber2;
      for(int m=0; m<BinNumber2; ++m)
        OxThetax.at(m)=HThetaxMin+deltaThetax*(m+0.5);
      deltaThetay = (HThetayMax - HThetayMin)/BinNumber2;
      for(int m=0; m<BinNumber2; ++m)
        OxThetay.at(m)=HThetayMin+deltaThetay*(m+0.5);
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
  //1st part
  std::ofstream foutne_theta1;
  foutne_theta1.open("ne_theta1.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_theta1<<m<<"   ";
    for(int j=0; j<cuba; ++j)
    {
      foutne_theta1<<std::setw(3)<<HTheta.at(j).at(m)<<" ";
    }
    foutne_theta1<<std::endl;
  }
  foutne_theta1.close();
  //2nd part
  std::ofstream foutne_theta2;
  foutne_theta2.open("ne_theta2.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_theta2<<m<<"   ";
    for(int j=cuba; j<cubn; ++j)
    {
      foutne_theta2<<std::setw(3)<<HTheta.at(j).at(m)<<" ";
    }
    foutne_theta2<<std::endl;
  }
  foutne_theta2.close();
  //Collect thetai:
  //1st part
  std::ofstream foutne_thetai;
  foutne_thetai.open("ne_thetai.dat");
  for(int m=0; m<BinNumber1; ++m)
  {
    foutne_thetai<<m<<"   ";
    foutne_thetai<<std::setw(10)<<OxTheta.at(m)<<" ";
    foutne_thetai<<std::endl;
  }
  foutne_thetai.close();
}

template <typename Floating> void DataHolder<Floating>::InitParticle() {
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) present(particles, LIFE, MAX_ELEMENT)
{
#endif
  t3::PDG_t initPDG = aParticleTable.makePDGfromZandA(1, 2);
  auto const m = aParticleTable.GetMass(initPDG);
  auto const Tkinls = TLS;
  auto const InitEls= m+Tkinls;
  const auto pls=std::sqrt(InitEls*InitEls-m*m);
  particles.at(LIFE) = Particle<Floating>(
     t3::LorentzVector<Floating>(Particle<Floating>::InitParticlex0*ag,
       Particle<Floating>::InitParticley0*ag,-cuba*ag,0.0),
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
    SumDGam_ += static_cast<decltype(SumDGam_)>(particles.at(i).de) *
                static_cast<decltype(SumDGam_)>(particles.at(i).wt);
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
    count01.at(b) = GL1;
    count23.at(b) = 0;
    count0.at(b) = GL1;
    count1.at(b) = 0;
    count2.at(b) = GL1;
    count3.at(b) = 0;
    if (b < n) {
      init.at(b) = b * dL;
      fin.at(b) = (b + 1) * dL;
    } else if (b == n) {
      init.at(b) = n * dL;
      fin.at(b) = n * dL + DL;
    } else if (b > n) {
      init.at(b) = n * dL + DL * (b - n);
      fin.at(b) = n * dL + DL * (b - n + 1);
    }
  }

#ifdef OPENACC
#pragma acc parallel loop copy(init,fin,count01,count23) present(particles,ind23)
#else
#pragma omp parallel for simd
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    for (unsigned int i = init.at(b); i < fin.at(b); ++i) {
      if (particles.at(i).ir < 2)
        ind23.at(b).at(count01.at(b)--) = i; // FIXME count01.at(b) maybe 0?
      else
        ind23.at(b).at(count23.at(b)++) = i;
    }
  }
#ifdef OPENACC
#pragma acc parallel loop copy(count0,count1,count2,count3,count01,count23,mini,ii23,init,fin) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    ii23.at(b) = count23.at(b) - 1;
    mini.at(b) = GL1 - count01.at(b);
    if (count23.at(b) < mini.at(b))
      mini.at(b) = count23.at(b);
    unsigned int js = 0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+ : js)
#endif
    for (unsigned int j = 0; j < mini.at(b); ++j)
      if (ind23.at(b).at(ii23.at(b) - j) > ind23.at(b).at(GL1 - j))
        ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
    for (unsigned int j = 0; j < js; ++j)
      std::swap(particles.at(ind23.at(b).at(ii23.at(b) - j)),
                particles.at(ind23.at(b).at(GL1 - j)));
    for (unsigned int i = init.at(b); i < fin.at(b); ++i) {
      if (particles.at(i).ir == 0)
        ind01.at(b).at(count0.at(b)--) = i;
      else if (particles.at(i).ir == 1)
        ind01.at(b).at(count1.at(b)++) = i;
      else if (particles.at(i).ir == 2)
        ind23.at(b).at(count2.at(b)--) = i;
      else
        ind23.at(b).at(count3.at(b)++) = i;
    }
  }

#ifdef OPENACC
#pragma acc parallel loop copy(GL1,count0,count1,count2,count3,mini,ii1,ii3) present(particles,ind01,ind23)
#else
#pragma omp parallel for
#pragma distribute_point
#endif
  for (unsigned int b = 0; b < Nbin; ++b) {
    ii1.at(b) = count1.at(b) - 1;
    mini.at(b) = GL1 - count0.at(b);
    if (count1.at(b) < mini.at(b))
      mini.at(b) = count1.at(b);
    unsigned int js = 0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+ : js)
#endif
    for (unsigned int j = 0; j < mini.at(b); ++j)
      if (ind01.at(b).at(ii1.at(b) - j) > ind01.at(b).at(GL1 - j))
        ++js;
#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
    for (unsigned int j = 0; j < js; ++j)
      std::swap(particles.at(ind01.at(b).at(ii1.at(b) - j)),
                particles.at(ind01.at(b).at(GL1 - j)));
    ii3.at(b) = count3.at(b) - 1;
    mini.at(b) = GL1 - count2.at(b);
    if (count3.at(b) < mini.at(b))
      mini.at(b) = count3.at(b);
    js = 0;
#ifdef OPENACC
#pragma acc loop vector reduction(+:js)
#else
#pragma omp simd reduction(+ : js)
#endif
    for (unsigned int j = 0; j < mini.at(b); ++j)
      if (ind23.at(b).at(ii3.at(b) - j) > ind23.at(b).at(GL1 - j))
        ++js;

#ifdef OPENACC
#pragma acc loop vector
#else
#pragma omp simd
#endif
    for (unsigned int j = 0; j < js; ++j)
      std::swap(particles.at(ind23.at(b).at(ii3.at(b) - j)),
                particles.at(ind23.at(b).at(GL1 - j)));
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
    count0.at(b) = GL1 - count0.at(b);
    count2.at(b) = GL1 - count2.at(b);
    POSITION0 += count0.at(b);
    POSITION1_ += count1.at(b);
    POSITION2_ += count2.at(b);
    POSITION3_ += count3.at(b);
    POSITION23_ += count23.at(b);
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
  pointer1.at(0) = pointer2.at(0) = pointer3.at(0) = 0;
  for (unsigned int b = 0; b < Nbin - 1; ++b) {
    pointer1.at(b + 1) = pointer1.at(b) + count1.at(b);
    pointer2.at(b + 1) = pointer2.at(b) + count2.at(b);
    pointer3.at(b + 1) = pointer3.at(b) + count3.at(b);
  }
#ifdef OPENACC
}
#endif

  // DO NOT parallelize or vectorize - undefined behavior
  for (unsigned int b = 0; b < Nbin; ++b) {
    auto &src1 = arr1.at(pointer1.at(b));
    auto const newsize =
        init.at(b) + count3.at(b) + count2.at(b) + count1.at(b) + 1;

#ifdef OPENACC
    #pragma acc host_data use_device(particles,arr1,arr2,arr3)
    {
      cudaMemcpy(&src1, &particles.at(init.at(b) + count3.at(b) + count2.at(b)), count1.at(b) * sizep, cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr2.at(pointer2.at(b)),
                   &particles.at(init.at(b) + count3.at(b)),
                   count2.at(b) * sizep, cudaMemcpyDeviceToDevice);
      cudaMemcpy(&arr3.at(pointer3.at(b)), &particles.at(init.at(b)),
                   count3.at(b) * sizep, cudaMemcpyDeviceToDevice);
    }
#else
    memcpy(&src1, &particles.at(init.at(b) + count3.at(b) + count2.at(b)), count1.at(b) * sizep);
    memcpy(&arr2.at(pointer2.at(b)),
                 &particles.at(init.at(b) + count3.at(b)),
                 count2.at(b) * sizep);
    memcpy(&arr3.at(pointer3.at(b)), &particles.at(init.at(b)),
                 count3.at(b) * sizep);
#endif

  }

#ifdef OPENACC
#pragma acc host_data use_device(particles,arr1,arr2,arr3)
 {
   cudaMemcpy(&particles.at(0), &arr3.at(0), POSITION3 * sizep, cudaMemcpyDeviceToDevice);
   cudaMemcpy(&particles.at(POSITION3), &arr2.at(0), POSITION2 * sizep, cudaMemcpyDeviceToDevice);
   cudaMemcpy(&particles.at(POSITION23), &arr1.at(0), POSITION1 * sizep, cudaMemcpyDeviceToDevice);
 }
#else
 memcpy(&particles.at(0), &arr3.at(0), POSITION3 * sizep);
 memcpy(&particles.at(POSITION3), &arr2.at(0), POSITION2 * sizep);
 memcpy(&particles.at(POSITION23), &arr1.at(0), POSITION1 * sizep);
#endif

 if(!histogram)
 {
#ifdef OPENACC
#pragma acc parallel loop present(particles)
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (unsigned int i = 0; i < POSITION3; ++i) {
    unsigned int i2 = i + i;
    auto &particlei = particles.at(i);
    auto &newParticle = particles.at(i + SHIFT);
    newParticle = particlei;
    particlei.rs.seed(static_cast<RandGen::result_type>(MAX_ELEMENT + i2));
    newParticle.rs.seed(
        static_cast<RandGen::result_type>(MAX_ELEMENT + i2 + 1));
    if (particlei.pdg == aParticleTable.makePDGfromZandA(1,2)) // inc d
    {
      particlei.pdg = aParticleTable.makePDGfromZandA(1,2);
      newParticle.pdg = aParticleTable.makePDGfromZandA(1,2);
    }
    else if (particlei.pdg == aParticleTable.makePDGfromZandA(22,48))// inc Ti48
    {
      particlei.pdg = aParticleTable.makePDGfromZandA(22,48);
      newParticle.pdg = aParticleTable.makePDGfromZandA(22,48);
    }//
    else if (particlei.pdg == t3::PDG_t(22)) // gamma PDG=22
    {
      particlei.pdg = t3::PDG_t(11);
      newParticle.pdg = t3::PDG_t(-11);
    } else if (particlei.pdg == t3::PDG_t(11)) // electron PDG=11
    {
      particlei.pdg = t3::PDG_t(22);
      newParticle.pdg = t3::PDG_t(11);
    } else if (particlei.pdg == t3::PDG_t(-11)) // positron PDG=-11
    {
      particlei.pdg = t3::PDG_t(22);
      newParticle.pdg = t3::PDG_t(22);
    } else if (particlei.pdg == t3::PDG_t(2112)) // neutrons PDG=2112
    {
      particlei.pdg = t3::PDG_t(2112);
      newParticle.pdg = t3::PDG_t(2112);
    }
  }
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(LIFE, MAX_ELEMENT) present(particles)
{
#endif
  MAX_ELEMENT += POSITION3 + POSITION3;
#ifdef OPENACC
}
#endif

#ifdef OPENACC
#pragma acc parallel loop copy(POSITION3, POSITION23, SHIFT, MAX_ELEMENT) present(particles)
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (unsigned int i = POSITION3; i < POSITION23; ++i) {
    auto const &particlei = particles.at(i);
    auto &newParticle = particles.at(i + SHIFT);
    newParticle = particlei;
    newParticle.rs.seed(static_cast<RandGen::result_type>(MAX_ELEMENT + i));
    if (particles.at(i).pdg == aParticleTable.makePDGfromZandA(1,2)) // inc d
      newParticle.pdg = aParticleTable.makePDGfromZandA(1,2);
    else if (particles.at(i).pdg == aParticleTable.makePDGfromZandA(22,48))// inc Ti48
      newParticle.pdg = aParticleTable.makePDGfromZandA(22,48);
    //
    else if (particles.at(i).pdg == t3::PDG_t(22)) //gamma PDG=22
      newParticle.pdg = t3::PDG_t(11); //new particle - electron PDG=11
    else
      newParticle.pdg = t3::PDG_t(22); //new particle - gamma PDG=22
  }
#ifdef OPENACC
#pragma acc parallel num_gangs(1) vector_length(1) copy(LIFE,MAX_ELEMENT) present(particles) 
{
#endif
  MAX_ELEMENT += POSITION23 - POSITION3;
#ifdef OPENACC
}
#endif

 }//end of (!histogram)
}

template <typename Floating> void DataHolder<Floating>::Propagate() {

  static std::array<Floating, Np> csIsotropic;
  static std::array<Floating, Np> csDuplicate;
  static std::array<Floating, Np> csRutherford;
  static std::array<Floating, Np> csMultipleScattering;
  
  auto inputlskinE = t3::makeAOSMemberView(
      particles.data(),
      [this](Particle<Floating> const &particle)  {//[this] is necessary to access aParticleTable.
        return particle.GetEtot() - aParticleTable.GetMass(particle.pdg);
      });
  auto inputPDG = t3::makeAOSMemberView(
      particles.data(),
      [](Particle<Floating> const &particle)  {
        return particle.pdg;
      });

  multiplescatteringProcess.GetCS(inputlskinE, inputPDG, t3::MatID_t(2u), LIFE,
                                   csMultipleScattering.data());

  constexpr Floating da = ag * 1e-10;
  constexpr Floating rho0 = 1.0;
  constexpr auto lMAX = std::numeric_limits<Floating>::max();
  Floating const ntid2 = aMaterialTable.GetConcentrations(t3::MatID_t(2u));
  Floating const rho_tid2 = aMaterialTable.GetDensity(t3::MatID_t(2u));
  Floating DE=0.0;
  Floating AVDE=0.0;

#ifdef OPENACC
#pragma acc parallel loop copy(da, rho0, lMAX, ntid2, rho_tid2) present(particles, LIFE)
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (size_t i = 0; i < LIFE; ++i)
  {
    auto &particlei = particles.at(i);
    auto En = particlei.GetEtot();

    if (particlei.ir > 0)
    {
      auto const csDuplicateI = csDuplicate.at(i);
      auto const csRutherfordi = csRutherford.at(i);
      auto const csMultipleScatteringi = csMultipleScattering.at(i);

      Floating  l1x =
          (particlei.vx() == 0.)
              ? lMAX
              : ((particlei.vx() > 0.)
                     ? ((particlei.ix() + 1) * ag - particlei.r.x())
                     : (particlei.ix() * ag - particlei.r.x())) /
                        particlei.vx() +
                    da;
      Floating  l1y =
          (particlei.vy() == 0.)
              ? lMAX
              : ((particlei.vy() > 0.)
                     ? ((particlei.jy() + 1) * ag - particlei.r.y())
                     : (particlei.jy() * ag - particlei.r.y())) /
                        particlei.vy() +
                    da;
      Floating  l1z =
          (particlei.vz() == 0.)
              ? lMAX
              : ((particlei.vz() > 0.)
                     ? ((particlei.kz() + 1) * ag - particlei.r.z())
                     : (particlei.kz() * ag - particlei.r.z())) /
                        particlei.vz() +
                    da;

      Floating const dEdx = 2.0 * MeV * cm * cm / gr;
      Floating const dEdxfull = dEdx*rho_tid2;
      Floating const range = (particlei.GetEtot() - aParticleTable.GetMass(particlei.pdg))/dEdxfull;
      Floating l0 = range;
      l0=1.0*cm;
      Floating const csMScatteringi = csMultipleScatteringi;
      Floating const lambdar = 1.0/ntid2/csMScatteringi;
      auto const R = particlei.GenerateCanonical();
      Floating l2 = std::abs(lambdar * std::log(R));
      Floating const l1 = std::min({l1x, l1y, l1z});
      int const irc = (l0 < l2 && l0 < l1) ? 0 : ((l2 < l1) ? 2 : 1);
      Floating const l = std::min({l1, l2, l0});
      Floating const dl = std::abs(l);
      auto const indexz = particlei.kz();

      if(irc == 1 && std::abs(l-l1z)<da)
      {
          auto const modparticlei = particlei.p.R();
          auto const costheta = particlei.p.z()/modparticlei;
          auto const theta = std::acos(costheta);
          const int id = particlei.id;
          auto const thetax=std::atan( particlei.p.x()/particlei.p.z() );//=px/pz.
          auto const thetay=std::atan( particlei.p.y()/particlei.p.z() );//=py/pz.
          auto const x=particlei.r.x() - Particle<Floating>::InitParticlex0 * Particle<Floating>::ag;
          auto const y=particlei.r.y() - Particle<Floating>::InitParticley0 * Particle<Floating>::ag;
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
              ++HTheta.at(indexz+cuba).at(ThetaBin);
          }
//-------------------------------------------
          const int ThetaxBin=(thetax-HThetaxMin)/deltaThetax;
          if(ThetaxBin>=0 && ThetaxBin<BinNumber2)
#ifdef OPENACC
            #pragma acc atomic update
#else
            #pragma omp atomic update
#endif
              ++HThetax.at(indexz+cuba).at(ThetaxBin);          
//-------------------------------------------
          const int ThetayBin=(thetay-HThetayMin)/deltaThetay;
          if(ThetayBin>=0 && ThetayBin<BinNumber2)
#ifdef OPENACC
            #pragma acc atomic update
#else
            #pragma omp atomic update
#endif
              ++HThetay.at(indexz+cuba).at(ThetayBin);                   
//*/              
      }      
      particlei.r.SetPxPyPzE (particlei.r.x()+particlei.vx() * (dl + da),
                              particlei.r.y()+particlei.vy() * (dl + da),
                              particlei.r.z()+particlei.vz() * (dl + da), 0.0);//t=0.0

      bool const out = (particlei.ix() >= cuba || particlei.jy() >= cuba ||
                        particlei.kz() >= cuba || particlei.ix() < -cuba ||
                        particlei.jy() < -cuba || particlei.kz() < -cuba);
      Floating loss = 0.;   
      if (aParticleTable.IsNucleus(particlei.pdg))
      {
        loss = dEdxfull * dl;
        if (loss >= particlei.GetEtot() - aParticleTable.GetMass(particlei.pdg))
        {
          particlei.de += particlei.GetEtot() - aParticleTable.GetMass(particlei.pdg);
          particlei.SetEtot(aParticleTable.GetMass(particlei.pdg), aParticleTable);
          particlei.ir = 0;
        }
        else
        {
          Floating const old_energy = particlei.GetEtot();
          Floating const new_energy = old_energy - loss;
          particlei.SetEtot(new_energy, aParticleTable);
          bool check=false;
          bool cond=indexz*ag<particlei.r.z() && particlei.r.z()<=(indexz+1)*ag;
          particlei.de += loss;
          particlei.ir = irc;
        }
      }  
      if (out)
      {
        particlei.ir = 0;
      }
    }//End if ir>0
  }//End of i-particle loop
  AVDE/=Np;
}

template <typename Floating> void DataHolder<Floating>::React() {
  
  auto start1=std::chrono::steady_clock::now();
  auto returnZeroMat = t3::makeConstantView(t3::MatID_t(2));
  auto inputP = t3::makeAOSMemberView(
      particles.data(), [](Particle<Floating> const &particle) {
        return particle.p;
      });
  auto inputPDG = t3::makeAOSMemberView(
      particles.data(), [](Particle<Floating> const &particle) {
        return particle.pdg;
      });
  auto inputRNG = t3::makeAOSMemberView(
      particles.data(),
      [](Particle<Floating> &particle) -> decltype(particle.rs) & {
        return particle.rs;
      });
  auto const outPDGoutput = std::make_tuple(outPDG1.data(), outPDG2.data());
  auto const outPoutput = std::make_tuple(outP1.data(), outP2.data());

  multiplescatteringProcess.GetFS(inputP, inputPDG, returnZeroMat, inputRNG,
                         POSITION23, outPDGoutput, outPoutput);
  auto start2=std::chrono::steady_clock::now();
  auto end2=std::chrono::steady_clock::now();
#ifdef OPENACC
#pragma acc parallel loop
#else
#pragma omp parallel for simd schedule(dynamic)
#pragma distribute_point
#endif
  for (unsigned int i = 0; i < POSITION23; ++i) {
    auto const &P1i = outP1.at(i);
    auto const pdg1 = outPDG1.at(i);
    auto &particlei = particles.at(i);
    particlei.p = P1i;
    particlei.ir = 1;
    particlei.pdg = pdg1;
  }
  auto end3=std::chrono::steady_clock::now();
}

template <typename Floating> void DataHolder<Floating>::Inject() {
  auto begin = std::chrono::steady_clock::now();
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
    auto sgi = particles.at(j).GetEtot() * static_cast<double>(particles.at(j).wt)-sgCompensation;
    auto sgtemp = sg + sgi;
    sgCompensation = (sgtemp - sg) - sgi;
    sg = sgtemp;
    auto wti = static_cast<double>(particles.at(j).wt) - wtCompensation;
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
  auto end = std::chrono::steady_clock::now();
  auto elapsed_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
 }
} // end of namespace dataholder

#include "T3MaterialTable.h"

#endif // T3DATAHOLDER_H
