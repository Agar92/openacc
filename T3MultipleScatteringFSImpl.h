#pragma once
#ifndef T3MULTIPLESCATTERINGFSIMPL_H
#define T3MULTIPLESCATTERINGFSIMPL_H

#include "T3Defs.h"
#include "T3LorentzVector.h"
#include "T3ParticleTable.h"
//#include "T3MaterialTable.h"
#include "T3MultipleScatteringCSImpl.h"

namespace t3 {
using namespace units;
//<typename Floating = double> this means that the default template type is double.
//If to write RutherfordCS<>(), the type of Floating will be double.
template <bool generateRecoil = true, typename Floating = double> class MultipleScatteringFS {
public:
  template <typename RandomEngine>
  inline auto GetFS(LorentzVector<Floating> const &p/*4-momentum of the inc particle*/, PDG_t incPDG,
                    MatID_t matID, RandomEngine &engine, Floating * csBorderDataFS,
                    int ind, const ParticleTable & aParticleTable,
                    const MaterialTable & aMaterialTable) const;
private:
  MultipleScatteringCS<Floating> aCS;
};

template <bool generateRecoil, typename Floating>
template <typename RandomEngine>
auto MultipleScatteringFS<generateRecoil, Floating>::GetFS(
    LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
    RandomEngine &engine, Floating * csBorderDataFS, int ind,
    const ParticleTable & aParticleTable,
    const MaterialTable & aMaterialTable) const {
  //generateSubCanonical() returns random number from 0.0 to 1.0.
  auto const generateSubCanonical = [&engine]() {
    return GenerateSubCanonical<Floating>(engine);
  };
  //if the particle is a gamma or a neutron, the Rutherford cross section is 0.0.
  /*if (incPDG == PDG_t(22) || incPDG == PDG_t(2112))
  {
    if constexpr (!generateRecoil) return Pair<Floating>(PDG_t(0), LorentzVector<Floating>(0.0,0.0,0.0,0.0));
    else return Four<Floating>(PDG_t(0), PDG_t(0), LorentzVector<Floating>(0.0,0.0,0.0,0.0), LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }*/  

//Get STEP of csBorderDataCS/FS:
  const int csBorderSTEP=aMaterialTable.GetcsBorderSTEP();

//new code:  
  const PDG_t deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
  const PDG_t titanPDG    = aParticleTable.makePDGfromZandA(22, 48);
  auto const elsfull=p.energy();
  auto const m=aParticleTable.GetMass(incPDG);
  auto const Tls=elsfull-m;
  auto const numi = aMaterialTable.GetNumberOfIsotopes(matID);
  auto sigma_average = 0.0;
  const size_t init=ind * csBorderSTEP;
  const size_t fin =init+numi;
  //\\//csBorder[init] = 0.0;
  csBorderDataFS[init] = 0.0;
  for(size_t i=0; i<numi; ++i)
  {
    sigma_average += aMaterialTable.GetFractions(matID,i)
      * aCS.GetCS(Tls, incPDG, aMaterialTable.GetIsotopes(matID,i), aParticleTable);
    //\\//csBorder[init+i+1] = sigma_average;
    csBorderDataFS[init+i+1] = sigma_average;
  }
  auto const aR = GenerateSubCanonical<Floating>(engine);
  auto const raR = aR * sigma_average;
  PDG_t isotopePDG=PDG_t(-1);
  for(size_t i=0; i<numi; ++i)
  {
    if(csBorderDataFS[init+numi-1-i]<raR && raR<=csBorderDataFS[init+numi-i])
      isotopePDG = aMaterialTable.GetIsotopes(matID,numi-1-i);
  }
  auto tmin=0.0;
  Floating const M=aParticleTable.GetMass(isotopePDG);
  auto Ede=0.0;
  if(isotopePDG == deuteronPDG)   tmin=2 * aCS.Edisplace_deuteron * M;
  else if(isotopePDG == titanPDG) tmin=2 * aCS.Edisplace_titan * M;
//end of new code.

  
//Old code:
  float part=generateSubCanonical();
  float gamma1;
  float gamma2;    
  float gamma=p.E();
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
  ThreeVector<Floating> e1;
  ThreeVector<Floating> e2;
  float xy=p.x()*p.y();
  float xz=p.x()*p.z();
  float yz=p.y()*p.z();
  float x2=p.x()*p.x();
  float y2=p.y()*p.y();
  float z2=p.z()*p.z();
  if(p.x() < p.y())
  {
    if(p.x() < p.z())
    {
      e1=ThreeVector<Floating>(0., p.z(), -p.y());
      e2=ThreeVector<Floating>(y2+z2, -xy, -xz);
    }
    else
    {
      e1=ThreeVector<Floating>(p.y(), -p.x(), 0.);
      e2=ThreeVector<Floating>(-xz, -yz, y2+x2);
    }
  }
  else
  {
    if(p.y() < p.z())
    {
      e1=ThreeVector<Floating>(p.z(), 0., -p.x());
      e2=ThreeVector<Floating>(xy, -x2-z2, yz);
    }
    else
    {
      e1=ThreeVector<Floating>(p.y(), -p.x(), 0.);
      e2=ThreeVector<Floating>(-xz, -yz, y2+x2);
    }
  }
  e1.Unit();
  e2.Unit();
  float phi=generateSubCanonical()*2*M_PI;
  // first particle
  float Theta1=1.;
  if(gamma>Dan) Theta1=Dan/gamma;
  float cost1=cosf(Theta1);
  float sint1=0.;
  if(cost1<1.f&&cost1>-1.f) sint1=sqrtf(1.f-cost1*cost1);
  ThreeVector<Floating> v1;
  float sf=sinf(phi);
  float cf=cosf(phi);
  float ss=sint1*sf;
  float sc=sint1*cf;
  float P=p.E();
  v1.SetX(p.x()/P*cost1+e1.x()*ss+e2.x()*sc);
  v1.SetY(p.y()/P*cost1+e1.y()*ss+e2.y()*sc);
  v1.SetZ(p.z()/P*cost1+e1.z()*ss+e2.z()*sc);
  v1.Unit();
  LorentzVector<Floating> p1;
  //p.SetE(gamma1);???
  p1.SetE(gamma1);
  float Theta2=1.;
  if(gamma>Dan) Theta2=Dan/gamma;
  float cost2=cosf(Theta2);
  float sint2=0.;
  if(cost2<1.f&&cost2>-1.f) sint2=sqrtf(1.f-cost2*cost2);
  ThreeVector<FloatingType> v2;
  ss=sint2*sf;
  sc=sint2*cf;
  v2.SetX(p.x()/P*cost2-e1.x()*ss-e2.x()*sc);
  v2.SetY(p.y()/P*cost2-e1.y()*ss-e2.y()*sc);
  v2.SetZ(p.z()/P*cost2-e1.z()*ss-e2.z()*sc);
  v2.Unit();
  p1.SetX(v1.x()*p1.E());
  p1.SetY(v1.y()*p1.E());
  p1.SetZ(v1.z()*p1.E());
  LorentzVector<Floating> p2;
  p2.SetE(gamma2);
  p2.SetX(v2.x()*p2.E());
  p2.SetY(v2.y()*p2.E());
  p2.SetZ(v2.z()*p2.E());
  
  const auto outPDGinc=incPDG;
  const auto outPDGtarget=incPDG;//??????
  const auto outPinc=p1;
  const auto outPtarget=p2;
  
  if constexpr (!generateRecoil)
        return Pair<Floating>(outPDGinc, outPinc);
  else
        return Four<Floating>(outPDGinc, outPDGtarget, outPinc, outPtarget);

  }
  
}// namespace t3
#endif // T3MULTIPLESCATTERINGFSIMPL_H
