#pragma once
#ifndef T3MULTIPLESCATTERINGFSIMPL_H
#define T3MULTIPLESCATTERINGFSIMPL_H

#include "T3Defs.h"
#include "T3ThreeVector.h"
#include "T3LorentzVector.h"
#include "T3Util.h"
#include "T3ParticleTable.h"
#include "T3MaterialTable.h"
#include "T3MultipleScatteringCSImpl.h"
#include <random>

namespace t3 {
using namespace units;
template <bool generateRecoil = true, typename Floating = double> class MultipleScatteringFS {
public:
  template <typename RandomEngine>
  inline auto GetFS(LorentzVector<Floating> const &p, PDG_t incPDG,
                    MatID_t matID, RandomEngine &engine) const;
private:
  MultipleScatteringCS<Floating> aCS;
  ParticleTable aParticleTable;
  MaterialTable aMaterialTable;
};

template <bool generateRecoil, typename Floating>
template <typename RandomEngine>
auto MultipleScatteringFS<generateRecoil, Floating>::GetFS(
    LorentzVector<Floating> const &p, PDG_t incPDG, MatID_t matID,
    RandomEngine &engine) const {
  auto const generateSubCanonical = [&engine]() {
    return GenerateSubCanonical<Floating>(engine);
  };
  if (incPDG == PDG_t(22) || incPDG == PDG_t(2112))
  {
    if constexpr (!generateRecoil) return std::make_tuple(PDG_t(0), LorentzVector<Floating>(0.0,0.0,0.0,0.0));
    else return std::make_tuple(PDG_t(0), PDG_t(0), LorentzVector<Floating>(0.0,0.0,0.0,0.0), LorentzVector<Floating>(0.0,0.0,0.0,0.0));
  }
  //\\//static const auto aParticleTable = ParticleTable();
  //\\//static const auto aMaterialTable = MaterialTable();
  const PDG_t deuteronPDG = aParticleTable.makePDGfromZandA(1, 2);
  const PDG_t titanPDG    = aParticleTable.makePDGfromZandA(22, 48);
  auto const elsfull=p.energy();
  auto const m=p.mass();

  auto const Tls=elsfull-m;
  auto const numi = aMaterialTable.GetNumberOfIsotopes(matID);
  t3::vector<Floating> csBorder;
  csBorder.resize(numi+1);
  auto sigma_average = 0.0;
  csBorder.at(0) = 0.0;
  for(size_t i=0; i<numi; ++i)
  {
    sigma_average += aMaterialTable.GetFractions(matID).at(i)
            * aCS.GetCS(Tls, incPDG, aMaterialTable.GetIsotopes(matID).at(i));
    csBorder.at(i+1) = sigma_average;
  }
  auto const aR = generateSubCanonical();
  auto const raR = aR * sigma_average;
  PDG_t isotopePDG=PDG_t(-1);
  for(size_t i=0; i<numi; ++i)
  {
    if(csBorder.at(numi-1-i)<raR && raR<=csBorder.at(numi-i))
        isotopePDG = aMaterialTable.GetIsotopes(matID).at(numi-1-i);
  }
  auto tmin=0.0;
  Floating const M=aParticleTable.GetMass(isotopePDG);
  auto Ede=0.0;
  if(isotopePDG == deuteronPDG)
  {
    Ede=MultipleScatteringCS<Floating>::Edisplace_deuteron;
    tmin=2 * MultipleScatteringCS<Floating>::Edisplace_deuteron * M;
  }
  else if(isotopePDG == titanPDG)
  {
    Ede=MultipleScatteringCS<Floating>::Edisplace_titan;
    tmin=2 * MultipleScatteringCS<Floating>::Edisplace_titan * M;
  }
  auto const c=1.0;
  auto const alpha=1.0/137;
  auto const ec=std::sqrt(alpha);
  auto const me=aParticleTable.GetMass(PDG_t(11));
  auto const h=1.0;
  auto const a0=h*h/me/ec/ec;
  auto const hc=200.0 * MeV * fm;
  auto const CTF = 1.0/2.0 * std::pow(3*M_PI/4,2.0/3.0);
  auto const z=aParticleTable.GetZ(incPDG);
  auto const Z=aParticleTable.GetZ(isotopePDG);
  auto const s=m*m+2*elsfull*M+M*M;
  auto const mr=m*M/std::sqrt(s);
  auto const pls2=Tls*(2*m+Tls);
  auto const rat=m/M;
  auto const rat2=rat*rat;
  auto const pcm2=pls2/(1.0+rat2+2.0*std::sqrt(rat2+pls2/M/M));
  auto const pcm=std::sqrt(pcm2);
  auto const aI=CTF*a0/std::sqrt(std::pow(z, 2.0/3.0) + std::pow(Z, 2.0/3.0));
  auto const beta_r=1.0;
  auto const c12=alpha*z*Z/beta_r;
  auto const c122=c12*c12;
  auto const const1=h/2/aI/pcm;
  auto const const2=const1*const1;
  auto const AS=const2*(1.13+3.76*c122);
  ThreeVector<Floating> InitMomentum = p.Vect();
  auto const plsinc=InitMomentum.R();
  auto const Elsinc=std::sqrt(m*m+plsinc*plsinc);
  auto const coef=M/std::sqrt(s);
  ThreeVector<Floating> Pcm = InitMomentum*coef;
  auto costhetacm=0.0;
  auto const R=generateSubCanonical();
  auto const tmax_real = 4 * pcm2;
  auto tmax = 0.0;
  if(tmax_real<=tmin) tmax = tmax_real;
  else                tmax = tmin;
  auto const rc1=4*pcm2*AS;
  auto const rc2=tmax/(tmax+rc1);
  auto const t=rc1*R*rc2/(1.0-R*rc2);
  auto sinhalfthetacm = std::sqrt(t/4/pcm2);
  auto thetacm = 2 * std::asin(sinhalfthetacm);
  costhetacm = std::cos(thetacm);
  auto const xy=Pcm.x()*Pcm.y(), xz=Pcm.x()*Pcm.z(), yz=Pcm.y()*Pcm.z();
  auto const x2=Pcm.x()*Pcm.x(), y2=Pcm.y()*Pcm.y(), z2=Pcm.z()*Pcm.z();
  ThreeVector<Floating> e1, e2;
  if(Pcm.x() < Pcm.y())
  {
    if(Pcm.x() < Pcm.z())
    {
        e1={0., Pcm.z(), -Pcm.y()};
        e2={y2+z2, -xy, -xz};
    }
    else
    {
        e1={Pcm.y(), -Pcm.x(), 0.};
        e2={-xz, -yz, y2+x2};
    }
  }
  else
  {
    if(Pcm.y() < Pcm.z())
    {
        e1={Pcm.z(), 0., -Pcm.x()};
        e2={xy, -x2-z2, yz};
    }
    else
    {
        e1={Pcm.y(), -Pcm.x(), 0.};
        e2={-xz, -yz, y2+x2};
    }
  }
  e1=e1.Unit();
  e2=e2.Unit();
  auto const phi=generateSubCanonical() * TWOPI;
  auto const sinthetacm = std::sin(thetacm);
  auto const cosPhi = std::cos(phi);
  auto const sinPhi = std::sin(phi);
  auto const pcminc=Pcm.R();
  auto const ss=pcminc*sinthetacm*sinPhi, sc=pcminc*sinthetacm*cosPhi;
  ThreeVector<Floating> Pcmfin;
  Pcmfin.SetX(Pcm.x()*costhetacm+e1.x()*ss+e2.x()*sc);
  Pcmfin.SetY(Pcm.y()*costhetacm+e1.y()*ss+e2.y()*sc);
  Pcmfin.SetZ(Pcm.z()*costhetacm+e1.z()*ss+e2.z()*sc);
  auto const Ecminc=std::sqrt(m*m+pcminc*pcminc);
  auto const Vcm=InitMomentum/(Elsinc+M);
  auto const AbsVcm=Vcm.R();
  Floating const gamma=1.0/std::sqrt(1-AbsVcm*AbsVcm);
  ThreeVector<Floating> Plsfin1=gamma*(Pcmfin+Vcm*Ecminc);
  ThreeVector<Floating> Plsfin2=InitMomentum-Plsfin1;     
  auto const Absplsfin1=Plsfin1.R();
  auto const Elsfin1=std::sqrt(m*m+Absplsfin1*Absplsfin1);
  auto const Absplsfin2=Plsfin2.R();
  auto const Elsfin2=std::sqrt(M*M+Absplsfin2*Absplsfin2);
  auto const outPinc    = LorentzVector<Floating>(Plsfin1.x(), Plsfin1.y(), Plsfin1.z(), Elsfin1);
  auto const outPtarget = LorentzVector<Floating>(Plsfin2.x(), Plsfin2.y(), Plsfin2.z(), Elsfin2);
  auto const outPDGinc = incPDG;
  auto const outPDGtarget = isotopePDG;
  auto const Els1=p.E();
  auto const de=Elsfin1-Els1;
  if constexpr (!generateRecoil)
          return std::make_tuple(std::move(outPDGinc), std::move(outPinc));
  else
          return std::make_tuple(std::move(outPDGinc), std::move(outPDGtarget), std::move(outPinc), std::move(outPtarget));

  }
}// namespace t3

#endif // T3MULTIPLESCATTERINGFSIMPL_H
