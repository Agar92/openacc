#pragma once

#ifndef T3LORENTZVECTOR_H
#define T3LORENTZVECTOR_H

#include <cmath>
#include "T3ThreeVector.h"

namespace t3
{

//************************************************************//
//These forward declarations are necessary for defining friend//
//functions in template class LorentzVector<T>:               //
//************************************************************//
template <typename T>
class LorentzVector;

template <typename T>
LorentzVector<T> operator*(const T a, const LorentzVector<T> & rhsvector);

template <typename T>
LorentzVector<T> operator*(const LorentzVector<T> & lhsvector, const T a);
/*
template <typename T>
LorentzVector<T> operator*(const int a, const LorentzVector<T> & rhsvector);

template <typename T>
LorentzVector<T> operator*(const LorentzVector<T> & lhsvector, const int a);
*/
//***************************************************************//
//These combinations of symbols may be reserved and used by the compiler.
// __x __X _X (x, X - any symbol) - everywhere
// _x (x -any symbol)             - in global scope
//So do not use anywhere in the program __x __X _X (x - any variable)
//and do not create global variable _x (x - any symbol). Or if the
//compiler may rewrite them with its varibles, and an error will be.
//***************************************************************//

template <typename T>
class LorentzVector
{
private:
//(x, y, z, t) or (px, py, pz, E).  
  T fx, fy, fz, fE;
public:
//contsructors:
//\\//LorentzVector()=delete;
  LorentzVector():fx(0.0), fy(0.0), fz(0.0), fE(0.0){}
  LorentzVector(T x, T y, T z, T E):fx(x),fy(y),fz(z),fE(E){}
  LorentzVector(LorentzVector const & vector);
//destructor:
  ~LorentzVector(){}
//return components and vector lenghth:
//getters:  
  T x() const {return fx;}
  T y() const {return fy;}
  T z() const {return fz;}
  T t() const {return fE;}
  T E() const {return fE;}
//returns mass (m=SQRT(E^2-p^2)):  
  T mass() const {return std::sqrt(fE*fE-fx*fx-fy*fy-fz*fz);}
//returns energy:  
  T energy() const {return fE;}
//returns the absolute value of 4-vector:  
  T R() const {return std::sqrt(fx*fx+fy*fy+fz*fz);}
//returns mass (m=SQRT(E^2-p^2)):  
  T m() const {return std::sqrt(fE*fE - fx*fx - fy*fy - fz*fz);}
//makes 3-vector from 4-vector:
  ThreeVector<T> Vect() const
  {
    ThreeVector<T> vector3=ThreeVector<T>(this->fx, this->fy, this->fz);
    return vector3;
  }
//setters:
  void SetX(T X){fx=X;}
  void SetY(T Y){fy=Y;}
  void SetZ(T Z){fz=Z;}
  void SetE(T E){fE=E;}
  void SetPxPyPzE(T x, T y, T z, T E)
  {
    fx=x, fy=y, fz=z, fE=E;
  }
//operator overloading:  
  LorentzVector & operator=(const LorentzVector & rhsvector);
  LorentzVector operator+(const LorentzVector & rhsvector);
  LorentzVector operator-(const LorentzVector & rhsvector);
  LorentzVector & operator+=(const LorentzVector & rhsvector);
  LorentzVector & operator-=(const LorentzVector & rhsvector);
  friend LorentzVector operator*<>(const T a, const LorentzVector<T> & rhsvector);
  friend LorentzVector operator*<>(const LorentzVector<T> & lhsvector, const T a);
//friend LorentzVector operator*<>(const int a, const LorentzVector<T> & rhsvector);
//friend LorentzVector operator*<>(const LorentzVector<T> & lhsvector, const int a);
  LorentzVector & operator*=(const T a);
//LorentzVector & operator*=(const int a);
  LorentzVector operator/(const T a);
//LorentzVector operator/(const int a);
  LorentzVector & operator/=(const T a);
//LorentzVector & operator/=(const int a);
};

template <typename T>
LorentzVector<T>::LorentzVector(LorentzVector const & vector):fx(vector.fx), fy(vector.fy), fz(vector.fz), fE(vector.fE){}

template <typename T>
LorentzVector<T> & LorentzVector<T>::operator=(const LorentzVector & rhsvector)
{
  if(this==&rhsvector)
  {
    return *this;
  }
  this->fx=rhsvector.fx, this->fy=rhsvector.fy, this->fz=rhsvector.fz, this->fE=rhsvector.fE;
  return *this;
}

template <typename T>
LorentzVector<T> LorentzVector<T>::operator+(LorentzVector const & rhsvector)
{
  LorentzVector<T> result(this->fx+rhsvector.fx,this->fy+rhsvector.fy,this->fz+rhsvector.fz,this->fE+rhsvector.fE);
  return result;
}

template <typename T>
LorentzVector<T> LorentzVector<T>::operator-(LorentzVector const & rhsvector)
{
  LorentzVector<T> result(this->fx-rhsvector.fx,this->fy-rhsvector.fy,this->fz-rhsvector.fz,this->fE-rhsvector.fE);
  return result;
}

template <typename T>
LorentzVector<T> & LorentzVector<T>::operator+=(LorentzVector const & rhsvector)
{
  this->fx+=rhsvector.fx, this->fy+=rhsvector.fy, this->fz+=rhsvector.fz,this->fE+=rhsvector.fE;
  return *this;
}

template <typename T>
LorentzVector<T> & LorentzVector<T>::operator-=(LorentzVector const & rhsvector)
{
  this->fx-=rhsvector.fx, this->fy-=rhsvector.fy, this->fz-=rhsvector.fz,this->fE-=rhsvector.fE;
  return *this;
}

template <typename T>
LorentzVector<T> operator*(const T a, const LorentzVector<T> & rhsvector)
{
  LorentzVector<T> result(rhsvector.fx*a, rhsvector.fy*a, rhsvector.fz*a, rhsvector.fE*a);
  return result;
}

template <typename T>
LorentzVector<T> operator*(const LorentzVector<T> & lhsvector, const T a)
{
  LorentzVector<T> result(lhsvector.fx*a, lhsvector.fy*a, lhsvector.fz*a, lhsvector.fE*a);
  return result;
}
/*
template <typename T>
LorentzVector<T> operator*(const int a, const LorentzVector<T> & rhsvector)
{
  LorentzVector<T> result(rhsvector.fx*a, rhsvector.fy*a, rhsvector.fz*a, rhsvector.fE*a);
  return result;
}
*/
/*
template <typename T>
LorentzVector<T> operator*(const LorentzVector<T> & lhsvector, const int a)
{
  LorentzVector<T> result(lhsvector.fx*a, lhsvector.fy*a, lhsvector.fz*a, lhsvector.fE*a);
  return result;
}
*/
template <typename T>
LorentzVector<T> & LorentzVector<T>::operator*=(const T a)
{
  this->fx*=a, this->fy*=a, this->fz*=a, this->fE*=a;
  return *this;
}
/*
template <typename T>
LorentzVector<T> & LorentzVector<T>::operator*=(const int a)
{
  this->fx*=a, this->fy*=a, this->fz*=a, this->fE*=a;
  return *this;
}
*/
template <typename T>
LorentzVector<T> LorentzVector<T>::operator/(const T a)
{
  LorentzVector<T> result(this->fx/a, this->fy/a, this->fz/a, this->fE/a);
  return result;
}
/*
template <typename T>
LorentzVector<T> LorentzVector<T>::operator/(const int a)
{
  LorentzVector<T> result(this->fx/a,this->fy/a,this->fz/a,this->fE/a);
  return result;
}
*/
template <typename T>
LorentzVector<T> & LorentzVector<T>::operator/=(const T a)
{
  this->fx/=a, this->fy/=a, this->fz/=a, this->fE/=a;
  return *this;
}
/*
template <typename T>
LorentzVector<T> & LorentzVector<T>::operator/=(const int a)
{
  this->fx/=a, this->fy/=a, this->fz/=a, this->fE/=a;
  return *this;
}
*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const LorentzVector<T> & vector)
{
  os<<vector.x()<<" "<<vector.y()<<" "<<vector.z()<<" "<<vector.E()<<" ";
  return os;
}

}

#endif  //T3LORENTZVECTOR
