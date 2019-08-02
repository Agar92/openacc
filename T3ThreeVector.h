#pragma once

#ifndef T3THREEVECTOR_H
#define T3THREEVECTOR_H

#include <cmath>

namespace t3
{

//************************************************************//
//These forward declarations are necessary for defining friend//
//functions in template class ThreeVector<T>:                 //
//************************************************************//
template <typename T>
class ThreeVector;

template <typename T>
ThreeVector<T> operator*(const T a, ThreeVector<T> const rhsside);

template <typename T>
ThreeVector<T> operator*(ThreeVector<T> const lhsside, const T a);
  
template <typename T>
class ThreeVector
{
private:
  T fx, fy, fz;
public:
//contsructors:
  ThreeVector():fx(0.0), fy(0.0), fz(0.0){}
  ThreeVector(T x, T y, T z):fx(x),fy(y),fz(z){}
  ThreeVector(ThreeVector const & vector);
//destructor:
  ~ThreeVector(){}
//return components and vector lenghth:
//getters:
  T x() const {return fx;}
  T y() const {return fy;}
  T z() const {return fz;}  
  T R() const {return std::sqrt(fx*fx+fy*fy+fz*fz);}
//normalizes 3-vector:
  ThreeVector<T>  Unit()
  {
    const T vector_length=std::sqrt(fx*fx+fy*fy+fz*fz);
    ThreeVector<T> result(this->fx/vector_length, this->fy/vector_length, this->fz/vector_length);
    return result;
  }
//setters:
  void SetX(T X){fx=X;}
  void SetY(T Y){fy=Y;}
  void SetZ(T Z){fz=Z;}
//operator overloading:  
  ThreeVector & operator=(ThreeVector const & rhsvector);
  ThreeVector operator+(ThreeVector const rhsvector);
  ThreeVector operator-(ThreeVector const rhsvector);
  ThreeVector operator+=(const ThreeVector rhsvector);
  ThreeVector operator-=(const ThreeVector rhsvector);
  friend ThreeVector operator*<>(const T a, ThreeVector const rhsside);
  friend ThreeVector operator*<>(ThreeVector const lhsside, const T a);
  ThreeVector operator*=(const T a);
  ThreeVector operator/(const T a);
  ThreeVector operator/=(const T a);
};

template <typename T>
ThreeVector<T>::ThreeVector(ThreeVector const & vector):fx(vector.fx), fy(vector.fy), fz(vector.fz){}

template <typename T>
ThreeVector<T> & ThreeVector<T>::operator=(const ThreeVector & rhsvector)
{
  if(this==&rhsvector)
  {
    return *this;
  }
  this->fx=rhsvector.fx, this->fy=rhsvector.fy, this->fz=rhsvector.fz;
  return *this;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator+(ThreeVector const rhsvector)
{
  ThreeVector<T> result(this->fx+rhsvector.fx,this->fy+rhsvector.fy,this->fz+rhsvector.fz);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator-(ThreeVector const rhsvector)
{
  ThreeVector<T> result(this->fx-rhsvector.fx,this->fy-rhsvector.fy,this->fz-rhsvector.fz);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator+=(ThreeVector const rhsvector)
{
  ThreeVector<Floating> result(this->fx+rhsvector.fx, this->fy+rhsvector.fy, this->fz+rhsvector.fz);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator-=(ThreeVector const rhsvector)
{
  ThreeVector<T> result(this->fx-rhsvector.fx, this->fy-rhsvector.fy, this->fz-rhsvector.fz);
  return result;
}

template <typename T>
ThreeVector<T> operator*(const T a, ThreeVector<T> const rhsside)
{
  ThreeVector<T> result(rhsside.fx*a, rhsside.fy*a, rhsside.fz*a);
  return result;
}

template <typename T>
ThreeVector<T> operator*(ThreeVector<T> const lhsside, const T a)
{
  ThreeVector<T> result(lhsside.fx*a, lhsside.fy*a, lhsside.fz*a);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator*=(const T a)
{
  ThreeVector<T> result(this->fx*a, this->fy*a, this->fz*a);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator/(const T a)
{
  ThreeVector<T> result(this->fx/a, this->fy/a, this->fz/a);
  return result;
}

template <typename T>
ThreeVector<T> ThreeVector<T>::operator/=(const T a)
{
  ThreeVector<T> result(this->fx=a, this->fy=a, this->fz=a);
  return result;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const ThreeVector<T> & vector)
{
  os<<vector.x()<<" "<<vector.y()<<" "<<vector.z()<<" ";
  return os;
}
  
}

#endif  //T3THREEVECTOR
