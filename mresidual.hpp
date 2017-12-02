#ifndef mresidual_hpp
#define mresidual_hpp

#include "control.hpp"

namespace mech {

template <typename T> struct FirstPK;
template <typename T> struct Displacement;

template <typename T>
struct MResidual : public Integrator {

  MResidual(RCP<FirstPK<T>> P, RCP<Integrator> u, Disc* d);
  void at_point(Vector const&, double ipw, double dv);

  RCP<Displacement<T>> u;
  RCP<FirstPK<T>> first_pk;

  int dim;
  Disc* disc;

};

}

#endif
