#ifndef presidual_hpp
#define presidual_hpp

#include "control.hpp"

namespace mech {

template <typename T> struct Displacement;
template <typename T> struct Kinematics;
template <typename T> struct Model;
template <typename T> struct Pressure;

template <typename T>
struct PResidual : public Integrator {

  PResidual(
      RCP<Integrator> u,
      RCP<Integrator> p,
      RCP<Kinematics<T>> k,
      RCP<Model<T>> m,
      Disc* d,
      Input* in);

  void set_elem_set(std::string const& p);
  void at_point(Vector const&, double ipw, double dv);
  void do_small_strain(double ipw, double dv);
  void do_large_strain(double ipw, double dv);

  RCP<Displacement<T>> u;
  RCP<Pressure<T>> p;
  RCP<Kinematics<T>> kin;
  RCP<Model<T>> model;

  int dim;
  double mu;
  double lambda;
  double kappa;

  Disc* disc;
  Input* input;

};

}

#endif
