#ifndef first_pk_hpp
#define first_pk_hpp

#include "control.hpp"
#include <MiniTensor.h>

namespace apf {
class MeshEntity;
}

namespace mech {

template <typename T> struct Model;
template <typename T> struct Pressure;
template <typename T> struct Kinematics;

template <typename T>
struct FirstPK : public Integrator {

  using TensorT = minitensor::Tensor<T>;

  FirstPK(
      RCP<Model<T>> m,
      RCP<Integrator> k,
      RCP<Integrator> p,
      Disc* d,
      bool save);

  void in_elem(apf::MeshElement* me);
  void at_point(Vector const&, double, double);
  void do_small_strain();
  void do_large_strain();
  void out_elem();

  TensorT& val() { return first_pk; }

  int dim;
  int ip;
  bool save_states;
  Disc* disc;

  TensorT first_pk;
  TensorT I;

  RCP<Model<T>> model;
  RCP<Kinematics<T>> kin;
  RCP<Pressure<T>> pressure;

  apf::MeshEntity* elem;

};

}

#endif
