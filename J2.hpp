#ifndef J2_hpp
#define J2_hpp

#include "model.hpp"

namespace mech {

template <typename T> struct Kinematics;

template <typename T>
struct J2 : public Model<T> {

  using TensorT = minitensor::Tensor<T>;

  J2(RCP<Kinematics<T>> k, Disc* d, Input* in, bool save);
  void set_elem_set(std::string const& set);
  void in_elem(apf::MeshElement* me);
  void at_point(Vector const&, double, double);
  void out_elem();

  TensorT const& get_sigma_dev() { return sigma_dev; }
  bool small_strain() { return false; }

  Disc* disc;
  Input* input;
  RCP<Kinematics<T>> kin;

  int dim;
  int ip;
  double mu;
  double kappa;
  double K;
  double Y;
  double sq23;
  bool save;

  T Jm23;
  T mubar;
  T smag;
  T f;
  T dgam;
  T eqps;

  TensorT Fp;
  TensorT Fpinv;
  TensorT Cpinv;
  TensorT be;
  TensorT s;
  TensorT N;
  TensorT Fpn;
  TensorT sigma_dev;

  apf::Element* Fpold_elem;
  apf::Element* eqpsold_elem;
  apf::MeshEntity* elem;

};

}

#endif
