#ifndef neohookean_hpp
#define neohookean_hpp

#include "model.hpp"

namespace mech {

template <typename T> struct Kinematics;

template <typename T>
struct Neohookean : public Model<T> {

  using TensorT = minitensor::Tensor<T>;

  Neohookean(RCP<Kinematics<T>> k, Disc* d, Input* in);
  void set_elem_set(std::string const& set);
  void at_point(Vector const&, double, double);

  TensorT const& get_sigma_dev() { return sigma_dev; }
  bool small_strain() { return false; }

  Disc* disc;
  Input* input;
  RCP<Kinematics<T>> kin;

  int dim;
  double mu;
  double kappa;

  TensorT b;
  TensorT sigma_dev;

};

}

#endif
