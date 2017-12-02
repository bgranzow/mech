#ifndef elastic_hpp
#define elastic_hpp

#include "model.hpp"

namespace mech {

template <typename T> struct Displacement;

template <typename T>
struct Elastic : public Model<T> {

  using TensorT = minitensor::Tensor<T>;

  Elastic(RCP<Integrator> disp, Disc* disc, Input* in);
  void set_elem_set(std::string const& set);
  void at_point(Vector const&, double, double);

  TensorT& get_sigma_dev() { return sigma_dev; }
  bool small_strain() { return true; }

  Disc* disc;
  Input* input;
  RCP<Displacement<T>> u;

  int dim;
  double mu;
  double lambda;

  TensorT sigma;
  TensorT sigma_dev;
  TensorT eps;
  TensorT I;

};

}

#endif
