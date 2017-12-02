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
  void in_elem(apf::MeshElement* me);
  void at_point(Vector const&, double, double);
  void out_elem();

  TensorT& get_cauchy() { return sigma; }
  TensorT& get_first_pk() { return sigma; }

  Disc* disc;
  Input* input;
  RCP<Displacement<T>> u;

  int dim;
  double mu;
  double lambda;

  TensorT sigma;
  TensorT eps;
  TensorT I;

  apf::MeshEntity* elem;
  apf::Element* sigma_elem;

};

}

#endif
