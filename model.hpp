#ifndef model_hpp
#define model_hpp

#include "control.hpp"
#include <MiniTensor.h>

namespace mech {

template <typename T>
struct Model : public Integrator {
  using TensorT = minitensor::Tensor<T>;
  Model();
  virtual ~Model();
  virtual void set_elem_set(std::string const&) {}
  virtual void in_elem(apf::MeshElement*) {}
  virtual void at_point(Vector const&, double, double) {}
  virtual TensorT& get_cauchy() = 0;
  virtual TensorT& get_first_pk() = 0;
};

}

#endif
