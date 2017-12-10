#ifndef model_hpp
#define model_hpp

#include "control.hpp"
#include <MiniTensor.h>

namespace apf {
class Element;
class MeshEntity;
}

namespace mech {

template <typename T>
struct Model : public Integrator {
  using TensorT = minitensor::Tensor<T>;
  Model();
  virtual ~Model();
  virtual void set_elem_set(std::string const&) {}
  virtual void in_elem(apf::MeshElement*) {}
  virtual void at_point(Vector const&, double, double) {}
  virtual TensorT const& get_sigma_dev() = 0;
  virtual bool small_strain() = 0;
};

}

#endif
