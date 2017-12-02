#ifndef kinematics_hpp
#define kinematics_hpp

#include "control.hpp"
#include <MiniTensor.h>

namespace mech {

template <typename T> struct Displacement;

template <typename T>
struct Kinematics : public Integrator {
  using TensorT = minitensor::Tensor<T>;
  Kinematics(RCP<Integrator> disp);
  void at_point(Vector const&, double, double);
  RCP<Displacement<T>> u;
  int dim;
  T J;
  TensorT F;
};

}

#endif
