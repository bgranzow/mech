#ifndef kinematics_hpp
#define kinematics_hpp

#include "control.hpp"
#include <MiniTensor.h>

namespace mech {

template <typename T> struct Displacement;

template <typename T>
struct Kinematics : public Integrator {
  using TensorT = minitensor::Tensor<T>;
  Kinematics(Integrator* disp);
  void at_point(Vector const&, double, double);
  T J;
  TensorT F;
  int dim;
  Displacement<T>* u;
};

}

#endif
