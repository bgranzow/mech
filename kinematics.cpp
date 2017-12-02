#include "displacement.hpp"
#include "kinematics.hpp"

namespace mech {

template <typename T>
Kinematics<T>::Kinematics(Integrator* disp) {
  u = static_cast<Displacement<T>*>(disp);
  dim = u->dim;
  J = 0.0;
  F.set_dimension(dim);
  this->name = "kinematics";
}

template <typename T>
void Kinematics<T>::at_point(Vector const&, double, double) {
  for (int i = 0; i < dim; ++i)
  for (int j = 0; j < dim; ++j)
    F(i, j) = u->grad(i, j);
  for (int i = 0; i < dim; ++i)
    F(i, i) += 1.0;
  J = minitensor::det(F);
}

template struct Kinematics<ST>;
template struct Kinematics<FADT>;

}
