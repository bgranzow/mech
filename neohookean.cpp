#include "disc.hpp"
#include "displacement.hpp"
#include "kinematics.hpp"
#include "neohookean.hpp"

namespace mech {

template <typename T>
Neohookean<T>::Neohookean(RCP<Kinematics<T>> k, Disc* d, Input* in) {
  kin = k;
  disc = d;
  input = in;
  dim = disc->dim;
  b.set_dimension(dim);
  sigma_dev.set_dimension(dim);
}

template <typename T>
void Neohookean<T>::set_elem_set(std::string const& set) {
  assert(input->mats.count(set));
  auto& mat = input->mats[set];
  double E = mat.E;
  double nu = mat.nu;
  mu = E / (2.0 * (1.0 + nu));
  kappa = E / (3.0 * (1.0 - 2.0 * nu));
}

template <typename T>
void Neohookean<T>::at_point(Vector const&, double, double) {
  auto F = kin->F;
  auto J = kin->J;
  T Jm13 = 1.0 / std::cbrt(J);
  T Jm23 = Jm13 * Jm13;
  T Jm53 = Jm23 * Jm23 * Jm13;
  b = F*minitensor::transpose(F);
  sigma_dev = mu * Jm53 * minitensor::dev(b);
}

template struct Neohookean<ST>;
template struct Neohookean<FADT>;

}
