#include "disc.hpp"
#include "displacement.hpp"
#include "elastic.hpp"

namespace mech {

template <typename T>
Elastic<T>::Elastic(RCP<Integrator> disp, Disc* d, Input* in) {
  u = rcp_static_cast<Displacement<T>>(disp);
  disc = d;
  input = in;
  dim = disc->dim;
  sigma.set_dimension(dim);
  eps.set_dimension(dim);
  I = minitensor::eye<T>(dim);
}

template <typename T>
void Elastic<T>::set_elem_set(std::string const& set) {
  assert(input->mats.count(set));
  auto& mat = input->mats[set];
  double E = mat.E;
  double nu = mat.nu;
  mu = E / (2.0 * (1.0 + nu));
  lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

template <typename T>
void Elastic<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
  sigma_elem = apf::createElement(disc->sigma, me);
}

template <typename T>
void Elastic<T>::at_point(Vector const&, double, double) {
  for (int i = 0; i < dim; ++i)
  for (int j = 0; j < dim; ++j)
    eps(i,j) = 0.5 * (u->grad(i,j) + u->grad(j,i));
  sigma = 2.0*mu*eps + lambda*minitensor::trace(eps)*I;
}

template <typename T>
void Elastic<T>::out_elem() {
  apf::destroyElement(sigma_elem);
  elem = 0;
}

template struct Elastic<ST>;
template struct Elastic<FADT>;

}
