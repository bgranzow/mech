#include "displacement.hpp"
#include "kinematics.hpp"
#include "model.hpp"
#include "pressure.hpp"
#include "presidual.hpp"

namespace mech {

template <typename T>
PResidual<T>::PResidual(
    RCP<Integrator> u_,
    RCP<Integrator> p_,
    RCP<Kinematics<T>> k_,
    RCP<Model<T>> m_,
    Disc* d_,
    Input* in_) {
  u = rcp_static_cast<Displacement<T>>(u_);
  p = rcp_static_cast<Pressure<T>>(p_);
  kin = k_;
  model = m_;
  disc = d_;
  dim = disc->dim;
  input = in_;
  this->name = "presidual";
}

template <typename T>
void PResidual<T>::set_elem_set(std::string const& set) {
  assert(input->mats.count(set));
  auto& mat = input->mats[set];
  double E = mat.E;
  double nu = mat.nu;
  mu = E / (2.0 * (1.0 + nu));
  lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  kappa = E / ( 3.0 * (1.0 - 2.0*nu));
}

template <typename T>
void PResidual<T>::do_small_strain(double ipw, double dv) {
  T pbar = 0.0;
  for (int i = 0; i < dim; ++i)
    pbar += u->grad(i, i);
  pbar *= (2.0*mu/dim + lambda);
  for (int n = 0; n < disc->num_p_elem_nodes; ++n)
    p->resid(n) += (p->val() - pbar) * p->BF[n] * ipw * dv;
}

template <typename T>
void PResidual<T>::do_large_strain(double ipw, double dv) {
  auto J = kin->J;
  T dUdJ = 0.5*(J - 1.0/J);
  for (int n = 0; n < disc->num_p_elem_nodes; ++n)
    p->resid(n) += (dUdJ - (p->val() / kappa)) * p->BF[n] * ipw * dv;
}

template <typename T>
void PResidual<T>::at_point(Vector const&, double ipw, double dv) {
  if (model->small_strain()) do_small_strain(ipw, dv);
  else do_large_strain(ipw, dv);
}

template struct PResidual<ST>;
template struct PResidual<FADT>;

}
