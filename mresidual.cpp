#include "displacement.hpp"
#include "first_pk.hpp"
#include "mresidual.hpp"

namespace mech {

template <typename T>
MResidual<T>::MResidual(RCP<FirstPK<T>> P_, RCP<Integrator> u_, Disc* d) {
  first_pk = P_;
  u = rcp_static_cast<Displacement<T>>(u_);
  disc = d;
  dim = disc->dim;
  this->name = "mresidual";
}

template <typename T>
void MResidual<T>::at_point(Vector const&, double ipw, double dv) {
  auto P = first_pk->val();
  for (int n = 0; n < disc->num_u_elem_nodes; ++n)
  for (int i = 0; i < dim; ++i)
  for (int j = 0; j < dim; ++j)
    u->resid(n, i) += P(i, j) * u->GBF[n][j] * ipw * dv;
}

template struct MResidual<ST>;
template struct MResidual<FADT>;

}
