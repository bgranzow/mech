#include "disc.hpp"
#include "displacement.hpp"
#include "J2.hpp"
#include "kinematics.hpp"
#include "states.hpp"

namespace mech {

template <typename T>
J2<T>::J2(RCP<Kinematics<T>> k, Disc* d, Input* in, bool save_) {
  kin = k;
  disc = d;
  input = in;
  dim = disc->dim;
  save = save_;
  Fp.set_dimension(dim);
  Fpinv.set_dimension(dim);
  Cpinv.set_dimension(dim);
  be.set_dimension(dim);
  s.set_dimension(dim);
  N.set_dimension(dim);
  Fpn.set_dimension(dim);
  sigma_dev.set_dimension(dim);
  sq23 = std::sqrt(2.0/3.0);
}

template <typename T>
void J2<T>::set_elem_set(std::string const& set) {
  assert(input->mats.count(set));
  auto& mat = input->mats[set];
  double E = mat.E;
  double nu = mat.nu;
  K = mat.K;
  Y = mat.Y;
  mu = E / (2.0 * (1.0 + nu));
  kappa = E / (3.0 * (1.0 - 2.0 * nu));
}

template <typename T>
void J2<T>::in_elem(apf::MeshElement* me) {
  ip = 0;
  elem = apf::getMeshEntity(me);
  Fpold_elem = apf::createElement(disc->Fp_old, me);
  eqpsold_elem = apf::createElement(disc->eqps_old, me);
}

template <typename T>
void J2<T>::at_point(Vector const& p, double, double) {

  // get def grad quantities
  auto F = kin->F;
  auto J = kin->J;
  Jm23 = std::pow(J, -2.0/3.0);

  // get plastic def grad quantities
  mech::get_tensor<T>(Fpold_elem, p, Fp);
  Fpinv = minitensor::inverse(Fp);

  // compute the trial state
  Cpinv = Fpinv * minitensor::transpose(Fpinv);
  be = Jm23 * F * Cpinv * minitensor::transpose(F);
  s = mu * minitensor::dev(be);
  mubar = minitensor::trace(be) * mu / dim;

  // check the yield condition
  smag = minitensor::norm(s);
  mech::get_scalar<T>(eqpsold_elem, p, eqps);
  f = smag - sq23 * (Y + K * eqps);

  // plastic increment
  if (f > 1.0e-12) {

    int iter = 0;
    bool converged = false;
    dgam = 0.0;

    T H = 0.0;
    T dH = 0.0;
    T alpha = 0.0;
    T res = 0.0;

    T X = 0.0;
    T R = f;
    T dRdX = -1.0 * mubar * (1.0 + H / (3.0 * mubar));

    while ((! converged) && (iter < 30)) {
      iter++;
      X = X - R / dRdX;
      alpha = eqps + sq23 * X;
      H = K * alpha;
      dH = K;
      R = smag - (2.0 * mubar * X + sq23 * (Y + H));
      dRdX = -2.0 * mubar * (1.0 + dH / (3.0 * mubar));
      res = std::abs(R);
      if ((res < 1.0e-11) || (res / Y < 1.0e-11) || (res / f < 1.0e-11))
        converged = true;
      if (iter == 30)
        fail("J2 : return mapping failed");
    }

    // updates
    dgam = X;
    N = (1.0 / smag) * s;
    s -= 2.0 * mubar * dgam * N;
    if (save) mech::set_scalar<T>(disc->eqps, elem, ip, alpha);

    // get Fpn
    Fpn = minitensor::exp(dgam * N) * Fp;
    if (save) mech::set_tensor<T>(disc->Fp, elem, ip, Fpn);
  }

  // otherwise elastic increment
  else
    if (save) mech::set_scalar(disc->eqps, elem, ip, eqps);

  sigma_dev = s/J;

  // increment integration point counter
  ip++;

}

template <typename T>
void J2<T>::out_elem() {
  elem = 0;
  apf::destroyElement(Fpold_elem);
  apf::destroyElement(eqpsold_elem);
}

template struct J2<ST>;
template struct J2<FADT>;

}
