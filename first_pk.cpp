#include "first_pk.hpp"
#include "kinematics.hpp"
#include "model.hpp"
#include "pressure.hpp"
#include "states.hpp"

namespace mech {

template <typename T>
FirstPK<T>::FirstPK(
    RCP<Model<T>> m,
    RCP<Integrator> k,
    RCP<Integrator> p,
    Disc* d,
    bool save) {
  model = m;
  kin = rcp_static_cast<Kinematics<T>>(k);
  pressure = rcp_static_cast<Pressure<T>>(p);
  disc = d;
  dim = disc->dim;
  save_states = save;
  first_pk.set_dimension(dim);
  I = minitensor::eye<T>(dim);
  this->name = "firstpk";
}

template <typename T>
void FirstPK<T>::in_elem(apf::MeshElement* me) {
  elem = apf::getMeshEntity(me);
  ip = 0;
}

template <typename T>
void FirstPK<T>::do_small_strain() {
  auto sigma_dev = model->get_sigma_dev();
  auto p = pressure->val();
  first_pk = sigma_dev + p * I;
}

template <typename T>
void FirstPK<T>::do_large_strain() {
  auto sigma_dev = model->get_sigma_dev();
  auto p = pressure->val();
  auto J = kin->J;
  auto F = kin->F;
  auto Finv = minitensor::inverse(F);
  first_pk = J * (sigma_dev + p*I) * minitensor::transpose(Finv);
}

template <typename T>
void FirstPK<T>::at_point(Vector const&, double, double) {
  if (model->small_strain()) do_small_strain();
  else do_large_strain();
  if (save_states)
    set_tensor<T>(disc->first_pk, elem, ip, first_pk);
  ip++;
}

template <typename T>
void FirstPK<T>::out_elem() {
  elem = 0;
}

template struct FirstPK<ST>;
template struct FirstPK<FADT>;

}
