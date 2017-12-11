#include "control.hpp"
#include "disc.hpp"
#include "states.hpp"
#include <apf.h>

namespace mech {

static ST get_val(ST v) {
  return v;
}

static ST get_val(FADT const& v) {
  return v.val();
}

static void zero(Tensor& v) {
  for (int i = 0; i < 3; ++i)
  for (int j = 0; j < 3; ++j)
    v[i][j] = 0.0;
}

template <typename T>
void get_scalar(apf::Element* e, Vector const& p, T& v) {
  auto val = apf::getScalar(e, p);
  v = (T)val;
}

template <typename T>
void get_tensor(apf::Element* e, Vector const& p, TensorT<T>& v) {
  Tensor val;
  apf::getMatrix(e, p, val);
  for (size_t i = 0; i < v.get_dimension(); ++i)
  for (size_t j = 0; j < v.get_dimension(); ++j)
    v(i,j) = (T)val[i][j];
}

template <typename T>
void set_scalar(apf::Field* f, apf::MeshEntity* e, int ip, T const& v) {
  apf::setScalar(f, e, ip, get_val(v));
}

template <typename T>
void set_tensor(apf::Field* f, apf::MeshEntity* e, int ip, TensorT<T> const& v) {
  Tensor val;
  zero(val);
  for (size_t i = 0; i < v.get_dimension(); ++i)
  for (size_t j = 0; j < v.get_dimension(); ++j)
    val[i][j] = get_val(v(i,j));
  apf::setMatrix(f, e, ip, val);
}

void update_states(Disc* d) {
  apf::copyData(d->Fp_old, d->Fp);
  apf::copyData(d->eqps_old, d->eqps);
}

template void get_scalar(apf::Element* e, Vector const& p, ST& v);
template void get_scalar(apf::Element* e, Vector const& p, FADT& v);

template void get_tensor(apf::Element* e, Vector const& p, TensorT<ST>& v);
template void get_tensor(apf::Element* e, Vector const& p, TensorT<FADT>& v);

template void set_scalar(apf::Field* f, apf::MeshEntity* e, int ip, ST const& v);
template void set_scalar(apf::Field* f, apf::MeshEntity* e, int ip, FADT const& v);

template void set_tensor(apf::Field* f, apf::MeshEntity* e, int ip, TensorT<ST> const& v);
template void set_tensor(apf::Field* f, apf::MeshEntity* e, int ip, TensorT<FADT> const& v);

}
