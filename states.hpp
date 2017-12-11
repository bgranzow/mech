#ifndef states_hpp
#define states_hpp

#include <MiniTensor.h>

namespace apf {
class Element;
class Field;
class MeshEntity;
}

namespace mech {

template <typename T>
using TensorT = minitensor::Tensor<T>;

template <typename T>
void get_scalar(apf::Element* e, Vector const& p, T& v);

template <typename T>
void get_tensor(apf::Element* e, Vector const& p, TensorT<T>& v);

template <typename T>
void set_scalar(apf::Field* f, apf::MeshEntity* e, int ip, T const& v);

template <typename T>
void set_tensor(apf::Field* f, apf::MeshEntity* e, int ip, TensorT<T> const& v);

void update_states(Disc* d);

}

#endif
