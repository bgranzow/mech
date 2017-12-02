#include "disc.hpp"
#include "residual.hpp"

namespace mech {

template <typename T>
Residual<T>::Residual(Disc* , Input*, int, bool) {
}

template <typename T>
void Residual<T>::set_elem_set(std::string const&) {
}

template <typename T>
void Residual<T>::gather(apf::MeshElement*) {
}

template <typename T>
void Residual<T>::at_point(Vector const&, double, double) {
}

template <typename T>
void Residual<T>::scatter(LinAlg*) {
}

template struct Residual<ST>;
template struct Residual<FADT>;

}
