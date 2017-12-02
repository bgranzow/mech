#include "model.hpp"

namespace mech {

template <typename T>
Model<T>::Model() {
  this->name = "model";
}

template <typename T>
Model<T>::~Model() {
}

template struct Model<ST>;
template struct Model<FADT>;

}
