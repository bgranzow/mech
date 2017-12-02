#ifndef residual_hpp
#define residual_hpp

#include "control.hpp"
#include <MiniTensor.h>

namespace mech {

struct Disc;

template <typename T>
struct Residual : public Integrator {

  Residual(Disc* d, Input* in, int mode, bool save);

  void set_elem_set(std::string const& set);
  void gather(apf::MeshElement* e);
  void at_point(Vector const& xi, double, double);
  void scatter(LinAlg* la);

  apf::NewArray<ST> uBF;
  apf::NewArray<Vector> uGBF;
  apf::NewArray<Vector> u_node_st;
  std::vector<minitensor::Vector<T>> u_node;
  minitensor::Vector<T> u;
  minitensor::Tensor<T> grad_u;

  apf::NewArray<ST> pBF;
  apf::NewArray<Vector> pGBF;
  apf::NewArray<ST> p_node_st;
  std::vector<T> p_node;
  T p;

  Disc* disc;
  Input* input;

  int mode;
  bool save_states;

  double E;
  double nu;
  double K;
  double Y;

};

}

#endif
