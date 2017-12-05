#ifndef displacement_hpp
#define displacement_hpp

#include "control.hpp"
#include "disc.hpp"
#include <functional>

namespace mech {

template <typename T> struct Displacement;

template <>
struct Displacement<ST> : public Integrator {
  Displacement(Disc* d, int m);
  ST& val(int i);
  ST& grad(int i, int j);
  ST& nodal(int n, int i);
  ST& resid(int n, int i);
  void gather(apf::MeshElement* me);
  void at_point(Vector const& p, double, double);
  void scatter_none(LinAlg* la);
  void scatter_primal(LinAlg* la);
  void scatter(LinAlg* la);
  std::function<void(Displacement<ST>*, LinAlg*)> op;
  int dim;
  Disc* disc;
  apf::MeshElement* elem;
  apf::NewArray<ST> BF;
  apf::NewArray<Vector> GBF;
  apf::NewArray<Vector> node_st;
  Vector value;
  Tensor gradient;
  std::vector<ST> residual;
};

template <>
struct Displacement<FADT> : public Integrator {
  Displacement(Disc* d, int m);
  FADT& val(int i);
  FADT& grad(int i, int j);
  FADT& nodal(int n, int i);
  FADT& resid(int n, int i);
  void gather(apf::MeshElement* me);
  void at_point(Vector const& p, double, double);
  void scatter_none(LinAlg* la);
  void scatter_primal(LinAlg* la);
  void scatter_adjoint(LinAlg* la);
  void scatter(LinAlg* la);
  std::function<void(Displacement<FADT>*, LinAlg*)> op;
  int dim;
  Disc* disc;
  apf::MeshElement* elem;
  apf::NewArray<ST> BF;
  apf::NewArray<Vector> GBF;
  apf::NewArray<Vector> node_st;
  std::vector<FADT> node_fadt;
  std::vector<FADT> value;
  std::vector<FADT> gradient;
  std::vector<FADT> residual;
};

}

#endif
