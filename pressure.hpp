#ifndef pressure_hpp
#define pressure_hpp

#include "control.hpp"
#include "disc.hpp"
#include <functional>

namespace mech {

template <typename T> struct Pressure;

template <>
struct Pressure<ST> : public Integrator {
  Pressure(Disc* d, int m);
  ST& val();
  ST& grad(int i);
  ST& nodal(int n);
  ST& resid(int n);
  void gather(apf::MeshElement* me);
  void at_point(Vector const& p, double, double);
  void scatter_none(LinAlg*);
  void scatter_primal(LinAlg* la);
  void scatter(LinAlg* la);
  std::function<void(Pressure<ST>*, LinAlg*)> op;
  int dim;
  Disc* disc;
  apf::MeshElement* elem;
  apf::NewArray<ST> BF;
  apf::NewArray<Vector> GBF;
  apf::NewArray<ST> node_st;
  ST value;
  Vector gradient;
  std::vector<ST> residual;
};

template <>
struct Pressure<FADT> : public Integrator {
  Pressure(Disc* d, int m);
  FADT& val();
  FADT& grad(int i);
  FADT& nodal(int n);
  FADT& resid(int n);
  void gather(apf::MeshElement* me);
  void at_point(Vector const& p, double, double);
  void scatter_none(LinAlg*);
  void scatter_primal(LinAlg* la);
  void scatter_adjoint(LinAlg* la);
  void scatter(LinAlg* la);
  std::function<void(Pressure<FADT>*, LinAlg*)> op;
  int dim;
  Disc* disc;
  apf::MeshElement* elem;
  apf::NewArray<ST> BF;
  apf::NewArray<Vector> GBF;
  apf::NewArray<ST> node_st;
  FADT value;
  std::vector<FADT> node_fadt;
  std::vector<FADT> gradient;
  std::vector<FADT> residual;
};

}

#endif
