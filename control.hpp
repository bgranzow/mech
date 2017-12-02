#ifndef control_hpp
#define control_hpp

#include <Sacado_Fad_SLFad.hpp>

namespace apf {
class Vector3;
class VectorElement;
using MeshElement = VectorElement;
}

namespace mech {

using ST = double;
using FADT = Sacado::Fad::SLFad<ST, 70>;

using GID = long;
using GIDs = std::vector<GID>;

using Vector = apf::Vector3;

typedef double (*function)(Vector const& x, double t);

enum Eval { PRIMAL, ADJOINT, NONE };

struct DBC {
  int eq;
  function val;
  std::string set;
};

struct Material {
  double E;
  double nu;
  double K;
  double Y;
};

struct Input {
  std::string geom_file;
  std::string mesh_file;
  std::string assoc_file;
  std::vector<DBC> dbcs;
  std::map<std::string, Material> mats;
};

struct LinAlg;

struct Integrator {
  virtual void set_time(double) {}
  virtual void pre_process(LinAlg*) {}
  virtual void set_elem_set(std::string const&) {}
  virtual void gather(apf::MeshElement*) {}
  virtual void at_point(Vector const&, double, double) {}
  virtual void scatter(LinAlg*) {}
  virtual void post_process(LinAlg*) {}
};

void init(int* argc, char*** argv);
void free();
void print(const char* msg, ...);
void fail(const char* why, ...) __attribute__((noreturn));
double time();

}

#endif
