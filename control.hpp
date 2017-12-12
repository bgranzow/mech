#ifndef control_hpp
#define control_hpp

#include <Sacado_Fad_SLFad.hpp>
#include <Teuchos_RCP.hpp>

namespace apf {
class Vector3;
class VectorElement;
class Matrix3x3;
using MeshElement = VectorElement;
}

namespace mech {

struct Disc;
struct FineDisc;
struct LinAlg;
struct Integrator;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_static_cast;

using ST = double;
using FADT = Sacado::Fad::SLFad<ST, 70>;

using GID = long;
using GIDs = std::vector<GID>;

using Vector = apf::Vector3;
using Tensor = apf::Matrix3x3;

using Evaluators = std::vector<RCP<Integrator>>;

typedef double (*function)(Vector const& x, double t);

enum ModelType { ELASTIC, NEOHOOKEAN, PLASTIC };
enum EvalType { PRIMAL, ADJOINT, NONE };

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
  int model;
  std::string geom_file;
  std::string mesh_file;
  std::string assoc_file;
  std::vector<DBC> dbcs;
  std::map<std::string, Material> mats;
};

struct Integrator {
  virtual ~Integrator();
  virtual void set_time(double) {}
  virtual void pre_process(LinAlg*) {}
  virtual void set_elem_set(std::string const&) {}
  virtual void gather(apf::MeshElement*) {}
  virtual void in_elem(apf::MeshElement*) {}
  virtual void at_point(Vector const&, double, double) {}
  virtual void out_elem() {}
  virtual void scatter(LinAlg*) {}
  virtual void post_process(LinAlg*) {}
  std::string name;
};

void init(int* argc, char*** argv);
void free();
void print(const char* msg, ...);
void fail(const char* why, ...) __attribute__((noreturn));
double time();

}

#endif
