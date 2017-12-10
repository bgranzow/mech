#include <control.hpp>
#include <disc.hpp>
#include <displacement.hpp>
#include <mechanics.hpp>
#include <pressure.hpp>
#include <linalg.hpp>

namespace mech {

double zero(Vector const&, double) { return 0; }
double tenth(Vector const&, double) { return 0.1; }

static void setup(Input* in, char** argv) {
  in->model = NEOHOOKEAN;
  in->geom_file = argv[1];
  in->mesh_file = argv[2];
  in->assoc_file = argv[3];
  Material mat1 = { 1000, .25, 100, 10 };
  in->mats["body"] = mat1;
  DBC dbc1 = { 0, zero, "xmin" };
  DBC dbc2 = { 1, zero, "ymin" };
  DBC dbc3 = { 0, tenth, "xmax" };
  in->dbcs = { dbc1, dbc2, dbc3 };
}

static void output(Disc* d) {
  auto p = apf::createFieldOn(d->mesh, "p2", apf::SCALAR);
  apf::projectField(p, d->p);
  std::vector<std::string> names = { "u2", "p2", "first_pk" };
  apf::writeVtkFiles("out_ex4", d->mesh, names);
  apf::destroyField(p);
}

static void run(char** argv) {
  Input in;
  Disc disc;
  setup(&in, argv);
  init_disc(&disc, &in);
  build_disc_data(&disc);
  solve_nonlinear_primal(&in, &disc, 0.0, 5, 1.0e-8);
  output(&disc);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
