#include <control.hpp>
#include <disc.hpp>
#include <mechanics.hpp>
#include <states.hpp>

namespace mech {

double zero(Vector const&, double) { return 0; }
double hundredth(Vector const&, double t) { return 0.01 * t; }

static void setup(Input* in, char** argv) {
  in->model = PLASTIC;
  in->geom_file = argv[1];
  in->mesh_file = argv[2];
  in->assoc_file = argv[3];
  Material mat1 = { 1000, .25, 100, 10 };
  in->mats["body"] = mat1;
  DBC dbc1 = { 0, zero, "xmin" };
  DBC dbc2 = { 1, zero, "ymin" };
  DBC dbc3 = { 0, hundredth, "xmax" };
  in->dbcs = { dbc1, dbc2, dbc3 };
}

static void output(Disc* d, std::string const& n) {
  auto p = apf::createFieldOn(d->mesh, "p2", apf::SCALAR);
  apf::projectField(p, d->p);
  apf::writeVtkFiles(n.c_str(), d->mesh);
  apf::destroyField(p);
}

static void run(char** argv) {
  Input in;
  Disc disc;
  setup(&in, argv);
  init_disc(&disc, &in);
  build_disc_data(&disc);
  double t = 0;
  double dt = 1.0;
  for (int i = 0; i < 5; ++i) {
    t += dt;
    auto name = "out_ex5_" + std::to_string(i);
    solve_nonlinear_primal(&in, &disc, t, 10, 1.0e-8);
    output(&disc, name);
    update_states(&disc);
  }
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
