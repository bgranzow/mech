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
  in->model = ELASTIC;
  in->geom_file = argv[1];
  in->mesh_file = argv[2];
  in->assoc_file = argv[3];
  Material mat1 = { 1000, .25, 100, 10 };
  in->mats["body"] = mat1;
  DBC dbc1 = { 0, zero, "xmin" };
  DBC dbc2 = { 0, zero, "ymin" };
  DBC dbc3 = { 0, tenth, "xmax" };
  in->dbcs = { dbc1, dbc2, dbc3 };
}

static void run(char** argv) {
  Input in;
  Disc disc;
  setup(&in, argv);
  init_disc(&disc, &in);
  build_disc_data(&disc);
  solve_linear_primal(&in, &disc, 0.0);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
