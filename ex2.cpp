#include <control.hpp>
#include <disc.hpp>
#include <residual.hpp>

namespace mech {

static double zero(Vector const&, double) { return 0; }
static double tenth(Vector const&, double) { return 0.1; }

template <typename T>
static void test_residual(Disc* d, Input* in) {
  Residual<T> resid(d, in, NONE, true);
}

static void set_input(char** argv, Input* in) {
  in->geom_file = argv[1];
  in->mesh_file = argv[2];
  in->assoc_file = argv[3];
  DBC dbc1 = { 0, zero, "xmin" };
  DBC dbc2 = { 0, zero, "ymin" };
  DBC dbc3 = { 0, tenth, "xmax" };
  Material mat1 = { 1000, .25, 100, 10 };
  in->dbcs = { dbc1, dbc2, dbc3 };
  in->mats["body"] = mat1;
}

static void run(char** argv) {
  Input in;
  Disc disc;
  set_input(argv, &in);
  init_disc(&disc, &in);
  build_disc_data(&disc);
  test_residual<FADT>(&disc, &in);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
