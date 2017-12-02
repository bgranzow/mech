#include <assembly.hpp>
#include <control.hpp>
#include <disc.hpp>
#include <displacement.hpp>
#include <elastic.hpp>
#include <kinematics.hpp>
#include <pressure.hpp>

namespace mech {

template <typename T>
static void test_assembly(Disc* d, Input* in) {
  Evaluators E;
  Displacement<T> disp(d, NONE);
  Pressure<T> press(d, NONE);
  E.push_back(disp);
  E.push_back(press);
  auto u = find_evaluator(E, "u");
  Kinematics<T> kin(&u);
  Elastic<T> model(&u, d, in);
  E.push_back(kin);
  E.push_back(model);
}

static void run(char** argv) {
  Input in;
  Disc disc;
  in.geom_file = argv[1];
  in.mesh_file = argv[2];
  in.assoc_file = argv[3];
  Material mat1 = { 1000, .25, 100, 10 };
  in.mats["body"] = mat1;
  init_disc(&disc, &in);
  build_disc_data(&disc);
  test_assembly<ST>(&disc, &in);
  test_assembly<FADT>(&disc, &in);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
