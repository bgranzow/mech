#include <control.hpp>
#include <disc.hpp>

namespace mech {

static void test_disc(Disc* d) {
  print("dim: %d", d->dim);
  print("elem type: %d", d->elem_type);
  print("num u elem nodes: %d", d->num_u_elem_nodes);
  print("num p elem nodes: %d", d->num_p_elem_nodes);
  print("num elem dofs: %d", d->num_elem_dofs);
  print("num owned dofs: %lu", d->num_owned_dofs);
  print("num total dofs: %lu", d->num_total_dofs);
  print("p dof offset: %lu", d->p_dof_offset);
}

static void run(char** argv) {
  Input in;
  Disc disc;
  in.geom_file = argv[1];
  in.mesh_file = argv[2];
  in.assoc_file = argv[3];
  init_disc(&disc, &in);
  build_disc_data(&disc);
  test_disc(&disc);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
