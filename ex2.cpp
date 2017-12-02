#include <control.hpp>
#include <disc.hpp>
#include <displacement.hpp>

namespace mech {

template <typename T>
static void test_displacement(Disc* d) {
  apf::Vector3 xi(0,0,0);
  Displacement<T> u(d, NONE);
  for (auto& set : d->elem_sets) {
    auto es_name = set.first;
    auto elems = set.second;
    u.set_elem_set(es_name);
    for (size_t i = 0; i < elems.size(); ++i) {
      auto elem = elems[i];
      auto me = apf::createMeshElement(d->mesh, elem);
      u.gather(me);
      for (int n = 0; n < d->num_u_elem_nodes; ++n)
      for (int i = 0; i < d->dim; ++i)
        std::cout << u.nodal(n, i) << std::endl;
      std::cout << std::endl;

      for (int ip = 0; ip < d->num_ips; ++ip) {
        apf::getIntPoint(me, d->q_degree, ip, xi);
        auto dv = apf::getDV(me, xi);
        auto w = apf::getIntWeight(me, d->q_degree, ip);
        u.at_point(xi, w, dv);
        for (int i = 0; i < d->dim; ++i)
          std::cout << u.val(i) << std::endl;
        std::cout << std::endl;
        for (int i = 0; i < d->dim; ++i)
        for (int j = 0; j < d->dim; ++j)
          std::cout << u.grad(i,j) << std::endl;
        std::cout << std::endl;
      }
      apf::destroyMeshElement(me);
    }
  }
}

static void run(char** argv) {
  Input in;
  Disc disc;
  in.geom_file = argv[1];
  in.mesh_file = argv[2];
  in.assoc_file = argv[3];
  init_disc(&disc, &in);
  build_disc_data(&disc);
  test_displacement<FADT>(&disc);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
