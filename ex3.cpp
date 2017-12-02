#include <control.hpp>
#include <disc.hpp>
#include <pressure.hpp>

namespace mech {

template <typename T>
static void test_pressure(Disc* d) {
  apf::Vector3 xi(0,0,0);
  Pressure<T> u(d, NONE);
  for (auto& set : d->elem_sets) {
    auto es_name = set.first;
    auto elems = set.second;
    u.set_elem_set(es_name);
    for (size_t i = 0; i < elems.size(); ++i) {
      auto elem = elems[i];
      auto me = apf::createMeshElement(d->mesh, elem);
      u.gather(me);
      for (int n = 0; n < d->num_p_elem_nodes; ++n)
        std::cout << u.nodal(n) << std::endl;
      std::cout << std::endl;
      for (int ip = 0; ip < d->num_ips; ++ip) {
        apf::getIntPoint(me, d->q_degree, ip, xi);
        auto dv = apf::getDV(me, xi);
        auto w = apf::getIntWeight(me, d->q_degree, ip);
        u.at_point(xi, w, dv);
        std::cout << u.val() << std::endl;
        std::cout << std::endl;
        for (int i = 0; i < d->dim; ++i)
          std::cout << u.grad(i) << std::endl;
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
  test_pressure<FADT>(&disc);
  free_disc_data(&disc);
  free_disc(&disc);
}

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
