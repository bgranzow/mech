#include "control.hpp"
#include "disc.hpp"
#include "linalg.hpp"
#include <apfShape.h>

namespace mech {

static void get_coord(Disc* d, apf::Node const& n, Vector& x) {
  Vector xi(0,0,0);
  auto me = apf::createMeshElement(d->mesh, n.entity);
  auto type = d->mesh->getType(n.entity);
  d->u_basis->getNodeXi(type, n.node, xi);
  apf::mapLocalToGlobal(me, xi, x);
  apf::destroyMeshElement(me);
}

void set_resid_dbcs(Input* in, Disc* d, LinAlg* la, double t) {
  Vector x(0,0,0);
  Vector disp(0,0,0);
  for (size_t i = 0; i < in->dbcs.size(); ++i) {
    auto& dbc = in->dbcs[i];
    auto set = dbc.set;
    auto idx = dbc.eq;
    assert(d->u_node_sets.count(set));
    auto nodes = d->u_node_sets[set];
    for (size_t node = 0; node < nodes.size(); ++node) {
      auto n = nodes[node];
      get_coord(d, n, x);
      apf::getVector(d->u, n.entity, n.node, disp);
      double sol = disp[idx];
      double v = dbc.val(x, t);
      GID row = get_u_gid(d, n, idx);
      set_to_residual(la, row, sol - v);
    }
  }
}

void set_jacob_dbcs(Input* in, Disc* d, LinAlg* la, double t) {
  Vector x(0,0,0);
  Vector disp(0,0,0);
  for (size_t i = 0; i < in->dbcs.size(); ++i) {
    auto& dbc = in->dbcs[i];
    auto set = dbc.set;
    auto idx = dbc.eq;
    assert(d->u_node_sets.count(set));
    auto nodes = d->u_node_sets[set];
    for (size_t node = 0; node < nodes.size(); ++node) {
      auto n = nodes[node];
      get_coord(d, n, x);
      apf::getVector(d->u, n.entity, n.node, disp);
      double sol = disp[idx];
      double v = dbc.val(x, t);
      GID row = get_u_gid(d, n, idx);
      set_to_residual(la, row, sol - v);
      diag_jacobian_row(la, row);
    }
  }
}

}
