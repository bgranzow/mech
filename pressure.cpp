#include "pressure.hpp"

namespace mech {

Pressure<ST>::Pressure(Disc* d, int) {
  disc = d;
  dim = disc->dim;
  elem = 0;
  residual.resize(disc->num_p_elem_nodes);
}

ST& Pressure<ST>::val() {
  return value;
}

ST& Pressure<ST>::grad(int i) {
  return gradient[i];
}

ST& Pressure<ST>::nodal(int n) {
  return node_st[n];
}

ST& Pressure<ST>::resid(int n) {
  return residual[n];
}

void Pressure<ST>::gather(apf::MeshElement* me) {
  elem = me;
  auto p_elem = apf::createElement(disc->p, elem);
  apf::getScalarNodes(p_elem, node_st);
  apf::destroyElement(p_elem);
  for (size_t i = 0; i < residual.size(); ++i)
    residual[i] = 0;
}

void Pressure<ST>::at_point(apf::Vector3 const& p, double, double) {
  apf::getBF(disc->p_basis, elem, p, BF);
  apf::getGradBF(disc->p_basis, elem, p, GBF);
  value = nodal(0) * BF[0];
  for (int n = 1; n < disc->num_p_elem_nodes; ++n)
    value += nodal(n) * BF[n];
  for (int i = 0; i < dim; ++i) {
    gradient[i] = nodal(0) * GBF[0][i];
    for (int n = 1; n < disc->num_p_elem_nodes; ++n)
      gradient[i] += nodal(n) * GBF[n][i];
  }
}

void Pressure<ST>::scatter(LinAlg*) {
  elem = 0;
}

Pressure<FADT>::Pressure(Disc* d, int) {
  disc = d;
  dim = disc->dim;
  elem = 0;
  gradient.resize(dim);
  node_fadt.resize(disc->num_p_elem_nodes);
  residual.resize(disc->num_p_elem_nodes);
}

FADT& Pressure<FADT>::val() {
  return value;
}

FADT& Pressure<FADT>::grad(int i) {
  return gradient[i];
}

FADT& Pressure<FADT>::nodal(int n) {
  return node_fadt[n];
}

FADT& Pressure<FADT>::resid(int n) {
  return residual[n];
}

void Pressure<FADT>::gather(apf::MeshElement* me) {
  elem = me;
  auto p_elem = apf::createElement(disc->p, elem);
  apf::getScalarNodes(p_elem, node_st);
  apf::destroyElement(p_elem);
  for (int n = 0; n < disc->num_p_elem_nodes; ++n) {
    int eq = disc->num_u_elem_nodes*dim + n;
    nodal(n).diff(eq, disc->num_elem_dofs);
    nodal(n).val() = node_st[n];
    resid(n) = 0.0;
  }
}

void Pressure<FADT>::at_point(apf::Vector3 const& p, double, double) {
  apf::getBF(disc->p_basis, elem, p, BF);
  apf::getGradBF(disc->p_basis, elem, p, GBF);
  value = nodal(0) * BF[0];
  for (int n = 1; n < disc->num_p_elem_nodes; ++n)
    value += nodal(n) * BF[n];
  for (int i = 0; i < dim; ++i) {
    gradient[i] = nodal(0) * GBF[0][i];
    for (int n = 1; n < disc->num_p_elem_nodes; ++n)
      gradient[i] += nodal(n) * GBF[n][i];
  }
}

void Pressure<FADT>::scatter(LinAlg*) {
  elem = 0;
}

}
