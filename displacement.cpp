#include "displacement.hpp"

namespace mech {

Displacement<ST>::Displacement(Disc* d, int) {
  disc = d;
  dim = disc->dim;
  elem = 0;
  residual.resize(dim * disc->num_u_elem_nodes);
  this->name = "u";
}

ST& Displacement<ST>::val(int i) {
  return value[i];
}

ST& Displacement<ST>::grad(int i, int j) {
  return gradient[i][j];
}

ST& Displacement<ST>::nodal(int n, int i) {
  return node_st[n][i];
}

ST& Displacement<ST>::resid(int n, int i) {
  return residual[n * dim + i];
}

void Displacement<ST>::gather(apf::MeshElement* me) {
  elem = me;
  auto u_elem = apf::createElement(disc->u, elem);
  apf::getVectorNodes(u_elem, node_st);
  apf::destroyElement(u_elem);
  for (size_t i = 0; i < residual.size(); ++i)
    residual[i] = 0;
}

void Displacement<ST>::at_point(apf::Vector3 const& p, double, double) {
  apf::getBF(disc->u_basis, elem, p, BF);
  apf::getGradBF(disc->u_basis, elem, p, GBF);
  for (int i = 0; i < dim; ++i) {
    value[i] = nodal(0, i) * BF[0];
    for (int n = 1; n < disc->num_u_elem_nodes; ++n)
      value[i] += nodal(n, i) * BF[n];
    for (int j = 0; j < dim; ++j) {
      gradient[i][j] = nodal(0, i) * GBF[0][j];
      for (int n = 1; n < disc->num_u_elem_nodes; ++n)
        gradient[i][j] += nodal(n, i) * GBF[n][j];
    }
  }
}

void Displacement<ST>::scatter(LinAlg*) {
  elem = 0;
}

Displacement<FADT>::Displacement(Disc* d, int) {
  disc = d;
  dim = disc->dim;
  elem = 0;
  value.resize(dim);
  gradient.resize(dim * dim);
  node_fadt.resize(dim * disc->num_u_elem_nodes);
  residual.resize(dim * disc->num_u_elem_nodes);
  this->name = "u";
}

FADT& Displacement<FADT>::val(int i) {
  return value[i];
}

FADT& Displacement<FADT>::grad(int i, int j) {
  return gradient[i * dim + j];
}

FADT& Displacement<FADT>::nodal(int n, int i) {
  return node_fadt[n * dim + i];
}

FADT& Displacement<FADT>::resid(int n, int i) {
  return residual[n * dim + i];
}

void Displacement<FADT>::gather(apf::MeshElement* me) {
  elem = me;
  auto u_elem = apf::createElement(disc->u, elem);
  apf::getVectorNodes(u_elem, node_st);
  apf::destroyElement(u_elem);
  for (int n = 0; n < disc->num_u_elem_nodes; ++n) {
    for (int d = 0; d < dim; ++d) {
      int eq = n*dim + d;
      nodal(n, d).diff(eq, disc->num_elem_dofs);
      nodal(n, d).val() = node_st[n][d];
      resid(n, d) = 0.0;
    }
  }
}

void Displacement<FADT>::at_point(apf::Vector3 const& p, double, double) {
  apf::getBF(disc->u_basis, elem, p, BF);
  apf::getGradBF(disc->u_basis, elem, p, GBF);
  for (int i = 0; i < dim; ++i) {
    val(i) = nodal(0, i) * BF[0];
    for (int n = 1; n < disc->num_u_elem_nodes; ++n)
      val(i) += nodal(n, i) * BF[n];
    for (int j = 0; j < dim; ++j) {
      grad(i, j) = nodal(0, i) * GBF[0][j];
      for (int n = 1; n < disc->num_u_elem_nodes; ++n)
        grad(i, j) += nodal(n, i) * GBF[n][j];
    }
  }
}

void Displacement<FADT>::scatter(LinAlg*) {
  elem = 0;
}

}
