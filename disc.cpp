#include "disc.hpp"
#include <apfMDS.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <PCU.h>

namespace mech {

static void zero_all(Disc* d) {
  d->dim = 0;
  d->elem_type = 0;
  d->num_owned_dofs = 0;
  d->num_u_elem_nodes = 0;
  d->num_p_elem_nodes = 0;
  d->num_elem_dofs = 0;
  d->num_owned_dofs = 0;
  d->num_total_dofs = 0;
  d->p_dof_offset = 0;
  d->mesh = 0;
  d->sets = 0;
  d->u_basis = 0;
  d->p_basis = 0;
  d->u = 0;
  d->p = 0;
  d->first_pk = 0;
  d->eqps = 0;
  d->eqps_old = 0;
  d->Fp = 0;
  d->Fp_old = 0;
  d->u_nmbr = 0;
  d->p_nmbr = 0;
}

static void load_mesh(Disc* d, Input* in) {
  gmi_register_mesh();
  d->mesh = apf::loadMdsMesh(in->geom_file.c_str(), in->mesh_file.c_str());
  d->sets = apf::create_sets(d->mesh, in->assoc_file.c_str());
  d->dim = d->mesh->getDimension();
  d->elem_type = apf::getFirstType(d->mesh, d->dim);
}

static void init_dofs(Disc* d) {
  d->u_basis = apf::getSerendipity();
  d->p_basis = apf::getLagrange(1);
  d->u = apf::createField(d->mesh, "u2", apf::VECTOR, d->u_basis);
  d->p = apf::createField(d->mesh, "p1", apf::SCALAR, d->p_basis);
  apf::zeroField(d->u);
  apf::zeroField(d->p);
}

static void init_elem_info(Disc* d) {
  auto u_es = d->u_basis->getEntityShape(d->elem_type);
  auto p_es = d->p_basis->getEntityShape(d->elem_type);
  d->num_u_elem_nodes = u_es->countNodes();
  d->num_p_elem_nodes = p_es->countNodes();
  d->num_elem_dofs = d->num_u_elem_nodes*d->dim + d->num_p_elem_nodes;
  d->q_degree = 2;
  d->num_ips = (d->dim == 2) ? 3 : 4;
}

static void identitize(Disc* d) {
  apf::MeshEntity* elem;
  apf::Matrix3x3 I2(1,0,0,0,1,0,0,0,0);
  apf::Matrix3x3 I3(1,0,0,0,1,0,0,0,0);
  auto I = (d->dim == 2) ? I2 : I3;
  auto it = d->mesh->begin(d->dim);
  while ((elem = d->mesh->iterate(it))) {
    for (int ip = 0; ip < d->num_ips; ++ip) {
      apf::setMatrix(d->Fp, elem, ip, I);
      apf::setMatrix(d->Fp_old, elem, ip, I);
    }
  }
  d->mesh->end(it);
}

static void init_states(Disc* d) {
  auto basis = apf::getIPFitShape(d->dim, 2);
  d->first_pk = apf::createField(d->mesh, "first_pk", apf::MATRIX, basis);
  d->eqps = apf::createField(d->mesh, "eqps", apf::SCALAR, basis);
  d->eqps_old = apf::createField(d->mesh, "eqps_old", apf::SCALAR, basis);
  d->Fp = apf::createField(d->mesh, "Fp", apf::MATRIX, basis);
  d->Fp_old = apf::createField(d->mesh, "Fp_old", apf::MATRIX, basis);
  apf::zeroField(d->first_pk);
  apf::zeroField(d->eqps);
  apf::zeroField(d->eqps_old);
  identitize(d);
}

void init_disc(Disc* d, Input* in) {
  zero_all(d);
  load_mesh(d, in);
  init_dofs(d);
  init_elem_info(d);
  init_states(d);
}

void free_disc(Disc* d) {
  delete d->sets;
  apf::destroyField(d->u);
  apf::destroyField(d->p);
  d->mesh->destroyNative();
  apf::destroyMesh(d->mesh);
  zero_all(d);
}

static void init_numbers(Disc* d) {
  auto un = apf::numberOwnedNodes(d->mesh, "un", d->u_basis);
  auto pn = apf::numberOwnedNodes(d->mesh, "pn", d->p_basis);
  auto num_owned_u_nodes = apf::countNodes(un);
  auto num_owned_p_nodes = apf::countNodes(pn);
  auto num_total_u_nodes = PCU_Add_Long(num_owned_u_nodes);
  auto num_total_p_nodes = PCU_Add_Long(num_owned_p_nodes);
  d->num_owned_dofs = num_owned_u_nodes*d->dim + num_owned_p_nodes;
  d->num_total_dofs = num_total_u_nodes*d->dim + num_total_p_nodes;
  d->p_dof_offset = num_total_u_nodes*d->dim;
  d->u_nmbr = apf::makeGlobal(un);
  d->p_nmbr = apf::makeGlobal(pn);
  apf::synchronize(d->u_nmbr);
  apf::synchronize(d->p_nmbr);
}

static void init_sets(Disc* d) {
  d->elem_sets = apf::get_elem_sets(d->mesh, d->sets);
  d->side_sets = apf::get_side_sets(d->mesh, d->sets);
  d->u_node_sets = apf::get_node_sets(d->mesh, d->sets, d->u_nmbr);
  d->p_node_sets = apf::get_node_sets(d->mesh, d->sets, d->p_nmbr);
}

void build_disc_data(Disc* d) {
  auto t0 = time();
  init_numbers(d);
  init_sets(d);
  auto t1 = time();
  print("disc data built in %f seconds", t1 - t0);
}

static void free_sets(Disc* d) {
  d->elem_sets.clear();
  d->side_sets.clear();
  d->u_node_sets.clear();
  d->p_node_sets.clear();
}

static void free_numbers(Disc* d) {
  apf::destroyGlobalNumbering(d->u_nmbr);
  apf::destroyGlobalNumbering(d->p_nmbr);
  d->u_nmbr = 0;
  d->p_nmbr = 0;
  d->num_owned_dofs = 0;
  d->num_total_dofs = 0;
  d->p_dof_offset = 0;
}

void free_disc_data(Disc* d) {
  free_sets(d);
  free_numbers(d);
}

GID get_u_gid(Disc* d, apf::MeshEntity* e, int n, int i) {
  apf::NewArray<GID> node_ids;
  apf::getElementNumbers(d->u_nmbr, e, node_ids);
  return node_ids[n] * d->dim + i;
}

GID get_p_gid(Disc* d, apf::MeshEntity* e, int n) {
  apf::NewArray<GID> node_ids;
  apf::getElementNumbers(d->p_nmbr, e, node_ids);
  return d->p_dof_offset + node_ids[n];
}

GID get_u_gid(Disc* d, apf::Node const& n, int i) {
  GID nmbr = apf::getNumber(d->u_nmbr, n);
  return nmbr * d->dim + i;
}

void get_gids(Disc* d, apf::MeshEntity* e, GIDs& ids) {
  apf::NewArray<long> u_node_ids;
  apf::NewArray<long> p_node_ids;
  apf::getElementNumbers(d->u_nmbr, e, u_node_ids);
  apf::getElementNumbers(d->p_nmbr, e, p_node_ids);
  ids.resize(d->num_elem_dofs);
  int eq = 0;
  for (int n = 0; n < d->num_u_elem_nodes; ++n)
  for (int i = 0; i < d->dim; ++i)
    ids[eq++] = u_node_ids[n] * d->dim + i;
  for (int n = 0; n < d->num_p_elem_nodes; ++n)
    ids[eq++] = d->p_dof_offset + p_node_ids[n];
}

}
