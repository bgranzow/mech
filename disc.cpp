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
  apf::reorderMdsMesh(d->mesh);
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

static void number_owned(apf::GlobalNumbering* nmbr,
    apf::FieldShape* basis, int comps, long& i) {
  apf::MeshEntity* ent;
  auto m = apf::getMesh(nmbr);
  for (int dim = 0; dim < 4; ++dim) {
    if (! basis->hasNodesIn(dim)) continue;
    auto it = m->begin(dim);
    while ((ent = m->iterate(it))) {
      if (! m->isOwned(ent)) continue;
      auto type = m->getType(ent);
      int nnodes = basis->countNodesOn(type);
      for (int n = 0; n < nnodes; ++n) {
        apf::number(nmbr, ent, n, i);
        i += comps;
      }
    }
    m->end(it);
  }
}

static void offset_nmbr(apf::GlobalNumbering* nmbr, long offset) {
  apf::DynamicArray<apf::Node> nodes;
  apf::getNodes(nmbr, nodes);
  for (size_t n = 0; n < nodes.size(); ++n) {
    auto node = nodes[n];
    auto ent = node.entity;
    auto lnode = node.node;
    auto old = apf::getNumber(nmbr, ent, lnode);
    apf::number(nmbr, ent, lnode, old + offset);
  }
}

static void init_numbers(Disc* d) {
  long i = 0;
  d->u_nmbr = apf::createGlobalNumbering(d->mesh, "un", d->u_basis);
  d->p_nmbr = apf::createGlobalNumbering(d->mesh, "pn", d->p_basis);
  number_owned(d->u_nmbr, d->u_basis, d->dim, i);
  number_owned(d->p_nmbr, d->p_basis, 1, i);
  auto num_owned_u_dofs = apf::countNodes(d->u_nmbr) * d->dim;
  auto num_owned_p_dofs = apf::countNodes(d->p_nmbr);
  d->num_owned_dofs = num_owned_u_dofs + num_owned_p_dofs;
  d->num_total_dofs = PCU_Add_Long(d->num_owned_dofs);
  auto offset = PCU_Exscan_Long(d->num_owned_dofs);
  offset_nmbr(d->u_nmbr, offset);
  offset_nmbr(d->p_nmbr, offset);
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
}

void free_disc_data(Disc* d) {
  free_sets(d);
  free_numbers(d);
}

GID get_u_gid(Disc* d, apf::MeshEntity* e, int n, int i) {
  apf::NewArray<GID> node_ids;
  apf::getElementNumbers(d->u_nmbr, e, node_ids);
  return node_ids[n] + i;
}

GID get_p_gid(Disc* d, apf::MeshEntity* e, int n) {
  apf::NewArray<GID> node_ids;
  apf::getElementNumbers(d->p_nmbr, e, node_ids);
  return node_ids[n];
}

GID get_u_gid(Disc* d, apf::Node const& n, int i) {
  GID nmbr = apf::getNumber(d->u_nmbr, n);
  return nmbr + i;
}

GID get_p_gid(Disc* d, apf::Node const& n) {
  return apf::getNumber(d->p_nmbr, n);
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
    ids[eq++] = u_node_ids[n] + i;
  for (int n = 0; n < d->num_p_elem_nodes; ++n)
    ids[eq++] = p_node_ids[n];
}

}
