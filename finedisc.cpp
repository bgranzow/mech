#include "finedisc.hpp"
#include <apfShape.h>

namespace mech {

static void zero_all(FineDisc* fd) {
  fd->dim = 0;
  fd->elem_type = 0;
  fd->num_owned_dofs = 0;
  fd->num_u_elem_nodes = 0;
  fd->num_p_elem_nodes = 0;
  fd->num_elem_dofs = 0;
  fd->num_owned_dofs = 0;
  fd->num_total_dofs = 0;
  fd->mesh = 0;
  fd->sets = 0;
  fd->u_basis = 0;
  fd->p_basis = 0;
  fd->u = 0;
  fd->p = 0;
  fd->zu = 0;
  fd->zp = 0;
  fd->first_pk = 0;
  fd->eqps = 0;
  fd->eqps_old = 0;
  fd->Fp = 0;
  fd->Fp_old = 0;
  fd->u_nmbr = 0;
  fd->p_nmbr = 0;
}

static void init_mesh(FineDisc* fd, Disc* d) {
  fd->mesh = d->mesh;
  fd->sets = d->sets;
  fd->dim = d->dim;
  fd->elem_type = d->elem_type;
}

static void init_dofs(FineDisc* fd, Disc* d) {
  fd->u_basis = apf::getLagrange(3);
  fd->p_basis = apf::getSerendipity();
  fd->u = apf::createField(fd->mesh, "u3", apf::VECTOR, fd->u_basis);
  fd->p = apf::createField(fd->mesh, "p2", apf::SCALAR, fd->p_basis);
  fd->zu = apf::createField(fd->mesh, "zu3", apf::VECTOR, fd->u_basis);
  fd->zp = apf::createField(fd->mesh, "zp2", apf::SCALAR, fd->p_basis);
  apf::projectField(fd->u, d->u);
  apf::projectField(fd->p, d->p);
  apf::zeroField(fd->zu);
  apf::zeroField(fd->zp);
}

static void init_elem_info(FineDisc* fd) {
  auto u_es = fd->u_basis->getEntityShape(fd->elem_type);
  auto p_es = fd->p_basis->getEntityShape(fd->elem_type);
  fd->num_u_elem_nodes = u_es->countNodes();
  fd->num_p_elem_nodes = p_es->countNodes();
  fd->num_elem_dofs = fd->num_u_elem_nodes*fd->dim + fd->num_p_elem_nodes;
  fd->q_degree = 4;
  fd->num_ips = (fd->dim == 2) ? 6 : 11;
}

static void init_states(FineDisc* fd, Disc* d) {
  fd->first_pk = d->first_pk;
  fd->eqps = d->eqps;
  fd->eqps_old = d->eqps_old;
  fd->Fp = d->Fp;
  fd->Fp_old = d->Fp_old;
}

void init_fine_disc(FineDisc* fd, Disc* d) {
  zero_all(fd);
  init_mesh(fd, d);
  init_dofs(fd, d);
  init_elem_info(fd);
  init_states(fd, d);
}

void free_fine_disc(FineDisc* fd) {
  apf::destroyField(fd->u);
  apf::destroyField(fd->p);
  apf::destroyField(fd->zu);
  apf::destroyField(fd->zp);
  zero_all(fd);
}


}
