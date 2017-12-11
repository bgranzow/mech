#include "disc.hpp"
#include "qoi.hpp"
#include "linalg.hpp"
#include <PCU.h>

namespace mech {

QoI<ST>::QoI(Disc* d) {
  disc = d;
  elem = 0;
  qoi_value = 0.0;
  elem_value = 0.0;
}

QoI<ST>::~QoI() {
}

void QoI<ST>::pre_process(LinAlg*) {
  qoi_value = 0.0;
}

void QoI<ST>::gather(apf::MeshElement* me) {
  elem = me;
  elem_value = 0.0;
}

void QoI<ST>::scatter(LinAlg*) {
  qoi_value += elem_value;
  elem = 0;
}

void QoI<ST>::post_process(LinAlg*) {
  disc = 0;
  PCU_Add_Doubles(&qoi_value, 1);
}

QoI<FADT>::QoI(Disc* d) {
  disc = d;
  elem = 0;
  qoi_value = 0.0;
  elem_value = 0.0;
}

QoI<FADT>::~QoI() {
}

void QoI<FADT>::pre_process(LinAlg*) {
  qoi_value = 0.0;
}

void QoI<FADT>::gather(apf::MeshElement* me) {
  elem = me;
  auto num_dofs = disc->num_elem_dofs;
  elem_value.diff(0, num_dofs);
  elem_value.fastAccessDx(0) = 0.0;
}

void QoI<FADT>::scatter(LinAlg* la) {
  GIDs rows;
  auto ent = apf::getMeshEntity(elem);
  auto num_dofs = disc->num_elem_dofs;
  get_gids(disc, ent, rows);
  for (int dof = 0; dof < num_dofs; ++dof) {
    GID row = rows[dof];
    auto val = get_elem_value().fastAccessDx(dof);
    add_to_functional(la, row, val);
  }
  qoi_value += elem_value.val();
  elem = 0;
}

void QoI<FADT>::post_process(LinAlg*) {
  disc = 0;
  PCU_Add_Doubles(&qoi_value, 1);
}

}
