#include "disc.hpp"
#include "linalg.hpp"

#include <PCU.h>

#define CALL(function) assert(0 == (function))

namespace mech {

LinAlg::LinAlg(Disc* d) {
  PetscInt n = d->num_owned_dofs;
  PetscInt N = d->num_total_dofs;
  CALL(VecCreateMPI(PETSC_COMM_WORLD, n, N, &f));
  CALL(VecCreateMPI(PETSC_COMM_WORLD, n, N, &dx));
  CALL(VecCreateMPI(PETSC_COMM_WORLD, n, N, &q));
  CALL(MatCreate(PETSC_COMM_WORLD, &J));
  CALL(MatSetType(J, MATMPIAIJ));
  CALL(MatSetSizes(J, n, n, N, N));
  CALL(MatSetUp(J));
  CALL(MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
}

LinAlg::~LinAlg() {
  CALL(MatDestroy(&J));
  CALL(VecDestroy(&f));
  CALL(VecDestroy(&dx));
  CALL(VecDestroy(&q));
}

void add_to_residual(LinAlg* la, GID row, double val) {
  PetscInt r[1] = { PetscInt(row) };
  PetscScalar v[1] = { val };
  CALL(VecSetValues(la->f, 1, r, v, ADD_VALUES));
}

void set_to_residual(LinAlg* la, GID row, double val) {
  PetscInt r[1] = { PetscInt(row) };
  PetscScalar v[1] { val };
  CALL(VecSetValues(la->f, 1, r, v, INSERT_VALUES));
}

void add_to_jacobian(LinAlg* la, GID row, GIDs const& cols, FADT const& val) {
  int sz = cols.size();
  PetscInt r[1] = { PetscInt(row) };
  PetscInt c[sz];
  PetscScalar v[sz];
  for (int i = 0; i < sz; ++i) {
    c[i] = cols[i]; 
    v[i] = val.fastAccessDx(i);
  }
  CALL(MatSetValues(la->J, 1, r, sz, c, v, ADD_VALUES));
}

void diag_jacobian_rows(LinAlg* la, GIDs const& rows) {
  PetscInt sz = rows.size();
  PetscInt rs[sz];
  for (int i = 0; i < sz; ++i)
    rs[i] = rows[i];
  CALL(MatZeroRows(la->J, sz, rs, 1.0, PETSC_NULL, PETSC_NULL));
}

void zero_residual(LinAlg* la) {
  CALL(VecSet(la->f, 0.0));
}

void zero_jacobian(LinAlg* la) {
  CALL(MatZeroEntries(la->J));
}

void synchronize(LinAlg* la) {
  CALL(VecAssemblyBegin(la->f));
  CALL(VecAssemblyEnd(la->f));
  CALL(VecAssemblyBegin(la->q));
  CALL(VecAssemblyEnd(la->q));
  CALL(MatAssemblyBegin(la->J, MAT_FLUSH_ASSEMBLY));
  CALL(MatAssemblyEnd(la->J, MAT_FLUSH_ASSEMBLY));
}

void finalize(LinAlg* la) {
  CALL(VecAssemblyBegin(la->f));
  CALL(VecAssemblyEnd(la->f));
  CALL(VecAssemblyBegin(la->q));
  CALL(VecAssemblyEnd(la->q));
  CALL(MatAssemblyBegin(la->J, MAT_FINAL_ASSEMBLY));
  CALL(MatAssemblyEnd(la->J, MAT_FINAL_ASSEMBLY));
}

void solve(LinAlg* la) {
  double t0 = time();
  KSP ksp;
  PC pc;
  CALL(VecScale(la->f, -1.0));
  CALL(KSPCreate(PETSC_COMM_WORLD, &ksp));
  CALL(KSPSetTolerances(ksp, 1.0e-10, 1.0e-10,
        PETSC_DEFAULT, 2000));
  CALL(KSPSetOperators(ksp, la->J, la->J));
  CALL(KSPSetFromOptions(ksp));
  CALL(KSPGetPC(ksp, &pc));
  CALL(PCSetType(pc, PCLU));
  CALL(KSPSolve(ksp, la->f, la->dx));
  double t1 = time();
  print(" > linear system solved in %f seconds", t1-t0);
}

void add_to_primal(LinAlg* la, Disc* d) {
  PetscScalar* dx;
  CALL(VecGetArray(la->dx, &dx));
  int row = 0;
  Vector u(0,0,0);
  Vector du(0,0,0);
  apf::DynamicArray<apf::Node> u_nodes;
  apf::DynamicArray<apf::Node> p_nodes;
  apf::getNodes(d->u_nmbr, u_nodes);
  apf::getNodes(d->p_nmbr, p_nodes);
  for (size_t n = 0; n < u_nodes.size(); ++n) {
    auto node = u_nodes[n];
    auto ent = node.entity;
    auto lnode = node.node;
    if (! d->mesh->isOwned(node.entity)) continue;
    for (int i = 0; i < d->dim; ++i)
      du[i] = dx[row++];
    apf::getVector(d->u, ent, lnode, u);
    apf::setVector(d->u, ent, lnode, u + du);
  }
  for (size_t n = 0; n < p_nodes.size(); ++n) {
    auto node = p_nodes[n];
    auto ent = node.entity;
    auto lnode = node.node;
    if (! d->mesh->isOwned(ent)) continue;
    auto dp = dx[row++];
    auto p = apf::getScalar(d->p, ent, lnode);
    apf::setScalar(d->p, ent, lnode, p + dp);
  }
  apf::synchronize(d->u);
  apf::synchronize(d->p);
}

}
