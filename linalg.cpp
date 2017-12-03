#include "disc.hpp"
#include "linalg.hpp"

#define CALL(function) assert(0 == (function))

namespace mech {

LinAlg::LinAlg(Disc* d) {
  PetscInt n = d->num_owned_dofs;
  PetscInt N = d->num_total_dofs;
  CALL(VecCreateMPI(PETSC_COMM_WORLD, n, N, &f));
  CALL(VecCreateMPI(PETSC_COMM_WORLD, n, N, &dx));
  CALL(VecCreateMPI(PETSC_COMM_WORLD, n, N, &q));
  CALL(MatCreateAIJ(PETSC_COMM_WORLD, n, n, N, N,
        500, PETSC_NULL, 500, PETSC_NULL, &J));
  CALL(MatSetOption(J, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
}

LinAlg::~LinAlg() {
  CALL(MatDestroy(&J));
  CALL(VecDestroy(&f));
  CALL(VecDestroy(&dx));
  CALL(VecDestroy(&q));
}

void add_to_residual(LinAlg* la, GID row, double val) {
  PetscInt r[1] = { row };
  PetscScalar v[1] = { val };
  CALL(VecSetValues(la->f, 1, r, v, ADD_VALUES));
}

void add_to_jacobian(LinAlg* la, GID row, GIDs cols, FADT const& val) {
  int sz = cols.size();
  PetscInt r[1] = { row };
  PetscInt c[sz];
  PetscScalar v[sz];
  for (int i = 0; i < sz; ++i) {
    c[i] = cols[i]; 
    v[i] = val.fastAccessDx(i);
  }
  CALL(MatSetValues(la->J, 1, r, sz, c, v, ADD_VALUES));
}

void diag_mat_row(LinAlg* la, GID row) {
  PetscInt r[1] = { row };
  CALL(MatZeroRows(la->J, 1, r, 1.0, PETSC_NULL, PETSC_NULL));
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
  CALL(MatAssemblyBegin(la->J, MAT_FINAL_ASSEMBLY));
  CALL(MatAssemblyEnd(la->J, MAT_FINAL_ASSEMBLY));
}

}
