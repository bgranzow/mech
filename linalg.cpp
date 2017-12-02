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
  CALL(MatSetOption(J, MAT_NEWNONZERO_ALLOCATION_ERR, PETSC_FALSE));
}

LinAlg::~LinAlg() {
  CALL(MatDestroy(&J));
  CALL(VecDestroy(&f));
  CALL(VecDestroy(&dx));
  CALL(VecDestroy(&q));
}

}
