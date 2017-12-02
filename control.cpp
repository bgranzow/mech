#include "control.hpp"
#include <cstdarg>
#include <cstdlib>
#include <PCU.h>
#include <petsc.h>

namespace mech {

void init(int* argc, char*** argv) {
  MPI_Init(argc, argv);
  PCU_Comm_Init();
  PetscInitialize(argc, argv, NULL, "intro");
}

void free() {
  PetscFinalize();
  PCU_Comm_Free();
  MPI_Finalize();
}

void print(const char* message, ...) {
  if (PCU_Comm_Self()) return void();
  va_list ap;
  va_start(ap, message);
  vfprintf(stdout, message, ap);
  va_end(ap);
  printf("\n");
}

void fail(const char* why, ...) {
  va_list ap;
  va_start(ap, why);
  vfprintf(stderr, why, ap);
  va_end(ap);
  printf("\n");
  abort();
}

double time() {
  return PCU_Time();
}

Integrator::~Integrator() {
}

}
