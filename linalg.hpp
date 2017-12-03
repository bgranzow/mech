#ifndef linalg_hpp
#define linalg_hpp

#include "control.hpp"
#include <petsc.h>

namespace mech {

struct LinAlg {
  LinAlg(Disc* d);
  ~LinAlg();
  Mat J;
  Vec f;
  Vec dx;
  Vec q;
};

void add_to_jacobian(LinAlg* la, GID row, GIDs cols, FADT const& resid);
void add_to_residual(LinAlg* la, GID row, double val);
void set_to_residual(LinAlg* la, GID row, double val);
void diag_jacobian_row(LinAlg* la, GID row);
void zero_residual(LinAlg* la);
void zero_jacobian(LinAlg* la);
void synchronize(LinAlg* la);

}

#endif
