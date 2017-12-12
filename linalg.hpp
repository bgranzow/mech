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

void add_to_jacobian(LinAlg* la, GID row, GIDs const& cols, FADT const& resid);
void add_to_residual(LinAlg* la, GID row, double val);
void set_to_residual(LinAlg* la, GID row, double val);
void add_to_functional(LinAlg* la, GID row, double val);
void set_to_functional(LinAlg* la, GID row, double val);
void diag_jacobian_rows(LinAlg* la, GIDs const& rows);
void zero_residual(LinAlg* la);
void zero_functional(LinAlg* la);
void zero_jacobian(LinAlg* la);
void synchronize(LinAlg* la);
void finalize(LinAlg* la);
void solve_primal_sys(LinAlg* la);
void solve_adjoint_sys(LinAlg* la);
void add_to_primal(LinAlg* la, Disc* d);
double get_resid_norm(LinAlg* la);

}

#endif
