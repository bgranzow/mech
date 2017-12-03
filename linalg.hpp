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
void synchronize(LinAlg* la);

}

#endif
