#ifndef linalg.hpp
#define linalg.hpp

#include "control.hpp"
#include <petsc.h>

namespace mech {

struct LinAlg {
  LinAlg(Disc* d);
  Mat J;
  Vec f;
  Vec dx;
  Vec q
};

void add_to_jacobian(LinAlg* la, GID row, GIDs cols, FADT const& resid);

}
