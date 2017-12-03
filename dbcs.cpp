#include "control.hpp"
#include "disc.hpp"
#include "linalg.hpp"

namespace mech {

void set_resid_dbcs(Input* in, LinAlg* la) {
  (void)in;
  (void)la;
}

void set_jacob_dbcs(Input* in, LinAlg* la) {
  for (size_t i = 0; i < in->dbcs.size(); ++i) {
    auto& dbc = in->dbcs[i];
    std::cout << dbc.set << std::endl;
  }
  (void)la;
}

}
