#ifndef assembly_hpp
#define assembly_hpp

#include "control.hpp"

namespace mech {

Integrator find_evaluator(Evaluators& E, std::string const& n);
void set_time(Evaluators& E, double t_now);
void assemble(Evaluators& E, Disc* d, LinAlg* la);

}

#endif
