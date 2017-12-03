#ifndef mechanics_hpp
#define mechanics_hpp

#include "control.hpp"

namespace mech {

RCP<Integrator> find_evaluator(Evaluators& E, std::string const& n);

void set_time(Evaluators& E, double t_now);

void assemble(Evaluators& E, Disc* d, LinAlg* la);

template <typename T>
void build_resid(Evaluators&E, Input* in, Disc* d, bool save);

void solve_linear_primal(Input* in, Disc* d, double t);

void set_resid_dbcs(Input* in, Disc* d, LinAlg* la, double t);

void set_jacob_dbcs(Input* in, Disc* d, LinAlg* la, double t);

}

#endif
