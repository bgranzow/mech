#ifndef mechanics_hpp
#define mechanics_hpp

#include "control.hpp"

namespace mech {

RCP<Integrator> find_evaluator(Evaluators& E, std::string const& n);
void set_time(Evaluators& E, double t_now);
void assemble(Evaluators& E, Disc* d, LinAlg* la);

template <typename T>
void build_resid(Evaluators&E, Disc* d, Input* in, bool save);

}

#endif
