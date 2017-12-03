#include "control.hpp"
#include "displacement.hpp"
#include "disc.hpp"
#include "linalg.hpp"
#include "mechanics.hpp"
#include "pressure.hpp"

namespace mech {

struct Primal {
  Primal(Input* in, Disc* d);
  Evaluators resid;
  Evaluators jacob;
  LinAlg la;
  Input* input;
  Disc* disc;
};

Primal::Primal(Input* in, Disc* d) : la(d), input(in), disc(d) {}

static void construct_resid(Primal* primal) {
  auto u = rcp(new Displacement<ST>(primal->disc, PRIMAL));
  auto p = rcp(new Pressure<ST>(primal->disc, PRIMAL));
  primal->resid.push_back(u);
  primal->resid.push_back(p);
  build_resid<ST>(primal->resid, primal->input, primal->disc, true);
}

static void construct_jacob(Primal* primal) {
  auto u = rcp(new Displacement<FADT>(primal->disc, PRIMAL));
  auto p = rcp(new Pressure<FADT>(primal->disc, PRIMAL));
  primal->jacob.push_back(u);
  primal->jacob.push_back(p);
  build_resid<FADT>(primal->jacob, primal->input, primal->disc, true);
}

static void construct_primal(Primal* primal, double t) {
  construct_resid(primal);
  construct_jacob(primal);
  set_time(primal->resid, t);
  set_time(primal->jacob, t);
}

#if 0
static void compute_residual(Primal* primal, double t) {
  auto t0 = time();
  zero_residual(&(primal->la));
  assemble(primal->resid, primal->disc, &(primal->la));
  set_resid_dbcs(primal->input, &(primal->la));
  auto t1 = time();
  print(" > residual computed in %f seconds", t1 - t0);
}
#endif

static void compute_jacobian(Primal* primal, double t) {
  auto t0 = time();
  zero_residual(&(primal->la));
  zero_jacobian(&(primal->la));
  assemble(primal->jacob, primal->disc, &(primal->la));
  synchronize(&(primal->la));
  set_jacob_dbcs(primal->input, primal->disc, &(primal->la), t);
  synchronize(&(primal->la));
  auto t1 = time();
  print(" > jacobian computed in %f seconds", t1 - t0);
}

void solve_linear_primal(Input* in, Disc* d, double t) {
  print("*** primal solve");
  print("*** at time: %f", t);
  print("*** dofs: %lu", d->num_total_dofs);
  Primal primal(in, d);
  construct_primal(&primal, t);
  compute_jacobian(&primal, t);
}

}
