#include "control.hpp"
#include "displacement.hpp"
#include "finedisc.hpp"
#include "linalg.hpp"
#include "mechanics.hpp"
#include "pressure.hpp"
#include "qoi.hpp"

namespace mech {

struct Adjoint {
  Adjoint(Input* in, FineDisc* d, RCP<QoI<FADT>> q, double t);
  Evaluators evals;
  LinAlg la;
  Input* input;
  FineDisc* disc;
  RCP<QoI<FADT>> qoi;
  double time;
};

Adjoint::Adjoint(Input* in, FineDisc* d, RCP<QoI<FADT>> q, double t) :
  la(d), input(in), disc(d), qoi(q), time(t) {}

static void construct_evaluators(Adjoint* adjoint) {
  auto u = rcp(new Displacement<FADT>(adjoint->disc, ADJOINT));
  auto p = rcp(new Pressure<FADT>(adjoint->disc, ADJOINT));
  adjoint->evals.push_back(u);
  adjoint->evals.push_back(p);
  build_resid<FADT>(adjoint->evals, adjoint->input, adjoint->disc, false);
  adjoint->qoi->set_fields(adjoint->evals);
  adjoint->evals.push_back(adjoint->qoi);
  set_time(adjoint->evals, adjoint->time);
}

static void compute_adjoint(Adjoint* adjoint) {
  auto t0 = time();
  zero_residual(&(adjoint->la));
  zero_functional(&(adjoint->la));
  zero_jacobian(&(adjoint->la));
  assemble(adjoint->evals, adjoint->disc, &(adjoint->la));
  finalize(&(adjoint->la));
  set_jacob_dbcs(adjoint->input, adjoint->disc, &(adjoint->la), adjoint->time);
  finalize(&(adjoint->la));
  auto t1 = time();
  print(" > adjoint computed in %f seconds", t1 - t0);
}

void solve_adjoint(Input* in, FineDisc* d, RCP<QoI<FADT>> qoi, double t) {
  print("*** adjoint solve");
  print("*** at time: %f", t);
  print("*** dofs: %lu", d->num_total_dofs);
  Adjoint adjoint(in, d, qoi, t);
  construct_evaluators(&adjoint);
  compute_adjoint(&adjoint);
  solve_adjoint_sys(&(adjoint.la));
}

}
