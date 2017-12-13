#include "control.hpp"
#include "displacement.hpp"
#include "disc.hpp"
#include "mechanics.hpp"
#include "pressure.hpp"
#include "qoi.hpp"

namespace mech {

struct Functional {
  Functional(Input* in, Disc* d, RCP<QoI<ST>> q, double t);
  Evaluators evals;
  Input* input;
  Disc* disc;
  RCP<QoI<ST>> qoi;
  double time;
};

Functional::Functional(Input* in, Disc* d, RCP<QoI<ST>> q, double t) :
  input(in), disc(d), qoi(q), time(t) {}

static void construct_evaluators(Functional* f) {
  auto u = rcp(new Displacement<ST>(f->disc, NONE));
  auto p = rcp(new Pressure<ST>(f->disc, NONE));
  f->evals.push_back(u);
  f->evals.push_back(p);
  build_resid<ST>(f->evals, f->input, f->disc, false);
  f->qoi->set_fields(f->evals);
  f->evals.push_back(f->qoi);
  set_time(f->evals, f->time);
}

static void do_evaluation(Functional* f) {
  auto t0 = time();
  set_time(f->evals, f->time);
  assemble(f->evals, f->disc, NULL);
  auto t1 = time();
  print(" > functional computed in %f seconds", t1 - t0);
}

static void print_value(Functional* f) {
  auto n = f->qoi->name;
  auto J = f->qoi->get_qoi_value();
  print(" > functional : %s", n.c_str());
  print(" > J(uH) = %.15e", J);
}

void compute_functional(Input* in, Disc* d, RCP<QoI<ST>> qoi, double t) {
  print("*** functional evaluation");
  print("*** at time: %f", t);
  Functional functional(in, d, qoi, t);
  construct_evaluators(&functional);
  do_evaluation(&functional);
  print_value(&functional);
}

}
