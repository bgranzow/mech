#include <control.hpp>
#include <displacement.hpp>
#include <finedisc.hpp>
#include <mechanics.hpp>
#include <qoi.hpp>

namespace mech {

double zero(Vector const&, double) { return 0; }
double tenth(Vector const&, double) { return 0.1; }

template <typename T>
struct AvgDisp : public QoI<T> {
  AvgDisp(Disc* d) : QoI<T>(d) {
    this->disc = d;
    dim = this->disc->dim;
    this->name = "avg disp";
  }
  void set_fields(Evaluators& E) {
    auto disp = find_evaluator(E, "u");
    u = rcp_static_cast<Displacement<T>>(disp);
  }
  void at_point(Vector const&, double w, double dv) {
    for (int i = 0; i < dim; ++i)
      this->elem_value += u->val(i) * w * dv;
    this->elem_value /= dim;
  }
  int dim;
  RCP<Displacement<T>> u;
};

static void setup(Input* in, char** argv) {
  in->model = ELASTIC;
  in->geom_file = argv[1];
  in->mesh_file = argv[2];
  in->assoc_file = argv[3];
  Material mat1 = { 1000.0, 0.25, 100.0, 10.0 };
  in->mats["body"] = mat1;
  DBC dbc1 = { 0, zero, "xmin" };
  DBC dbc2 = { 1, zero, "ymin" };
  DBC dbc3 = { 0, tenth, "xmax" };
  in->dbcs = { dbc1, dbc2, dbc3 };
}

static void output(FineDisc* d, std::string const& n) {
  auto zu = apf::createFieldOn(d->mesh, "zu2", apf::VECTOR);
  apf::projectField(zu, d->zu);
  apf::writeVtkFiles(n.c_str(), d->mesh);
  apf::destroyField(zu);
}

static void do_primal_solve(Input* in, Disc* d) {
  build_disc_data(d);
  solve_linear_primal(in, d, 0.0);
  free_disc_data(d);
}

static void do_functional_eval(Input* in, Disc* d) {
  build_disc_data(d);
  auto qoi = rcp(new AvgDisp<ST>(d));
  compute_functional(in, d, qoi, 0.0);
  free_disc_data(d);
}

static void do_adjoint_solve(Input* in, Disc* d) {
  FineDisc fd;
  init_fine_disc(&fd, d);
  build_disc_data(&fd);
  auto qoi = rcp(new AvgDisp<FADT>(&fd));
  solve_adjoint(in, &fd, qoi, 0.0);
  output(&fd, "out_ex7");
  free_disc_data(&fd);
  free_fine_disc(&fd);
}

static void run(char** argv) {
  Input in;
  Disc disc;
  setup(&in, argv);
  init_disc(&disc, &in);
  do_primal_solve(&in, &disc);
  do_functional_eval(&in, &disc);
  do_adjoint_solve(&in, &disc);
  free_disc(&disc);
}

template struct AvgDisp<ST>;
template struct AvgDisp<FADT>;

}

int main(int argc, char** argv) {
  mech::init(&argc, &argv);
  mech::run(argv);
  mech::free();
}
