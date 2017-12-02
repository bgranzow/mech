#include "disc.hpp"
#include "mechanics.hpp"

namespace mech {

RCP<Integrator> find_evaluator(Evaluators& E, std::string const& n) {
  for (size_t i = 0; i < E.size(); ++i)
    if (E[i]->name == n)
      return E[i];
  fail("integrator %s not found", n.c_str());
}

void set_time(Evaluators& E, double t_now) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->set_time(t_now);
}

void pre_process(Evaluators& E, LinAlg* la) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->pre_process(la);
}

void set_elem_set(Evaluators& E, std::string const& n) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->set_elem_set(n);
}

void gather(Evaluators& E, apf::MeshElement* me) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->gather(me);
}

void in_elem(Evaluators& E, apf::MeshElement* me) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->in_elem(me);
}

void at_point(Evaluators& E, Vector const& p, double w, double dv) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->at_point(p, w, dv);
}

void out_elem(Evaluators& E) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->out_elem();
}

void scatter(Evaluators& E, LinAlg* la) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->scatter(la);
}

void post_process(Evaluators& E, LinAlg* la) {
  for (size_t i = 0; i < E.size(); ++i)
    E[i]->post_process(la);
}

void assemble(Evaluators& E, Disc* d, LinAlg* la) {
  apf::Vector3 xi(0,0,0);
  pre_process(E, la);
  for (auto& set : d->elem_sets) {
    auto es_name = set.first;
    auto elems = set.second;
    set_elem_set(E, es_name);
    for (size_t i = 0; i < elems.size(); ++i) {
      auto elem = elems[i];
      auto me = apf::createMeshElement(d->mesh, elem);
      gather(E, me);
      in_elem(E, me);
      for (int ip = 0; ip < d->num_ips; ++ip) {
        apf::getIntPoint(me, d->q_degree, ip, xi);
        auto dv = apf::getDV(me, xi);
        auto w = apf::getIntWeight(me, d->q_degree, ip);
        at_point(E, xi, w, dv);
      }
      out_elem(E);
      scatter(E, la);
      apf::destroyMeshElement(me);
    }
  }
  post_process(E, la);
}

}
