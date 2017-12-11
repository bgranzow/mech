#include "elastic.hpp"
#include "J2.hpp"
#include "kinematics.hpp"
#include "mechanics.hpp"
#include "mresidual.hpp"
#include "neohookean.hpp"
#include "first_pk.hpp"
#include "presidual.hpp"

namespace mech {

template <typename T>
void build_resid(Evaluators& E, Input* in, Disc* d, bool save) {
  auto u = find_evaluator(E, "u");
  auto p = find_evaluator(E, "p");
  auto kin = rcp(new Kinematics<T>(u));
  RCP<Model<T>> cm;
  if (in->model == ELASTIC)
    cm = rcp(new Elastic<T>(u, d, in));
  else if (in->model == NEOHOOKEAN)
    cm = rcp(new Neohookean<T>(kin, d, in));
  else if (in->model == PLASTIC)
    cm = rcp(new J2<T>(kin, d, in, save));
  else
    fail("unknown model: %d", in->model);
  auto first_pk = rcp(new FirstPK<T>(cm, kin, p, d, save));
  auto mresid = rcp(new MResidual<T>(first_pk, u, d));
  auto presid = rcp(new PResidual<T>(u, p, kin, cm, d, in));
  E.push_back(kin);
  E.push_back(cm);
  E.push_back(first_pk);
  E.push_back(mresid);
  E.push_back(presid);
}

template void build_resid<ST>(Evaluators& E, Input* in, Disc* d, bool save);
template void build_resid<FADT>(Evaluators& E, Input* in, Disc* d, bool save);

}
