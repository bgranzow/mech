#ifndef qoi_hpp
#define qoi_hpp

#include "control.hpp"
#include "disc.hpp"

namespace mech {

template <typename T> struct QoI;

template <>
struct QoI<ST> : public Integrator {
  QoI(Disc* d);
  virtual ~QoI();
  ST const& get_qoi_value() const { return qoi_value; }
  ST const& get_elem_value() const { return elem_value; }
  virtual void set_time(double) {}
  virtual void set_elem_set(std::string const&) {}
  virtual void pre_process(LinAlg*);
  virtual void gather(apf::MeshElement* me);
  virtual void in_elem(apf::MeshElement*) {}
  virtual void at_point(Vector const&, double, double);
  virtual void out_elem() {}
  virtual void scatter(LinAlg* la);
  virtual void post_process(LinAlg* la);
  Disc* disc;
  apf::MeshElement* elem;
  ST qoi_value;
  ST elem_value;
};

template <>
struct QoI<FADT> : public Integrator {
  QoI(Disc* d);
  virtual ~QoI();
  ST const& get_qoi_value() const { return qoi_value; }
  FADT const& get_elem_value() const { return elem_value; }
  virtual void set_time(double) {}
  virtual void set_elem_set(std::string const&) {}
  virtual void pre_process(LinAlg*);
  virtual void gather(apf::MeshElement* me);
  virtual void in_elem(apf::MeshElement*) {}
  virtual void at_point(Vector const&, double, double);
  virtual void out_elem() {}
  virtual void scatter(LinAlg* la);
  virtual void post_process(LinAlg* la);
  Disc* disc;
  apf::MeshElement* elem;
  ST qoi_value;
  FADT elem_value;
};

}

#endif
