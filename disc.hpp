#ifndef disc_hpp
#define disc_hpp

#include "control.hpp"
#include <apfAlbany.h>
#include <apfMesh2.h>

namespace mech {

struct Disc {
  int dim;
  int elem_type;
  int num_u_elem_nodes;
  int num_p_elem_nodes;
  int num_elem_dofs;
  int q_degree;
  int num_ips;
  int num_owned_dofs;
  int num_total_dofs;
  int p_dof_offset;
  apf::Mesh2* mesh;
  apf::StkModels* sets;
  apf::FieldShape* u_basis;
  apf::FieldShape* p_basis;
  apf::Field* u;
  apf::Field* p;
  apf::Field* first_pk;
  apf::Field* eqps;
  apf::Field* eqps_old;
  apf::Field* Fp;
  apf::Field* Fp_old;
  apf::GlobalNumbering* u_nmbr;
  apf::GlobalNumbering* p_nmbr;
  apf::ElemSets elem_sets;
  apf::SideSets side_sets;
  apf::NodeSets u_node_sets;
  apf::NodeSets p_node_sets;
};

void init_disc(Disc* d, Input* in);
void free_disc(Disc* d);

void build_disc_data(Disc* d);
void free_disc_data(Disc* d);

GID get_u_gid(Disc* d, apf::MeshEntity* e, int n, int i);
GID get_p_gid(Disc* d, apf::MeshEntity* e, int n);
GID get_u_gid(Disc* d, apf::Node const& n, int i);
void get_gids(Disc* d, apf::MeshEntity* e, GIDs& ids);

}

#endif
