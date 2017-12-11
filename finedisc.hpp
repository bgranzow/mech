#ifndef finedisc_hpp
#define finedisc_hpp

#include "disc.hpp"

namespace mech {

struct FineDisc : public Disc {
  apf::Field* zu;
  apf::Field* zp;
};

void init_fine_disc(FineDisc* fd, Disc* d);
void free_fine_disc(FineDisc* fd);

}

#endif
