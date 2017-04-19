#ifndef _TIMME_ASSEMBLY_
#define _TIMME_ASSEMBLY_

#include "Cutoff.h"
#include "SiteID.h"
#include "Assembly.h"

////////////////////////////////////////////////////////////////////////////////
class TIMME_Assembly : public Assembly<Cutoff, SiteID> {

public:
  TIMME_Assembly()
  {};
  ~TIMME_Assembly()
  {};

};

#endif

