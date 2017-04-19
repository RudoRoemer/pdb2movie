#ifndef TARGETED_SITE_SELECTOR_H
#define TARGETED_SITE_SELECTOR_H

#include <vector>

#include "SiteID.h"
#include "SiteSelector.h"
//#include "Froda.h"

class Froda;

class TargetedSiteSelector : public SiteSelector {
public:
  TargetedSiteSelector(Froda * const);
  ~TargetedSiteSelector(){};
  bool includeSite(SiteID siteID);
  
private:
  Froda * froda; // FIXME - why won't this work? (circular includes?)
};

#endif // TARGETED_SITE_SELECTOR_H
