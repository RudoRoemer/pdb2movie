#ifndef COMPOSITE_SITE_SELECTOR_H
#define COMPOSITE_SITE_SELECTOR_H

#include <set>

#include "SiteID.h"
#include "SiteSelector.h"

class CompositeSiteSelector : public SiteSelector {
public:
  
  void addSiteSelector(SiteSelector *);
  bool includeSite(SiteID siteID);
  
private:
  std::set <SiteSelector *> setOfSiteSelectors;
};

#endif // COMPOSITE_SITE_SELECTOR_H
