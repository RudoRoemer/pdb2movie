#ifndef SITE_SELECTOR_H
#define SITE_SELECTOR_H

#include "SiteID.h"

class SiteSelector {
public:
  virtual bool includeSite(SiteID siteID);
  virtual ~SiteSelector();
};
#endif // SITE_SELECTOR_H
