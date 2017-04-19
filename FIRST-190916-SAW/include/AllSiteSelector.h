#ifndef _ALL_SITE_SELECTOR_
#define _ALL_SITE_SELECTOR_

#include "SiteID.h"
#include "SiteSelector.h"

// a SiteSelector that always returns true for includeSite(siteID), given any 
// valid any siteID
// 
////////////////////////////////////////////////////////////////////////////////
class AllSiteSelector : public SiteSelector {
  bool includeSite(SiteID siteID);

 public:
  AllSiteSelector(){};
  ~AllSiteSelector(){};
};

#endif // _ALL_SITE_SELECTOR_
