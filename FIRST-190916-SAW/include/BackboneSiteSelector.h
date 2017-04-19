#ifndef BACKBONE_SITE_SELECTOR_H
#define BACKBONE_SITE_SELECTOR_H

#include "MolFramework.h"
#include "SiteID.h"
#include "SiteSelector.h"

class BackboneSiteSelector : public SiteSelector {

public:
  BackboneSiteSelector(MolFramework *);
  ~BackboneSiteSelector(){};
  bool includeSite(SiteID siteID);
  
private:
    MolFramework * molFramework;
};

#endif // BACKBONE_SITE_SELECTOR_H
