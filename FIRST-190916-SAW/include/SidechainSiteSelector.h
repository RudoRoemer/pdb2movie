#ifndef SIDECHAIN_SITE_SELECTOR_H
#define SIDECHAIN_SITE_SELECTOR_H

#include "MolFramework.h"
#include "SiteID.h"
#include "SiteSelector.h"

class SidechainSiteSelector : public SiteSelector {

public:
  SidechainSiteSelector(MolFramework * mol_Framework);
  ~SidechainSiteSelector(){};
  bool includeSite(SiteID siteID);
  
private:
  MolFramework * molFramework;
};

#endif // SIDECHAIN_SITE_SELECTOR_H
