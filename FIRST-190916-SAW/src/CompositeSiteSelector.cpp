#include "CompositeSiteSelector.h"

void CompositeSiteSelector::addSiteSelector(SiteSelector *siteSelector) {
  setOfSiteSelectors.insert(siteSelector);
}

bool CompositeSiteSelector::includeSite(SiteID siteID) {
  for (std::set <SiteSelector *>::iterator selectorIterator = setOfSiteSelectors.begin();
       selectorIterator != setOfSiteSelectors.end();
       selectorIterator++) {
    SiteSelector * siteSelector = *selectorIterator;
    
    if (!siteSelector->includeSite(siteID)) {
      return false;
    } 
  }
  
  return true;
}

