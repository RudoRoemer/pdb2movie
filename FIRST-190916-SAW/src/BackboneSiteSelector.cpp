#include "BackboneSiteSelector.h"

BackboneSiteSelector::BackboneSiteSelector(MolFramework * molFramework) {
  this->molFramework = molFramework;
}

bool BackboneSiteSelector::includeSite(SiteID siteID) {
  if (siteID > molFramework->total_sites) {
    // siteID is out of range (too large)

    return false;
  }
  
  if (molFramework->isMainchain(siteID) > 0) {
    return true;
  }
  
  return false;
}
