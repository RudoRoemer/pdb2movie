#include "TargetedSiteSelector.h"

#include "Froda.h"

TargetedSiteSelector::TargetedSiteSelector(Froda * const froda) {
  this->froda = froda;
}

bool TargetedSiteSelector::includeSite(SiteID siteID) {
  if (froda->isTargeted) { // FIXME - implement
    
    if (froda->frodaAtom.at(siteID).myTarget > 0) {
      
      return true;
    } else {
      
      //std::cout << "not targeted: " << siteID<< std::endl;
    }
  }
  
  return false;
}
