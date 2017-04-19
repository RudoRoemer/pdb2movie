#include "Residue.h"

void Residue::addSite(Site_Info * site) {
  memberSites.push_back(site);
}

Residue::SiteIterator Residue::beginSites() {
  return memberSites.begin();
}

Residue::SiteIterator Residue::endSites() {
  return memberSites.end();
}

size_t Residue::size() {
  return memberSites.size();
}
