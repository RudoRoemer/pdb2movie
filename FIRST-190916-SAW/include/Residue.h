#ifndef _RESIDUE_
#define _RESIDUE_

#include <vector>
#include "Site_Info.h"

class Residue {
public: 
  typedef std::vector<Site_Info *> Sites;
  typedef Sites::iterator SiteIterator;

  void addSite(Site_Info * site);
  Sites::iterator beginSites();
  Sites::iterator endSites();
  size_t size();
  
private:
	Sites memberSites;
};
#endif // _RESIDUE_
