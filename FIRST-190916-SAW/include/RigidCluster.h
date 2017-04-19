#ifndef RIGIDCLUSTER_H_
#define RIGIDCLUSTER_H_

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <utility>

#include "AllSiteSelector.h"
#include "RigidLabel.h"
#include "SiteID.h"
#include "SiteSelector.h"

// TODO - does it make any sense to make this generic?
class RigidCluster {
 public:
	
  typedef std::vector<RigidCluster *>::iterator RigidClusterIterator;
  typedef std::set<SiteID>::iterator SiteIdIterator;
  typedef std::set<SiteID>::reverse_iterator SiteIdReverseIterator;

  RigidCluster();
  virtual ~RigidCluster(){};
	
  void initializeRigidCluster();

  RigidCluster(unsigned int siteID);
	
  RigidCluster(const RigidCluster &rigidCluster);
	
  RigidCluster & operator=(const RigidCluster &rigidCluster);
	
  RigidLabel rigidLabel;
  // TODO - deallocate all references (smart pointers / garbage collection would be very nice here :-)
		
  SiteIdIterator beginSiteIDs(); 
  SiteIdIterator endSiteIDs();

  SiteIdReverseIterator rbeginSiteIDs();
  SiteIdReverseIterator rendSiteIDs();

  RigidClusterIterator beginSubClusters(); 
	
  RigidClusterIterator endSubClusters();
	
  void insertSiteID(unsigned int siteID);
	
  void insertSiteIDs(const SiteIdIterator begin, SiteIdIterator end);
			
  void insertSubCluster(RigidCluster * subCluster);
	
  size_t size();
  size_t size(SiteSelector*);
  
  SiteIdIterator findSiteID(const SiteID& siteID);

  size_t eraseSiteID(const SiteID& siteID);
	 
  friend bool operator<(RigidCluster const & firstRigidCluster, const RigidCluster & secondRigidCluster);
  friend bool operator!=(RigidCluster const & firstRigidCluster, const RigidCluster & secondRigidCluster);
  friend bool operator==(RigidCluster const & firstRigidCluster, const RigidCluster & secondRigidCluster);
	
	// TODO - make these private
	unsigned int minLabel;
	unsigned int maxLabel;
	
//private:	// FIXME - make the following private
	std::vector<RigidCluster *> subClusters; // FIXME - go back to std::container<RigidCluster *> or (better, if it would work) std::container<RigidCluster >
	std::set<SiteID> memberSiteIDs;
  
};


#endif /*RIGIDCLUSTER_H_*/
