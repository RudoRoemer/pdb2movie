#include "RigidCluster.h"
#include <algorithm>

void RigidCluster::initializeRigidCluster() {	
}

RigidCluster::RigidCluster() {
	initializeRigidCluster();
}

RigidCluster::RigidCluster(unsigned int siteID) {
	initializeRigidCluster();

	insertSiteID(siteID);
};

RigidCluster::RigidCluster(const RigidCluster &rigidCluster) {
	subClusters = rigidCluster.subClusters;
	memberSiteIDs = rigidCluster.memberSiteIDs;
};

RigidCluster & RigidCluster::operator=(const RigidCluster &rigidCluster) {
	subClusters = rigidCluster.subClusters;
	memberSiteIDs = rigidCluster.memberSiteIDs;

	return *this;
}

RigidCluster::SiteIdIterator RigidCluster::beginSiteIDs() {
	return this->memberSiteIDs.begin();
}; 

RigidCluster::SiteIdIterator RigidCluster::endSiteIDs() {
	return this->memberSiteIDs.end();	
};

RigidCluster::SiteIdReverseIterator RigidCluster::rbeginSiteIDs() {
  return this->memberSiteIDs.rbegin();
};

RigidCluster::SiteIdReverseIterator RigidCluster::rendSiteIDs() {
  return this->memberSiteIDs.rend();
};

RigidCluster::RigidClusterIterator RigidCluster::beginSubClusters() {
	return this->subClusters.begin();
}; 

RigidCluster::RigidClusterIterator RigidCluster::endSubClusters() {
	return this->subClusters.end();	
};

void RigidCluster::insertSiteID(unsigned int siteID) {
	this->memberSiteIDs.insert(siteID);
};

void RigidCluster::insertSiteIDs(const RigidCluster::SiteIdIterator begin, RigidCluster::SiteIdIterator end) {
	this->memberSiteIDs.insert(begin, end);
};

void RigidCluster::insertSubCluster(RigidCluster * subCluster) {
  // remove any sites already contained in the subCluster  
  //set_difference(inserter(this->memberSiteIDs));
  
  
  set_difference(subCluster->memberSiteIDs.begin(), subCluster->memberSiteIDs.end(), 
                 memberSiteIDs.begin(), memberSiteIDs.end(),
                 inserter(subCluster->memberSiteIDs, subCluster->memberSiteIDs.begin()));
  
/* // TODO - remove any sites that are common to 
  this->memberSiteIDs.erase(subCluster->memberSiteIDs.begin(),
                            subCluster->memberSiteIDs.end());*/
  
  // add the subCluster to this cluster's set of SiteIDs
	this->subClusters.push_back(subCluster);
}

size_t RigidCluster::size() {
/*	size_t memberSiteIDsSize = this->memberSiteIDs.size();

	size_t totalSubClusterSize = 0;

	for (RigidCluster::RigidClusterIterator subclusterIterator = beginSubClusters();
			subclusterIterator != endSubClusters();
			subclusterIterator ++) {
		RigidCluster * subCluster = *subclusterIterator;

		totalSubClusterSize += subCluster->size();
	}

	size_t totalSize = memberSiteIDsSize + totalSubClusterSize;

	return totalSize;*/
  
  AllSiteSelector *allSiteSelector = new AllSiteSelector();
  size_t totalSize = size(allSiteSelector);

  delete allSiteSelector;
  return totalSize;
}

// Description:
// return a count of the member sites for which the function pointer selector returns true
// 
size_t RigidCluster::size(SiteSelector * siteSelector) {  
  size_t size = 0;
  
  for (std::set<SiteID>::iterator siteIDiterator = memberSiteIDs.begin();
       siteIDiterator != memberSiteIDs.end();
       siteIDiterator++) {
    SiteID siteID = *siteIDiterator;
    
    if (siteSelector->includeSite(siteID)) {
      size++;
    }
  }
  
	for (RigidCluster::RigidClusterIterator subclusterIterator = beginSubClusters();
       subclusterIterator != endSubClusters();
       subclusterIterator ++) {
		RigidCluster * subCluster = *subclusterIterator;
    
		size += subCluster->size(siteSelector);
	}
  
 	return size;
}

RigidCluster::SiteIdIterator RigidCluster::findSiteID(const SiteID& siteID) {
	return this->memberSiteIDs.find(siteID); // FIXME - find a way to identify site IDs 
};

size_t RigidCluster::eraseSiteID(const SiteID& siteID) {
	return this->memberSiteIDs.erase(siteID);
};

bool operator<(RigidCluster const & firstRigidCluster, const RigidCluster & secondRigidCluster) {
	return (firstRigidCluster.memberSiteIDs.size() < secondRigidCluster.memberSiteIDs.size());	
}

bool operator!=(RigidCluster const & firstRigidCluster, const RigidCluster & secondRigidCluster) {
	return (
			(firstRigidCluster.memberSiteIDs != secondRigidCluster.memberSiteIDs) &&
			(firstRigidCluster.subClusters   != secondRigidCluster.subClusters) );	
}

bool operator==(RigidCluster const & firstRigidCluster, const RigidCluster & secondRigidCluster) {
	return (
			(firstRigidCluster.memberSiteIDs == secondRigidCluster.memberSiteIDs) &&
			(firstRigidCluster.subClusters   == secondRigidCluster.subClusters) );	
}

