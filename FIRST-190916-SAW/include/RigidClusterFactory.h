#ifndef RIGIDCLUSTERFACTORY_H_
#define RIGIDCLUSTERFACTORY_H_

#include <map>
#include "MolFramework.h" 
#include "Site_Info.h"

#include "Chain.h"
#include "Residue.h"

#include "flexweb.h"
#include "RigidCluster.h"
#include "RigidClusterHierarchy.h"

#include "Cutoff.h"

class RigidClusterFactory { // TODO - depreciate and merge into RigidCluster.h
	
 public:
  RigidClusterFactory();
  ~RigidClusterFactory();
	
  // create and return a new RigidCluster containing the previous largestRigidClusters for firstSiteID and secondSiteID
	
		
  RigidCluster * joinSites(Cutoff cutoff, std::set<SiteID>);
  
  RigidCluster * joinSites(Cutoff cutoff,
                           SiteID firstSiteID, 
                           SiteID secondSiteID);
	  
  void insertRigidSubCluster(RigidCluster * rigidCluster,
                             Cutoff cutoff,
                             RigidCluster * rigidSubCluster);
    
  RigidCluster * joinRigidClusters(Cutoff cutoff,
                                   RigidCluster * firstRigidCluster, 
                                   RigidCluster * secondRigidCluster);
	
  RigidCluster * getLargestRigidClusterContainingSiteID(SiteID siteID);
  RigidCluster * getLargestRigidClusterContainingRigidCluster(RigidCluster *,
                                                              unsigned int recursionDepth = 0);
  	
  RigidCluster * decomposeIntoRigidClusters(std::multimap<Cutoff, SiteID * >);
  RigidCluster * decomposeIntoRigidClustersOLD(std::multimap<Cutoff, std::set<SiteID> >);
	
  RigidClusterHierarchy * getRigidClusterHierarchy();
  
private:
  RigidClusterHierarchy * rigidClusterHierarchy;
    							
};

#endif /*RIGIDCLUSTERFACTORY_H_*/
