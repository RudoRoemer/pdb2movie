#ifndef RIGID_CLUSTER_HIERARCHY_H_
#define RIGID_CLUSTER_HIERARCHY_H_

#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "RigidCluster.h"
#include "Cutoff.h"
#include "Color.h"

using namespace std;

class RigidClusterHierarchy {

public:
  std::map<Cutoff, std::set<RigidCluster *> > mapFromCutoffToSetOfLargestRigidClusters;
  
  std::map<unsigned int, unsigned int> mapFromSiteIdToRigidLabel;
  std::map<unsigned int, unsigned int> mapFromRigidLabelToSiteId;
  std::map<unsigned int, RigidCluster *> mapFromRigidLabelToCluster;

  std::map<unsigned int, ColorNumber> mapFromRigidClusterLabelToColorNumber;
  
  std::map<SiteID, SiteID> mapFromRigidLabelToMeanSiteID;
  std::map<SiteID, SiteID> mapFromRigidLabelToRMSDfromMeanSiteID;
  
  std::map<SiteID, SiteID> mapFromMaxLabelToMeanSiteID; // FIXME - removeme
  std::map<SiteID, SiteID> mapFromMaxLabelToRMSDfromMeanSiteID; // FIXME - removeme
  
  std::map<RigidCluster *, RigidCluster *> mapFromRigidClusterToLargestKnown; 

  std::set<std::pair<size_t, unsigned int> > setOfSizeAndRigidLabelPairs;
  std::set<RigidCluster *> setOfLargestRigidClusters;

  std::set<Cutoff> setOfCutoffs;
  //std::set<float> setOfAverageRs; 
  std::map<unsigned int, RigidCluster *> mapFromSiteIDToRigidClusterContainingSiteID;
  std::map<RigidCluster *, RigidCluster *> mapFromRigidSubClusterToParentRigidCluster;
  std::map<RigidCluster *, Cutoff> mapFromRigidClusterToCutoff;


  void printCutoffs(){

    ofstream oldCutoffs;
    oldCutoffs.open( "oldCutoffs.txt" );
    
    map< RigidCluster*, Cutoff>::iterator mapFromRigidClusterToCutoffIter = mapFromRigidClusterToCutoff.begin();
    while( mapFromRigidClusterToCutoffIter != mapFromRigidClusterToCutoff.end() ){
      oldCutoffs << showpoint << setw(16) << setprecision(16) << mapFromRigidClusterToCutoffIter->second << endl;      
      mapFromRigidClusterToCutoffIter++;
    }

  }

};

#endif /* RIGID_CLUSTER_HIERARCHY_H_ */
