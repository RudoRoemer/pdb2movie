#include <algorithm>
#include <iostream>
#include <limits>

#include "RigidClusterFactory.h"
#include "generalUtils.h"
#include "global_defs.h"
#include "Color.h"
#include "RigidLabel.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
RigidClusterFactory::RigidClusterFactory() {

  rigidClusterHierarchy = new RigidClusterHierarchy();
}

////////////////////////////////////////////////////////////////////////////////
RigidClusterFactory::~RigidClusterFactory() {
  // TODO - implement destructor 
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Create and return a new RigidCluster containing the previous 
//   largestRigidClusters for firstSiteID and secondSiteID
////////////////////////////////////////////////////////////////////////////////
RigidCluster *RigidClusterFactory::joinSites(Cutoff cutoff,
					     SiteID firstSiteID, 
					     SiteID secondSiteID) {  
  
  RigidCluster *firstRigidCluster  = getLargestRigidClusterContainingSiteID(firstSiteID);
  RigidCluster *secondRigidCluster = getLargestRigidClusterContainingSiteID(secondSiteID);
  
  return joinRigidClusters( cutoff,
                            firstRigidCluster, 
                            secondRigidCluster );
}

////////////////////////////////////////////////////////////////////////////////
bool sortClustersBySize_clusterFactory (RigidCluster * firstRigidCluster, 
                                        RigidCluster * secondRigidCluster) {

  return firstRigidCluster->size() > secondRigidCluster->size();
}

////////////////////////////////////////////////////////////////////////////////
RigidClusterHierarchy * RigidClusterFactory::getRigidClusterHierarchy() {
  return rigidClusterHierarchy;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
RigidCluster *RigidClusterFactory::decomposeIntoRigidClusters(std::multimap<Cutoff, SiteID *> multimapFromCutoffToJoinedSites) {

  RigidCluster *parentRigidCluster = NULL;

  // Output some text for Flexweb
  parameters.mapFromTaskNameToStatus["Decompose into rigid clusters"] = "Running...";
  outputStatus();
	
  unsigned int totalNumberOfCutoffs = multimapFromCutoffToJoinedSites.size(); // Number of unique keys, or energy cutoffs. 
  unsigned int currentNumberOfCutoffs = 0;
  
  float previousPercentComplete = 0;
  
  // cycle over each cutoff - rigid cluster pair.
  multimap<Cutoff, SiteID * >::iterator pairIterator = multimapFromCutoffToJoinedSites.begin();
  while( pairIterator != multimapFromCutoffToJoinedSites.end() ){    
    
    Cutoff cutoff = (*pairIterator).first;
    SiteID *setOfJoinedSiteIDs = (*pairIterator).second;

    RigidCluster *firstRigidCluster  = getLargestRigidClusterContainingSiteID(setOfJoinedSiteIDs[0]);
    RigidCluster *secondRigidCluster = getLargestRigidClusterContainingSiteID(setOfJoinedSiteIDs[1]);
    
    joinRigidClusters( cutoff,
		       firstRigidCluster, 
		       secondRigidCluster );
    
    pairIterator++;

    currentNumberOfCutoffs++;
    if(parameters.flexweb) {
      // compute and update percentComplete for flexweb status 
      float percentComplete = 100*(float) currentNumberOfCutoffs / totalNumberOfCutoffs;
      if (percentComplete - previousPercentComplete > 5) {
	stringstream percentCompleteString;
	percentCompleteString << "" << noshowpoint << fixed << setprecision(0) << percentComplete << "% Complete";
	parameters.mapFromTaskNameToStatus["Decompose into rigid clusters"] = percentCompleteString.str();
	outputStatus();		
	
	previousPercentComplete = percentComplete;
      }
    }
  }
  
  parameters.mapFromTaskNameToStatus["Decompose into rigid clusters"] = "Complete";
  outputStatus();
  
  return parentRigidCluster;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
RigidCluster *RigidClusterFactory::decomposeIntoRigidClustersOLD(std::multimap<Cutoff, std::set<SiteID> > multimapFromCutoffToJoinedSites) {

  RigidCluster *parentRigidCluster = NULL;

  // Output some text for Flexweb
  parameters.mapFromTaskNameToStatus["Decompose into rigid clusters"] = "Running...";
  outputStatus();
	
  unsigned int totalNumberOfCutoffs = multimapFromCutoffToJoinedSites.size(); // Number of unique keys, or energy cutoffs. 
  unsigned int currentNumberOfCutoffs = 0;
  
  float previousPercentComplete = 0;
  
  // cycle over each cutoff - rigid cluster pair.
  for( multimap<Cutoff, set<SiteID> >::iterator pairIterator = multimapFromCutoffToJoinedSites.begin();
        pairIterator != multimapFromCutoffToJoinedSites.end();
       pairIterator++ ) {
    /*  for (std::multimap<Cutoff, std::set<SiteID> >::reverse_iterator pairIterator = multimapFromCutoffToJoinedSites.rbegin();
	pairIterator != multimapFromCutoffToJoinedSites.rend();
	pairIterator++) {*/
    
    pair<Cutoff, set<SiteID> > pair = (*pairIterator);
    
    Cutoff cutoff = pair.first;
    
    set<SiteID> setOfJoinedSiteIDs = pair.second;
    if (setOfJoinedSiteIDs.size() == 0) {
      continue;
    }
    
    set<SiteID>::iterator joinedSiteIDIterator = setOfJoinedSiteIDs.begin();
    
    // use the first siteID as an anchor site and join each subsequent site to it
    SiteID anchorSiteID = *joinedSiteIDIterator;
    joinedSiteIDIterator++;
    //cout << "? Cutoff:" << cutoff << endl;
    //cout << "  " << anchorSiteID << endl;
    
      // Cycle over all the atoms in this rigid cluster. 
    //////////////////////////////////////////////////////////////////////
    while (joinedSiteIDIterator != setOfJoinedSiteIDs.end()) {

      SiteID currentSiteID = *joinedSiteIDIterator;
      RigidCluster *possibleParentRigidCluster = NULL;

      //cout << "  " << currentSiteID << endl;
     
      possibleParentRigidCluster = joinSites(cutoff, anchorSiteID, currentSiteID);
      
      
      /*      if (parentRigidCluster != NULL) {
	      if (parentRigidCluster->size() < possibleParentRigidCluster->size() ) {
	      parentRigidCluster = possibleParentRigidCluster;
	      }
	      
	      } else {
	      parentRigidCluster = possibleParentRigidCluster;
	      
	      }*/
      
      joinedSiteIDIterator++;
    }
    
    currentNumberOfCutoffs++;
    if(parameters.flexweb) {
      // compute and update percentComplete for flexweb status 
      float percentComplete = 100*(float) currentNumberOfCutoffs / totalNumberOfCutoffs;
      if (percentComplete - previousPercentComplete > 5) {
	stringstream percentCompleteString;
	percentCompleteString << "" << noshowpoint << fixed << setprecision(0) << percentComplete << "% Complete";
	parameters.mapFromTaskNameToStatus["Decompose into rigid clusters"] = percentCompleteString.str();
	outputStatus();		
	
	previousPercentComplete = percentComplete;
      }
    }
  }
  
  parameters.mapFromTaskNameToStatus["Decompose into rigid clusters"] = "Complete";
  outputStatus();
  
  return parentRigidCluster;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void RigidClusterFactory::insertRigidSubCluster(RigidCluster * rigidCluster,
                                                Cutoff cutoff,
                                                RigidCluster * rigidSubCluster) {
  
  rigidCluster->insertSubCluster(rigidSubCluster);
  
  rigidClusterHierarchy->setOfLargestRigidClusters.insert(rigidCluster);
  rigidClusterHierarchy->setOfLargestRigidClusters.erase(rigidSubCluster);

  rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster[rigidSubCluster]  = rigidCluster;

}

////////////////////////////////////////////////////////////////////////////////
RigidCluster *RigidClusterFactory::joinRigidClusters( Cutoff cutoff,
						      RigidCluster *firstRigidCluster, 
						      RigidCluster *secondRigidCluster ) {

  rigidClusterHierarchy->setOfCutoffs.insert(cutoff);
  
  RigidCluster *joinedRigidCluster = NULL;
  
  if( firstRigidCluster == secondRigidCluster ) {	
    // the rigid clusters we're trying to join are the same - we are already done 
    
    joinedRigidCluster = firstRigidCluster;
    
  } 
  else {
    joinedRigidCluster = new RigidCluster();
    
    insertRigidSubCluster(joinedRigidCluster,
			  cutoff,
			  firstRigidCluster);
    
    insertRigidSubCluster(joinedRigidCluster,
			  cutoff,
			  secondRigidCluster);
    
    // store a copy of the current setOfLargestRigidClusters 
    rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters[cutoff] = rigidClusterHierarchy->setOfLargestRigidClusters;	
    
    
    if( (rigidClusterHierarchy->mapFromRigidClusterToCutoff[joinedRigidCluster] < cutoff) || 
	(rigidClusterHierarchy->mapFromRigidClusterToCutoff[joinedRigidCluster] == 0.0)) {
      
      rigidClusterHierarchy->mapFromRigidClusterToCutoff[joinedRigidCluster] = cutoff;
    }
    
  }
  
  return joinedRigidCluster;
}

////////////////////////////////////////////////////////////////////////////////
RigidCluster *RigidClusterFactory::getLargestRigidClusterContainingRigidCluster(RigidCluster * rigidCluster, 
                                                                                 unsigned int recursionDepth) {

/*  unsigned int maximumRecursionDepth = 10;
  
  if (recursionDepth > maximumRecursionDepth) {
    std::cerr << " exceeded maximumRecursionDepth" << std::endl;
    return rigidCluster;
  }*/
  
	RigidCluster * largestRigidCluster = rigidCluster;
		
	if (rigidClusterHierarchy->mapFromRigidClusterToLargestKnown.find(rigidCluster) !=
      rigidClusterHierarchy->mapFromRigidClusterToLargestKnown.end()) {
		largestRigidCluster = rigidClusterHierarchy->mapFromRigidClusterToLargestKnown[rigidCluster];
	}
  
	if (rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster.find(largestRigidCluster) 
	    != rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster.end()) {
    RigidCluster * parentRigidCluster = rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster[rigidCluster];
    largestRigidCluster = getLargestRigidClusterContainingRigidCluster(parentRigidCluster, recursionDepth+1);
		
	}
  
	rigidClusterHierarchy->mapFromRigidClusterToLargestKnown[rigidCluster] = largestRigidCluster;
  
	return largestRigidCluster;
}

////////////////////////////////////////////////////////////////////////////////
RigidCluster * RigidClusterFactory::getLargestRigidClusterContainingSiteID(SiteID siteID) {	

  RigidCluster * largestRigidCluster = NULL;
	
  if (rigidClusterHierarchy->mapFromSiteIDToRigidClusterContainingSiteID.find(siteID) 
      != rigidClusterHierarchy->mapFromSiteIDToRigidClusterContainingSiteID.end()) {
		
    RigidCluster * rigidCluster = rigidClusterHierarchy->mapFromSiteIDToRigidClusterContainingSiteID[siteID];
    
    largestRigidCluster = getLargestRigidClusterContainingRigidCluster(rigidCluster);
    
  } else {
    // no rigid cluster exists - create one containing only siteID
    RigidCluster * rigidCluster = new RigidCluster(siteID);
    
    rigidClusterHierarchy->mapFromSiteIDToRigidClusterContainingSiteID[siteID] = rigidCluster;
    
    largestRigidCluster = rigidCluster;
  }
  
  return largestRigidCluster;
}
