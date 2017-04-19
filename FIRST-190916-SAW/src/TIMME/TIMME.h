#ifndef _TIMME_H_
#define _TIMME_H_

#include <vector>
#include <map> 
#include <set>
#include <string>
#include <algorithm>
#include "math.h"
#include "global_defs.h"
#include "generalUtils.h"
#include "MolFramework.h" 
#include "flexweb.h"
#include "Site_Info.h"
#include "VectorAlgebra.h"
#include "readPDB.h"
#include "RigidCluster.h"
#include "RigidClusterFactory.h"

////////////////////////////////////////////////////////////////////////////////
// Description:
//   TIMME.H - definition of the Tool for Identification of Mobility in Macromolecular Ensembles class
//
//   Scott Menor
//   Arizona State University
//   Department of Physics and Astronomy
//   Biophysics Theory Group
//   PO Box 871504
//   Tempe, AZ 85287
//   Scott.Menor@asu.edu
//
//   Copyright (c) 2005-2007, Arizona State University, All Rights Reserved.
////////////////////////////////////////////////////////////////////////////////

class TIMME
{
public:
	TIMME();
	
	void addStructure(MolFramework);
	
	// methods for ensemble analysis
	////////////////////////////////////////////////////////////////////////////////
	void validateStructures();
	void rigidClusterDecomposition();
	void rigidClusterDilution();
		
	void calculateMobility();
	void calculateMobility(unsigned int siteID);
	void calculateFlexibility();
	float calculateFlexibility(unsigned int siteID);
	
	bool allStructuresHaveEqualNumberOfSites();
	
	void saveOutput();
	void saveOutput(std::vector<MolFramework*> *vectorOfStructures, 
                  string outputFilename);
  
	void savePDBmodel(MolFramework* model,
                    std::map<unsigned int, float> &mapFromSiteIDtoOccupancy,
                    std::map<unsigned int, float> &mapFromSiteIDtoCharge,
                    ofstream &outputFileStream);
	
	void saveResidueData (string filename,
                        MolFramework* structure, 
                        std::map<unsigned int, float> &mapFromSiteIDtoFloat);
	
	void saveResidueRMSflexibility(string filename);
	
	float calculatePairStandardDeviationSquared(SiteID , SiteID);
  
	float calculatePairStandardDeviationSquared(std::set<SiteID> firstSiteIDs , std::set<SiteID> secondSiteIDs);
  
  //	float calculatePairStandardDeviationSquaredPerMeanSquared(unsigned int , unsigned int);
	float calculateStandardDeviationBetweenMeanCoordinatesAndSite(unsigned int siteIndex);
	float calculateStandardDeviationSquaredBetweenMeanCoordinatesAndSite(unsigned int siteIndex);
	
	float calculateResidueRMSflexibility(unsigned int residueID); 
	
	std::set<unsigned int> * getSetOfSiteIDsFromResidue(int residueSeqID);
	
	// perhaps this should be refactored into the molecular framework class? (or a generic utility class)
	std::set<unsigned int> getNearestNeighbors(std::set<unsigned int> siteIDs);  
	
	std::set<unsigned int> getNthNearestNeighbors(std::set<unsigned int> siteIDs, 
                                                int n);
	
	std::set<unsigned int> getNthNearestNeighbors(unsigned int siteID, 
                                                int n);
  
	std::set<unsigned int> getUpToNthNearestNeighbors(std::set<unsigned int> siteIDs, 
                                                    int n);
	
	std::set<unsigned int> getUpToNthNearestNeighbors(unsigned int siteID, 
                                                    int n);
  
	std::map<unsigned int, std::set<unsigned int> > mapFromResidueSeqIDtoSetOfSites;
	
	bool isSiteIncludedInAnalysis(unsigned int);
	
	void computeMeanCoordinates();
  
  bool considerPairOfSites(SiteID, SiteID);
	
	virtual ~TIMME();
	
private:
  float * period;
    MolFramework *referenceStructure;
	
	std::map<MolFramework*, float*> mapFromStructureToMeanCoordinates;
	
	std::vector<MolFramework*> *vectorOfStructures;
	
	// FIXME - DRY - remove these and move them to the authoratative referenceStructure
	std::map<SiteID, float> mapFromSiteIDtoFlexibility;	
	std::map<SiteID, float> mapFromSiteIDtoMobility;		
	
	std::map<int, std::map<string, SiteID> > mapFromSeqNumberToMapFromSiteNameToSite;
};

#endif //_TIMME_H_
