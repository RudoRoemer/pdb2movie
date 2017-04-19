#include <algorithm>
#include <map>

#include "RigidClusterAnalysis.h"
#include "TIMME.h"

#include "Assembly.h"
#include "TIMME_Assembly.h"
#include "OutputGraphML.h"
#include "TIMME_AssemblyOutputSVG.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
//
//   $Id: TIMME.cpp,v 1.1 2011/04/23 04:38:26 dwfarrel Exp $ 
//   definition of the Tool for Identification of Mobility in Macromolecular Ensembles class
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

// Preferred constructor.
////////////////////////////////////////////////////////////////////////////////
TIMME::TIMME() {
		
	vectorOfStructures = new std::vector<MolFramework*>;
	readPDB_Data((*vectorOfStructures), 
							 parameters.infile_name);
	
	if (vectorOfStructures->size() < 2) {
		std::cout << std::endl << "error: TIMME analysis requires two or more structures ";
		std::cout << "( only " << vectorOfStructures->size() << " was loaded)." << std::endl;
		
		exit(1);
	}
	
	referenceStructure = vectorOfStructures->at(0);
  
	for (std::vector<MolFramework*>::iterator structureIterator = vectorOfStructures->begin();
			 structureIterator != vectorOfStructures->end();
			 structureIterator ++) {
	
		MolFramework * structure = *structureIterator;
		
	  structure->includeHydrogenbonds = parameters.includeHydrogenbonds;
		
	  structure->includeHydrophobicTethers = parameters.includeHydrophobicTethers;
		
	  excludeSites(*structure);
	  setVdwRadii(*structure);
	  buildFramework(*structure); // TODO - figure out the easiest way to exclude hydrogen bonds and hydrophobic tethers
	}
	
	if( parameters.verbose ) {
		std::cout << vectorOfStructures->size() << " structures loaded" << std::endl;
	}
}

//
// Compute emprical rigid cluster decomposition hierarchy using a hard O(n^2) algorithm.
void TIMME::rigidClusterDecomposition() {

  TIMME_Assembly testDilution;
  //Assembly<Cutoff, SiteID> testDilution;

	std::cout << " Computing Rigid Cluster Hierarchy..." << std::endl;
	parameters.mapFromTaskNameToStatus["Compute rigid cluster hierarchy"] = "Running...";
	outputStatus();
	
	validateStructures();  // ensure that preconditions on structures are satisfied 
	
	unsigned int totalNumberOfPairs = referenceStructure->total_sites*(referenceStructure->total_sites-1);
	
	parameters.mapFromTaskNameToStatus["Compute pair deviations"] = "Running...";
	outputStatus();
	
	std::multimap<Cutoff, SiteID*> multimapFromCutoffToJoinedSites;
  
  ofstream * sigmaMatrixOutputFile = NULL;
  
  if (parameters.saveSigma) {
    // save the sparse or dense sigma matrix in the Matrix Market format
    // http://math.nist.gov/MatrixMarket/formats.html#MMformat
    string sigmaMatrixOutputFilename = 	referenceStructure->base_name + "_sigma.mtx";
    
    sigmaMatrixOutputFile = new ofstream(sigmaMatrixOutputFilename.c_str());
    
  }
  
  int totalNonzeroEntries = 0; // FIXME - precompute this to satisfy the mm format
  if (sigmaMatrixOutputFile != NULL) {
    *sigmaMatrixOutputFile << "%%MatrixMarket matrix coordinate real symmetric" << std::endl;
    *sigmaMatrixOutputFile << "%" << std::endl;
    
    *sigmaMatrixOutputFile << referenceStructure->total_sites << "\t" << referenceStructure->total_sites << "\t" << totalNonzeroEntries << std::endl;
  }
	
	ofstream *neighborListFile = NULL;
	
	if (parameters.outputNeighborList) {
		string neighborListFilename = referenceStructure->base_name + "_dilution.neighbors";
		neighborListFile = new ofstream(neighborListFilename.c_str());  
	}
  
	cout << "Computing pairwise deviations" << endl;	
	
	std::set<SiteID> nthNearestNeighbors;
	
	ofstream *bracingFile = NULL;
	if (true) { // FIXME - add a new parameter / command line option for this
		string bracingFilename = referenceStructure->base_name + "_bracing.py";
		bracingFile = new ofstream(bracingFilename.c_str());
		*bracingFile << "def showBraces(thresholdDeviation) :" << std::endl;
	}
	
	// PRECONDITION WARNING - assumes that all structures have the same number of sites and that every site is present in every model
	for (SiteID firstSiteID = 1; firstSiteID <= referenceStructure->total_sites; firstSiteID++) {			 


	  if( fmod( (double) firstSiteID, 100.0) == 0.0 )
	    cout << "Checking site:" << firstSiteID << endl;
		
		float plasticity = 0.0;
		
		if (parameters.outputNeighborList) {
			*neighborListFile << firstSiteID << "\t";
		}
		
		if (parameters.covalentTIMME) {
			nthNearestNeighbors = getUpToNthNearestNeighbors(firstSiteID, 
																											 parameters.timme_nThNeighborForPlasticity);
		}
		
//#pragma omp parallel for shared(multimapFromCutoffToJoinedSites)
		for (SiteID secondSiteID = 1; secondSiteID < firstSiteID; secondSiteID++) {
			bool considerPair = false;
			
			if (parameters.useLinearTimme) {
				for (std::vector<MolFramework*>::iterator structureIterator = vectorOfStructures->begin();
						 structureIterator != vectorOfStructures->end();
						 structureIterator++) {
					
					MolFramework * structure = *structureIterator;
					
					// compute pair distance in this structure
					float pairDistance = structure->getDistance(firstSiteID, 
																											secondSiteID);
					
					// if pair distance in structure is less than the threshold 
					// then consider the pair, and we are done
					if (pairDistance < parameters.timme_neighborCutoffDistance) {
						considerPair = true;
						break;
					}
				}				
			}
			
			else { // (!parameters.useLinearTimme) 
				// consider every possible pair in this case (exhaustive TIMME)
				considerPair = true;
			}
			
			// if we are analyzing a covalently bonded network, 
			// ignore up to n-th nearest neighbors 
			if (parameters.covalentTIMME) {
				if (nthNearestNeighbors.find(secondSiteID) != nthNearestNeighbors.end()) {
					considerPair = false;
				} 
			}
			
			if (considerPair) {
				
				float currentDeviationSquared = 0.0f;
				
				if (parameters.timme_bodyClusterSize > 0) {
					std::set<SiteID> firstBody  = getUpToNthNearestNeighbors(firstSiteID, parameters.timme_bodyClusterSize);
					std::set<SiteID> secondBody = getUpToNthNearestNeighbors(secondSiteID, parameters.timme_bodyClusterSize);
					
					currentDeviationSquared = calculatePairStandardDeviationSquared(firstBody, secondBody);
					
				} else {
					
					currentDeviationSquared = calculatePairStandardDeviationSquared(firstSiteID, secondSiteID);
				}
				
				float currentDeviation = sqrt(currentDeviationSquared);					

				if (currentDeviation > 0.0) {
					plasticity += 1.0/currentDeviation;				 
				}
				
				if (parameters.timme_complement) {
					currentDeviation *= -1.0;
				}
				
				{ // FIXME / TODO - replace this with a better way of selecting cutoffs
					currentDeviation = floor(currentDeviation*100.)/100.;
				}
				
				if (sigmaMatrixOutputFile != NULL) {
					*sigmaMatrixOutputFile << firstSiteID << "\t" << secondSiteID << "\t" << currentDeviation << std::endl;
				}

				SiteID *newPair = new SiteID[2];
				newPair[0] = firstSiteID;
				newPair[1] = secondSiteID;
				
				if (bracingFile != NULL) {
					if (currentDeviation < .15) {
						*bracingFile << "\tif (" << currentDeviation << " < thresholdDeviation) :" << std::endl;
						*bracingFile << "\t\tcmd.dist('sigmaBraces', 'id " << firstSiteID << "', 'id " <<  secondSiteID << "')" << std::endl;						
					}
				}
				
				// a pair of sites is added at each step. 
				if( parameters.timme_useAssembly ){
				  testDilution.insertPiece( currentDeviation, newPair );
				  delete newPair;
				}
				else{
				  multimapFromCutoffToJoinedSites.insert(std::pair<Cutoff, SiteID*> (currentDeviation, newPair));
				}
			}
		}
		
		mapFromSiteIDtoFlexibility[firstSiteID] = .01*plasticity;
				
		unsigned int currentNumberOfPairs = firstSiteID*(firstSiteID-1);
				
		float percentComplete = 100*(float)currentNumberOfPairs/totalNumberOfPairs;
		stringstream percentCompleteString;
		percentCompleteString << "" << noshowpoint << fixed << setprecision(0) << percentComplete << "% Complete";
		parameters.mapFromTaskNameToStatus["Compute pair deviations"] = percentCompleteString.str();
		outputStatus();
	}
	 
	if (sigmaMatrixOutputFile != NULL) {
		sigmaMatrixOutputFile->close();
	}
	
	if (bracingFile != NULL) {
		bracingFile->close();
	}

	if (parameters.outputNeighborList) {
		neighborListFile->close() ;
	}
	 
	 parameters.mapFromTaskNameToStatus["Compute pair deviations"] = "Complete";
	 outputStatus();

	 RigidClusterFactory rigidClusterFactory;
	 rigidClusterFactory.decomposeIntoRigidClusters(multimapFromCutoffToJoinedSites);

	 RigidClusterHierarchy * rigidClusterHierarchy = rigidClusterFactory.getRigidClusterHierarchy();
	 
	 // Assembly-specific options. 
	 //////////////////////////////////////////////////////////////////////
	 if( parameters.timme_useAssembly ){
	   testDilution.acceptVisitor( new OutputGraphML<Cutoff, SiteID> ( referenceStructure) );

	   TIMME_AssemblyOutputSVG<Cutoff, SiteID> *timme_assemblyOutputSVG = new TIMME_AssemblyOutputSVG<Cutoff, SiteID>( referenceStructure );
	   testDilution.acceptVisitor( timme_assemblyOutputSVG );
	 }
	 
	 std::cout << "   rigidClusterAnalysis" << std::endl;
	 RigidClusterAnalysis rigidClusterAnalysis(rigidClusterHierarchy); 
 	 rigidClusterAnalysis.setStructure(referenceStructure);
	 rigidClusterAnalysis.setXAxisType(parameters.dilutionPlotXAxisType);
	 
	 std::cout << "   assignRigidClusterLabels" << std::endl;
	 rigidClusterAnalysis.assignRigidClusterLabels(); // FIXME - is there a way to do this while joining sites?
	 
	 std::cout << "   computeRigidClusterColoring" << std::endl;
	 rigidClusterAnalysis.computeRigidClusterColoring();
	 
	 string dilutionFilename = referenceStructure->base_name + "_dilution";
	 std::cout << "   done" << std::endl;
	 
	 std::cout << " Saving Dilution..." << std::endl;

	 if ( !parameters.timme_useAssembly ){
		 rigidClusterAnalysis.saveDilution(dilutionFilename);
	 }
	 
	 std::cout << "   done" << std::endl;
	 
	 //calculateFlexibility(); 
	 //calculateMobility();
	 
	 // output the hierarchy as a PDB file
	 saveOutput();
	 
	 parameters.mapFromTaskNameToStatus["Compute rigid cluster hierarchy"] = "Complete";
	 outputStatus();
}

void TIMME::computeMeanCoordinates() {
	
	for (std::vector<MolFramework*>::iterator structureIterator = vectorOfStructures->begin();
			 structureIterator != vectorOfStructures->end();
			 structureIterator ++) { // TODO - foreach 
		
		float *sumOfCoordinates = new float[3];
		for (int i=0;i<3;i++) {
			sumOfCoordinates[i] = 0.0; // TODO - do we need to initialize this?
		}
		
		MolFramework *framework = (*structureIterator);
		
		for (unsigned int siteID = 0;
				 siteID < framework->total_sites;
				 siteID ++) {
			
			// TODO - get the coordinates for the site and sum over them 
			Site_Info site = framework->site_info[siteID]; // TODO - getSiteInfo
			
			for (int i=0;i<3;i++) {
				sumOfCoordinates[i] += site.coords[i];
			}
		}
		
		float *meanCoordinates = new float[3];
		for (int i=0;i<3;i++) {
			meanCoordinates[i] = sumOfCoordinates[i] / framework->total_sites;
		}
		
		mapFromStructureToMeanCoordinates[framework] = meanCoordinates;
	}
}

void TIMME::validateStructures() {
	
	if (vectorOfStructures->size() <= 1) {
		std::cerr << vectorOfStructures->size() + " structures were loaded - TIMME analysis requires at least 2, could not proceed" << std::endl;
		// TODO - use new logging system
		exit(-1);
	}
	
	if (!allStructuresHaveEqualNumberOfSites()) {
		std::cerr << "all structures do not have an equal number of sites - TIMME analysis could not proceed" << std::endl;
		// TODO - use new logging system
		exit(-1); // TODO - fail more gracefully ?
	}	
}

// borrowed from Brandon's output_RCD_PDB_format method in MolFramework_Output 
// TODO - refactor to MolFramework_Output
void TIMME::savePDBmodel(MolFramework* model,
                         std::map<unsigned int, float> &mapFromSiteIDtoOccupancy,
                         std::map<unsigned int, float> &mapFromSiteIDtoCharge,
                         ofstream &outputFileStream) {
	for(unsigned int siteID = 1; siteID <= model->total_sites; siteID++ ){
		// TODO - ensure that mapFromSiteIDtoOccupancy and mapFromSiteIDtoBfactor don't map to values outside of tmapFromSiteIDtoBfactor range (-999 to 999)
		
		if( siteID <= 99999 ){		
			Site_Info siteInfo = model->site_info[siteID];
			outputFileStream  << showpoint << setiosflags(ios::fixed) 
				<< setw(6) << siteInfo.record_name 
				<< setw(5) << siteInfo.orig_atom_number
				<< setw(5) << model->atomNamePDBFormat(siteInfo.atom_name, siteInfo.element_name)
				<< setw(4) << siteInfo.residue_name
				<< setw(2) << char(siteInfo.chain_ID)
				<< setw(4) << siteInfo.seq_number
				<< "    "
				<< setw(8) << setprecision(3) << siteInfo.coords[X]
				<< setw(8) << setprecision(3) << siteInfo.coords[Y]
				<< setw(8) << setprecision(3) << siteInfo.coords[Z]
				<< setw(6) << setprecision(2) << mapFromSiteIDtoOccupancy[siteID]
				<< setw(6) << setprecision(2) << (((float)referenceStructure->site_info[siteID].rigid_label)/100.) // FIXME - magic number and a bit ugly / cryptic
				<< "            "
				//			<< setw(4) << "      " // TODO - segID
				//			<< setw(2) << "  " // TODO - something like siteInfo.element
				//			<< setw(2) << setprecision(1) << mapFromSiteIDtoCharge[siteID] // TODO - only integers!
				<< endl; 
		}
	}
}

// save map data for each site number into a two-column, tab-delimited text file 
//
void saveSiteData (string filename,
                   MolFramework* structure, 
                   std::map<unsigned int, float> &mapFromSiteIDtoFloat) {
	
	ofstream outputFileStream(filename.c_str());
	
	for (std::map<unsigned int, float>::iterator siteIDtoFloatPairIterator = mapFromSiteIDtoFloat.begin();
			 siteIDtoFloatPairIterator != mapFromSiteIDtoFloat.end();
			 siteIDtoFloatPairIterator++) {
		
		unsigned int siteID = siteIDtoFloatPairIterator->first;
		float value = siteIDtoFloatPairIterator->second;
		
		outputFileStream << siteID << "\t" << value << std::endl; // TODO - we probably really want residue IDs here
	}
	
	outputFileStream.close();
}

// save map data for each residue sequence number into a three-column, tab-delimited text file 
// 
void TIMME::saveResidueData (string filename,
                             MolFramework* structure, 
                             std::map<unsigned int, float> &mapFromSiteIDtoFloat) {
	
  ofstream outputFileStream(filename.c_str());
	
  for (std::set<int>::iterator sequenceNumberIterator = structure->seq_numbers.begin();
       sequenceNumberIterator != structure->seq_numbers.end();
       sequenceNumberIterator++) {
		
    int residueSeqID = *sequenceNumberIterator;
		
    std::set<unsigned int> * setOfSiteIDs = getSetOfSiteIDsFromResidue(residueSeqID);
		
    for (std::set<unsigned int>::iterator siteIDiterator = setOfSiteIDs->begin();
				 siteIDiterator != setOfSiteIDs->end();
				 siteIDiterator++) {
			
      unsigned int siteID = *siteIDiterator;
			
      float floatValue = mapFromSiteIDtoFloat[siteID];
			
      outputFileStream << residueSeqID << "\t" << referenceStructure->site_info[siteID].atom_name << "\t" << floatValue << std::endl;	
    }			
  }
	
  outputFileStream.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Output the inferred rigid custer decomposition, flexibility and mobility 
////////////////////////////////////////////////////////////////////////////////
void TIMME::saveOutput(std::vector<MolFramework*> *vectorOfStructures, 
                       string outputPrefix) {
	string pdbOutputFilename = outputPrefix + ".pdb";
	
	ofstream outputPDBfile(pdbOutputFilename.c_str());
	
	unsigned int modelNumber = 1;
	for (std::vector<MolFramework*>::iterator structureIterator = vectorOfStructures->begin();
			 structureIterator != vectorOfStructures->end();
			 structureIterator ++) {
		
		MolFramework* model = (*structureIterator);
		
		outputPDBfile  << noshowpoint << setiosflags(ios::fixed);
		outputPDBfile << "MODEL     ";
		outputPDBfile << setw(4) << modelNumber << std::endl;
		
		// TODO - mapFromSiteIDtoMobility
		// TODO - add user options to choose where to store mobility, flexibility and rigidity
		savePDBmodel(model,
								 mapFromSiteIDtoFlexibility,
								 mapFromSiteIDtoMobility,
								 outputPDBfile);
		
		outputPDBfile << "ENDMDL" << endl;
		
		modelNumber ++;
	}
	
	outputPDBfile.close();
	
	saveResidueRMSflexibility(referenceStructure->base_name+".residueRMSflexibility");
	
	saveResidueData(referenceStructure->base_name+".residueFlexibility",
									referenceStructure,
									mapFromSiteIDtoFlexibility);
	
	saveResidueData(referenceStructure->base_name+".residueMobility",
									referenceStructure,
									mapFromSiteIDtoMobility);	
	
	output_FRODA_TIMME_Jmol_script(*referenceStructure); // TODO - removeme?
	
}

void TIMME::saveOutput() {
	//	string outputPrefix = referenceStructure->base_name + "_timme";
	string outputPrefix = referenceStructure->base_name + "_dilution";  // TODO - should we save this as _RCD?
	saveOutput(vectorOfStructures,
						 outputPrefix);
}

// 
//
//////////////////////////////////////////////////////////////////////
float TIMME::calculateStandardDeviationSquaredBetweenMeanCoordinatesAndSite(unsigned int siteIndex) {
	// calculate mean pair distance
	float sumOfPairDistances = 0.;
	
	for (std::vector<MolFramework*>::iterator structureIterator = vectorOfStructures->begin();
			 structureIterator != vectorOfStructures->end();
			 structureIterator ++) {
		
		MolFramework* structure = (*structureIterator);
		float *referenceCoordinates = mapFromStructureToMeanCoordinates[structure];
		
		Site_Info siteInfo = structure->site_info[siteIndex];
		float *siteCoordinates = siteInfo.coords;
		
		float pairDistance = VectorAlgebra::distance(referenceCoordinates, siteCoordinates);
		sumOfPairDistances += pairDistance;
	}
	
	int numberOfStructures = vectorOfStructures->size();
	float meanPairDistance = sumOfPairDistances/((float)numberOfStructures);
	
	// calculate sum of squared deviations from the mean pair distance
	float sumOfSquaredPairDeviations = 0.;
	
	for (std::vector<MolFramework*>::iterator structureIterator = vectorOfStructures->begin();
			 structureIterator != vectorOfStructures->end();
			 structureIterator ++) {
		
		MolFramework* structure = (*structureIterator);
		float *referenceCoordinates = mapFromStructureToMeanCoordinates[structure];
		
		Site_Info siteInfo = structure->site_info[siteIndex];
		float *siteCoordinates = siteInfo.coords;
		
		float pairDistance = VectorAlgebra::distanceSquared(referenceCoordinates, siteCoordinates); // FIXME - hotspot
		float squaredDifferenceFromMean = pow(pairDistance - meanPairDistance, 2);           // FIXME - hotspot
		sumOfSquaredPairDeviations += squaredDifferenceFromMean;
	}
	
	float standardDeviationSquared = sumOfSquaredPairDeviations/((float)numberOfStructures);
	return standardDeviationSquared;
}

float TIMME::calculateStandardDeviationBetweenMeanCoordinatesAndSite(unsigned int siteIndex) {	
	float standardDeviationSquared = calculateStandardDeviationSquaredBetweenMeanCoordinatesAndSite(siteIndex);
	
	float standardDeviation = sqrt(standardDeviationSquared);
	
	return standardDeviation;
}

// @depricated 
bool TIMME::considerPairOfSites(SiteID firstSiteID, 
                                SiteID secondSiteID) {
  
  // FIXME - make this a command line option or (better) automate it
  float maximumSeparationDistance = parameters.timme_neighborCutoffDistance;
  
  if (referenceStructure->getDistance(firstSiteID, secondSiteID) < maximumSeparationDistance) { 
    return true;
  }
  
  return false;
}

float TIMME::calculatePairStandardDeviationSquared(std::set<SiteID> firstSiteIDs,
                                                   std::set<SiteID> secondSiteIDs) {
  
  float sumOfSquaredDeviations = 0.0;
	
  for (std::set<SiteID>::iterator firstSiteID = firstSiteIDs.begin();
       firstSiteID != firstSiteIDs.end();
       firstSiteID++) {
    
    for (std::set<SiteID>::iterator secondSiteID = secondSiteIDs.begin();
         secondSiteID != secondSiteIDs.end();
         secondSiteID++) {
      sumOfSquaredDeviations += calculatePairStandardDeviationSquared(*firstSiteID, *secondSiteID);
			
    }
  }
	
  float standardDeviationSquared = sumOfSquaredDeviations / ((float)(firstSiteIDs.size() * secondSiteIDs.size()));
  
  return standardDeviationSquared;
}

// calculate the standard deviation of pair distances between two sites 
// 
//////////////////////////////////////////////////////////////////////
float TIMME::calculatePairStandardDeviationSquared(SiteID firstSiteIndex, 
                                                   SiteID secondSiteIndex) {
	
  bool cached = true;
  
  std::pair<SiteID, SiteID> parameters(firstSiteIndex, secondSiteIndex);
	
  static std::map < std::pair<SiteID, SiteID>, float> cache;
  
  static unsigned int calls = 0;
  static unsigned int hits = 0;
  calls ++;
  
  if (cached) {
    if (cache.find(parameters) != cache.end()) {
      hits ++;
			
      return cache[parameters];
    }
  }
  
	if (firstSiteIndex == secondSiteIndex) {
		return 0.0f;
	}
	
	int numberOfStructures = vectorOfStructures->size();
  
	float sumOfPairDistances = 0.0f;
	
	for (int structureNumber=0;structureNumber<numberOfStructures;structureNumber++) {
		float pairDistance = (*vectorOfStructures)[structureNumber]->getDistance(firstSiteIndex, secondSiteIndex);	
		
		sumOfPairDistances += pairDistance;
	}
	
	float meanPairDistance = sumOfPairDistances / (float)numberOfStructures;
	
	float sumOfSquaredDeviations = 0.0f;
	for (int structureNumber=0;structureNumber<numberOfStructures;structureNumber++) {
		float pairDistance = (*vectorOfStructures)[structureNumber]->getDistance(firstSiteIndex, secondSiteIndex);	
		
		sumOfSquaredDeviations += pow(pairDistance - meanPairDistance, 2);
	}
	
	float standardDeviationSquared = sumOfSquaredDeviations / (float)numberOfStructures;
	
	//	standardDeviationSquared /= pow(meanPairDistance, 2); // FIXME - removeme
	
  if (cached) {
    cache[parameters] = standardDeviationSquared;
  }
  
	return standardDeviationSquared;
}

bool TIMME::allStructuresHaveEqualNumberOfSites() {
	int numberOfStructures = vectorOfStructures->size();
	if (numberOfStructures>0) {
		MolFramework *firstStructure = vectorOfStructures->at(0);
		unsigned int numberOfSites = firstStructure->total_sites;
		
		for (int structureNumber=1;structureNumber<numberOfStructures;structureNumber++) {
			
			MolFramework *currentStructure = vectorOfStructures->at(structureNumber);
			
			if (currentStructure->total_sites != numberOfSites) {
				return false; // found a structure with a different number of sites :-( 
			}
		}
	}
	
	return true;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculate mobility (RMSD in Angstroms) of each site, relative to the rigid core  
////////////////////////////////////////////////////////////////////////////////
void TIMME::calculateMobility() {
	std::cout << " Computing Mobility..." << std::endl;
	parameters.mapFromTaskNameToStatus["Calculate mobility"] = "Running...";
	outputStatus();
	
	computeMeanCoordinates(); // TODO - compute center of mass (small but significant difference)
	
	// TODO - implement
	// iterate through all sites and compute the RMSD from the mean coordinates for each site
	for (unsigned int siteID = 0; 
			 siteID < referenceStructure->total_sites; 
			 siteID ++) {
		
		float siteMobility = calculateStandardDeviationBetweenMeanCoordinatesAndSite(siteID);
		mapFromSiteIDtoMobility[siteID] = siteMobility;
	}
	
	parameters.mapFromTaskNameToStatus["Calculate mobility"] = "Complete";
	outputStatus();
}

std::set<SiteID> TIMME::getUpToNthNearestNeighbors(SiteID siteID, 
                                                   int n) {
	std::set<SiteID> siteIDs;
	siteIDs.insert(siteID);
	
	return getUpToNthNearestNeighbors(siteIDs, 
                                    n);
}

std::set<SiteID> TIMME::getUpToNthNearestNeighbors(std::set<SiteID> siteIDs, 
                                                   int n) {
	
	if (n>0) {
		std::set<SiteID> nearestNeighbors = getNearestNeighbors(siteIDs);
		std::set<SiteID> upToNthNearestNeighbors = getUpToNthNearestNeighbors(nearestNeighbors, 
																																					n-1);
		upToNthNearestNeighbors.insert(siteIDs.begin(), siteIDs.end());
		
		return upToNthNearestNeighbors;
		
	} else {
		if (n==0) {
			return siteIDs; 
		} 
		
		// n < 0 - return empty set
		std::set<SiteID> emptySet; // FIXME - is there an STL empty set? 
		return emptySet;
	}
}

std::set<SiteID> TIMME::getNthNearestNeighbors(SiteID siteID, 
                                               int n) {
	std::set<SiteID> siteIDs;
	siteIDs.insert(siteID);
	
	return getNthNearestNeighbors(siteIDs, 
																n);
}

std::set<SiteID> TIMME::getNthNearestNeighbors(std::set<SiteID> siteIDs, 
                                               int n) {
	
	if (n>0) {
		std::set<SiteID> nearestNeighbors = getNearestNeighbors(siteIDs);
		std::set<SiteID> nThNearestNeighbors = getNthNearestNeighbors(nearestNeighbors, 
                                                                  n-1);
		// remove all siteIDs from nThNearestNeighbors (TODO - refactor this out to its own method)
		for (std::set<SiteID>::iterator siteIDiterator = siteIDs.begin();
				 siteIDiterator != siteIDs.end();
				 siteIDiterator++) {
			
			unsigned int excludedSiteID = (*siteIDiterator);
			
			if (nThNearestNeighbors.find(excludedSiteID) != nThNearestNeighbors.end()) {
				nThNearestNeighbors.erase(excludedSiteID);				
			}
		}
		
		return nThNearestNeighbors;
		
	} else {
		if (n==0) {
			return siteIDs; 
		} 
		
		// n < 0 - return empty set
		std::set<unsigned int> emptySet; // FIXME - is there an STL empty set? 
		return emptySet;
	}
}

// TODO - include only covalent bonds (h-bonds or tethers will artificially inflate flexibility)
std::set<SiteID> TIMME::getNearestNeighbors(std::set<SiteID> siteIDs) {
	std::set<SiteID> nearestNeighbors;
	
	for (std::set<SiteID>::iterator siteIDiterator = siteIDs.begin();
			 siteIDiterator != siteIDs.end();
			 siteIDiterator ++) {
		
		unsigned int siteID = *siteIDiterator;
		Site_Info siteInfo = referenceStructure->site_info[siteID];
		
		for (SiteID neighborNumber=0;
				 neighborNumber < siteInfo.neighbor_list.size();
				 neighborNumber++) {
			
			SiteID neighborSiteID =  siteInfo.neighbor_list[neighborNumber];
						
			if (siteIDs.find(neighborSiteID) == siteIDs.end()) {
				nearestNeighbors.insert(neighborSiteID);				 
			}
		}		
	}
	
	return nearestNeighbors;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculate flexibility (sigma in Angstroms) of each site, relative to adjacent sites not constrained  
////////////////////////////////////////////////////////////////////////////////
float TIMME::calculateFlexibility(SiteID siteID) {
	
	std::set<SiteID> nextNearestNeighbors = getNthNearestNeighbors(siteID, 
																																 parameters.timme_nThNeighborForFlexibility-1); 
	
	float sumOfSquaredPairStandardDeviations = 0;
	int numberOfNextNearestNeighbors = nextNearestNeighbors.size();
	
	for (std::set<unsigned int>::iterator firstNextNearestNeighborIterator = nextNearestNeighbors.begin();
			 firstNextNearestNeighborIterator != nextNearestNeighbors.end();
			 firstNextNearestNeighborIterator ++) {
		
		unsigned int firstNextNearestNeighborSiteID = (*firstNextNearestNeighborIterator);
		
		for (std::set<unsigned int>::iterator secondNextNearestNeighborIterator = firstNextNearestNeighborIterator; // TODO - don't include a pair of two identical nextNearestNeighbors in the analysis (peek ahead in the iterator?)
				 secondNextNearestNeighborIterator != nextNearestNeighbors.end();
				 secondNextNearestNeighborIterator ++) {
			
			unsigned int secondNextNearestNeighborSiteID = (*secondNextNearestNeighborIterator);
			
			if (firstNextNearestNeighborSiteID != secondNextNearestNeighborSiteID) { // TODO - fix iterators so this test isn't necessary
				sumOfSquaredPairStandardDeviations += calculatePairStandardDeviationSquared(firstNextNearestNeighborSiteID,
																																										secondNextNearestNeighborSiteID);
			}
		}
	}
	
	int numberOfUniquePairsOfNextNearestNeighbors = numberOfNextNearestNeighbors * (numberOfNextNearestNeighbors-1);
	float meanSumOfSquaredPairStandardDeviations = sumOfSquaredPairStandardDeviations / (float)numberOfUniquePairsOfNextNearestNeighbors;
	
	float rootMeanSquaredPairStandardDeviations = 0; // TODO - define a value for undefined flexibility (sites without enough next-next nearest neighbors; I like -1 but it makes other things messy and complicated)
	if (numberOfUniquePairsOfNextNearestNeighbors > 0) {
		rootMeanSquaredPairStandardDeviations = sqrt(meanSumOfSquaredPairStandardDeviations);
	}
	
	// TODO - assign the flexibility index to the referenceStructure
	referenceStructure->site_info[siteID].flexibility = rootMeanSquaredPairStandardDeviations;// TODO - remove?
		mapFromSiteIDtoFlexibility[siteID] = rootMeanSquaredPairStandardDeviations;
		
		return rootMeanSquaredPairStandardDeviations;
}

// save a two column formated tab-delimited text file containing residue number (first column) and residue root mean squared flexibility (second column)
// 
//
void TIMME::saveResidueRMSflexibility(string filename) {
	ofstream outputFileStream(filename.c_str());
		
	for (std::set<int>::iterator sequenceNumberIterator = referenceStructure->seq_numbers.begin();
			 sequenceNumberIterator != referenceStructure->seq_numbers.end();
			 sequenceNumberIterator++) {
		
		int residueSeqID = *sequenceNumberIterator;
		
		float residueRMSflexibility = calculateResidueRMSflexibility(residueSeqID);
		
		outputFileStream << residueSeqID << "\t" << residueRMSflexibility << std::endl;		
	}
	
	outputFileStream.close();
}

std::set<unsigned int> *TIMME::getSetOfSiteIDsFromResidue(int residueSeqID) {
	if (mapFromResidueSeqIDtoSetOfSites.size() == 0) {
		// first use of mapFromResidueSeqIDtoSetOfSites -> fill it up with data
		// TODO - only update for the specific residue of interest?
		for (unsigned int siteID = 1; siteID <= referenceStructure->total_sites; siteID++) {
			if (referenceStructure->isMainchain(siteID) == 1) { // TODO - brittle hard-coded constant
				
				unsigned int residueSeqNumber = referenceStructure->site_info[siteID].seq_number;
				mapFromResidueSeqIDtoSetOfSites[residueSeqNumber].insert(siteID);
			}
		}
	}
	
	return &mapFromResidueSeqIDtoSetOfSites[residueSeqID];
}

float TIMME::calculateResidueRMSflexibility(unsigned int residueSeqID) {
	// problem - need an easy way to get specific sites from a residue 
	
	// actual calculations
	float residueSSflexibility = 0;
	std::set<unsigned int> *setOfSiteIDs = getSetOfSiteIDsFromResidue(residueSeqID);
	
	for (std::set<unsigned int>::iterator siteIDiterator = setOfSiteIDs->begin();
			 siteIDiterator != setOfSiteIDs->end();
			 siteIDiterator++) {
		
		unsigned int siteID = *siteIDiterator;
		
		residueSSflexibility += calculateFlexibility(siteID);
	}
	
	float residueMSflexibility = residueSSflexibility / (float)setOfSiteIDs->size();
	
	float residueRMSflexibility = sqrt(residueMSflexibility);
	
	return residueRMSflexibility;
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Calculate flexibility (sigma in Angstroms) of each site, relative to adjacent sites not constrained  
////////////////////////////////////////////////////////////////////////////////
void TIMME::calculateFlexibility() {
	std::cout << " Calculating Flexibility..." << std::endl;
	parameters.mapFromTaskNameToStatus["Calculate flexibility"] = "Running...";
	outputStatus();
	
	for (unsigned int siteID = 1; siteID <= referenceStructure->total_sites; siteID++) {
		// find root mean squared sigma between siteID and all next-next-nearest neighbors 
		calculateFlexibility(siteID);		
	}
	
	parameters.mapFromTaskNameToStatus["Calculate flexibility"] = "Complete";
	outputStatus();
	
	std::cout << "   done" << std::endl;
}

// Default destructor.
////////////////////////////////////////////////////////////////////////////////
TIMME::~TIMME()
{
}
