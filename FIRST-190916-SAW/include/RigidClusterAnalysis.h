#ifndef _RIGID_CLUSTER_ANALYSIS_H
#define _RIGID_CLUSTER_ANALYSIS_H

#include "RigidClusterHierarchy.h"

#include "Plot.h"
#include "SvgPlot.h"
#include "PsPlot.h"
#include "Parameters.h"
#include "MolFramework.h" 
#include "SiteSelector.h"

#include "Cutoff.h"

enum PlotInteractivity {Interactive, NonInteractive};
typedef std::set<RigidCluster *> CollectionOfRigidClusters;

class RigidClusterAnalysis : public SiteSelector {
  public:

  RigidClusterAnalysis(RigidClusterHierarchy * rigidClusterHierarchy);
  void setStructure(MolFramework *structure);
  void setXAxisType(XAxisType xAxisType);
  
  void saveDilution(string filename);
	
  void saveDilutionStripePlot(Plot *plot, PlotInteractivity);
	
  void saveBlockDiagonalDilutionPlot(Plot *plot);

  void saveDilutionIntervals(string filename);
  void saveDilutionData(string filename); // AJR 03.29.06
  void saveJmolDilutionScript(string filename);
	
  void savePymolDilutionScript(string filename);
	
  void saveDilutionStatistics(string filename);
  
  float computeMeanRigidClusterSize(Cutoff cutoff);
  float computeLargestRigidClusterSize(Cutoff cutoff);
  
  float computeMeanRigidClusterSize(CollectionOfRigidClusters);
  float computeLargestRigidClusterSize(CollectionOfRigidClusters);

  unsigned int getRigidLabel(SiteID siteNumber, Cutoff cutoff);
	
  unsigned int getRigidLabel(RigidCluster * rigidCluster, Cutoff cutoff);
	
  void assignRigidClusterLabels();
  
  void assignRigidClusterLabels(RigidCluster *); 
  
  void assignRigidClusterLabels(RigidCluster *, 
                                unsigned int minLabel, 
                                unsigned int maxLabel); 
  
  void computeRigidClusterColoring();
  
  void saveRigidClusterIntervals(ofstream & outputStream,
                                 Cutoff cutoff, 
                                 CollectionOfRigidClusters setOfRigidClusters);
  	
  RigidClusterHierarchy * getRigidClusterHierarchy();
  
  void saveOutput(std::vector<MolFramework*> *vectorOfStructures, 
                  string outputFilename);
  
	void savePDBmodel(MolFramework* model,
                    std::map<SiteID, float> &mapFromSiteIDtoOccupancy,
                    std::map<SiteID, float> &mapFromSiteIDtoCharge,
                    ofstream &outputFileStream);
  
private:
  MolFramework *structure;
  XAxisType xAxisType;
      
//  typedef std::vector<Site_Info *> Sites;
  typedef std::vector<unsigned int> Sites;

  RigidClusterHierarchy * rigidClusterHierarchy;

  float stripePlotWidth;
  float stripePlotHorizontalPadding;
  float stripePlotHeight;
  bool showChainsInStripePlot;
  float verticalOffsetForPlotArea;
  float horizontalOffsetForPlotArea;
  float fractionalSpacingBetweenLines;
  float stripeThickness;
  float separationBetweenStripes;
  float plotAreaHeight;
    
  float yLabelWidth;
  float yTicsWidth;
  float widthPerSite;
  float widthPerResidue;
  float offsetBetweenChains;
  float chartTitleHeight;
  float chainSectionHeight;
  
  float separationBetweenResidues;
	
  float xLabelHeight;	
		
  CollectionOfRigidClusters cleanupCollectionOfRigidClusters(CollectionOfRigidClusters collectionOfRigidClusters);

  bool includeSite(SiteID siteID);
  static bool includeSite(MolFramework *structure, 
                             XAxisType xAxisType,
                             unsigned siteID);
  
  float getXCoordinate(SiteID siteID);
  float getXCoordinate(Residue &residue);
  float getXCoordinate(Site_Info &site);
  float getXCoordinate(Chain &chain);
  float getWidth(SiteID siteNumber);
  float getWidth(Residue &residue);
  float getWidth(Site_Info &site);
  float getWidth(Chain &chain);
  float getYCoordinateFromCutoff(Cutoff cutoff);
  float getStripeHeight(float cutoff, float previousCutoff, RigidLabel previousRigidLabel);

  SiteID getSumOfSiteIDs(RigidCluster *rigidCluster);
  SiteID getMeanSiteID(RigidCluster *rigidCluster);
  
  SiteID getSumOfSiteIDsquaredDeviations(RigidCluster *rigidCluster, SiteID meanValue); 
  SiteID getSiteIDrmsdFromMean(RigidCluster *rigidCluster); 

  void generateSetOfCutoffs();
  void identifyVisibleSites();
  
  void setupPlotDimensions();
  
  void dilutionPlotScript(Plot *plot);
  
  void drawChainBars(Plot *plot);
  void drawXTics(Plot *plot);
  void drawYTics(Plot *plot);
  
  void dilutionPlot(Plot *plot);
  void subclusterBlockDiagonal(RigidCluster * rigidCluster, Plot *plot);
  void blockDiagonalPlot(Plot *plot);
  void dilutionPlotOverlay(Plot *plot);
  
  float getTicSpacing(float intervalLength,
                      unsigned int numberOfTics);
  
  float getClosestCutoffBelowCutoff(Cutoff cutoff);
  
  void preflightForStripePlot();
  
  Sites sites;
  std::set<Cutoff> distinctCutoffs;
};

#endif //_RIGID_CLUSTER_ANALYSIS_H
