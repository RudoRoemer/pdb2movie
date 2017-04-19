#include <typeinfo>
#include <limits.h>
#include "Parameters.h"
#include "RigidClusterAnalysis.h"
#include "RigidLabel.h"
#include "flexweb.h"
#include "Color.h"

extern Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Constructor for a RigidClusterAnalysis of a RigidClusterHierarchy
// 
////////////////////////////////////////////////////////////////////////////////
RigidClusterAnalysis::RigidClusterAnalysis(RigidClusterHierarchy * rigidClusterHierarchy) {
  this->rigidClusterHierarchy = rigidClusterHierarchy;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// true if firstRigidCluster is larger than secondRigidCluster
// false, otherwise
// 
////////////////////////////////////////////////////////////////////////////////
bool sortClustersBySize (RigidCluster * firstRigidCluster, 
    RigidCluster * secondRigidCluster) {
  return firstRigidCluster->size() > secondRigidCluster->size();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Set the MolFramework structure to analyze. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::setStructure(MolFramework *structure) {
  this->structure = structure;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::setXAxisType(XAxisType xAxisType) {
  this->xAxisType = xAxisType;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Save data files and generated plots / molecular visualization representations 
// for the dilution hierarchy contained in a RigidClusterAnalysis object.
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveDilution(string prefix) {

  int last_slash = prefix.find_last_of("/\\");
  string basename = prefix.substr( last_slash+1 );

  preflightForStripePlot();

  {
    string dilutionStripePlotFilename = prefix + ".ps";
    Plot *plot = new PsPlot(dilutionStripePlotFilename);
    saveDilutionStripePlot(plot, NonInteractive);  
  }

  if (parameters.flexweb) {
    string dilutionStripePlotFilename = prefix + "_interactive.svg";
    Plot *plot = new SvgPlot(dilutionStripePlotFilename);
    saveDilutionStripePlot(plot, Interactive);

    string dilutionStripePlotPlaceholderFilename = prefix + ".html";

    // FIXME - remove the folowing (tight coupling between FlexWeb and FIRST):
    ofstream dilutionStripePlotPlaceholder(dilutionStripePlotPlaceholderFilename.c_str());

    // FIXME - uncomment the next line and remove the following
    // the IE 7 preview won't render the SVG with type='image/svg+xml' present :-(
    //    dilutionStripePlotPlaceholder << "<object type='image/svg+xml' data=\""<< dilutionStripePlotFilename << "\" width=\"" << stripePlotWidth << "\" height=\"" << (stripePlotHeight + 6.) << "\">" << std::endl; 

    string dilutionStripePlotBasename = basename + "_interactive.svg";

    dilutionStripePlotPlaceholder << "<object data=\""<< dilutionStripePlotBasename << "\" width=\"" << stripePlotWidth << "\" height=\"" << (stripePlotHeight + 6.) << "\">" << std::endl;
    dilutionStripePlotPlaceholder << "Your browser doesn't appear to support SVG." << std::endl;
    dilutionStripePlotPlaceholder << "This tile should work with Internet Explorer 7, Firefox 1.5 or Safari 2.1 or greater." << std::endl;
    dilutionStripePlotPlaceholder << "Older browsers should be able to use the <a href='http://www.adobe.com/svg/viewer/install/main.html'>Adobe SVG Plugin</a>";
    dilutionStripePlotPlaceholder << "</object>";
    dilutionStripePlotPlaceholder.close();
  }

  {
    string dilutionStripePlotFilename = prefix + ".svg";
    Plot *plot = new SvgPlot(dilutionStripePlotFilename);
    saveDilutionStripePlot(plot, NonInteractive);

    string dilutionBlockDiagonalFilename = prefix + ".blockDiagonal.svg";
    Plot *blockDiagonalPlot = new SvgPlot(dilutionBlockDiagonalFilename);
    saveBlockDiagonalDilutionPlot(blockDiagonalPlot);
  }  

  if (!parameters.run_timme) {
    saveDilutionData(prefix); // AJRader 03.30.06
  }

  saveDilutionIntervals(prefix);
  savePymolDilutionScript(prefix);
  saveDilutionStatistics(prefix);

  if (parameters.flexweb) {
    saveJmolDilutionScript(prefix);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Save the rigid label intervals for each rigid cluster at a given cutoff 
// of a rigid cluster hierarchy. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveDilutionIntervals(string prefix) {
  string dilutionCutoffsFilename = prefix + "_cutoffs.js"; 
  ofstream dilutionCutoffsFile(dilutionCutoffsFilename.c_str());

  string dilutionIntervalsFilename = prefix + "_intervals_pymol"; // FIXME - use a language agnostic extension (this is here to appease Python's import)
  ofstream dilutionIntervalsFile(dilutionIntervalsFilename.c_str());

  dilutionIntervalsFile << "mapFromCutoffToNewRigidClusterIntervals = {};" << std::endl;
  dilutionCutoffsFile << "cutoffs = [];" << std::endl;
  CollectionOfRigidClusters previousSetOfLargestRigidClusters;
  for (std::map<float, CollectionOfRigidClusters >::reverse_iterator cutoffToCollectionOfLargestRigidClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rbegin();
      cutoffToCollectionOfLargestRigidClusters != rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rend();
      cutoffToCollectionOfLargestRigidClusters++) {
    std::pair<float, CollectionOfRigidClusters > pair = (*cutoffToCollectionOfLargestRigidClusters);

    float cutoff = pair.first;
    dilutionCutoffsFile << "cutoffs.push(" << cutoff << ");" << std::endl;

    CollectionOfRigidClusters setOfLargestRigidClusters = pair.second;
    CollectionOfRigidClusters setOfNewRigidClusters;
    set_difference(setOfLargestRigidClusters.begin(), setOfLargestRigidClusters.end(),
        previousSetOfLargestRigidClusters.begin(), previousSetOfLargestRigidClusters.end(),
        inserter(setOfNewRigidClusters, setOfNewRigidClusters.begin()));

    dilutionIntervalsFile << "mapFromCutoffToNewRigidClusterIntervals" << "[" << cutoff << "] = ";
    saveRigidClusterIntervals(dilutionIntervalsFile, cutoff, setOfNewRigidClusters);
    setOfNewRigidClusters.clear();

    dilutionIntervalsFile << ";" << std::endl;

    previousSetOfLargestRigidClusters = setOfLargestRigidClusters;
  }

  float maximumCutoff = (*rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rbegin()).first;

  dilutionIntervalsFile << "mapFromRigidLabelToLargestInterval = ";

  saveRigidClusterIntervals(dilutionIntervalsFile, maximumCutoff, rigidClusterHierarchy->setOfLargestRigidClusters);
  dilutionIntervalsFile << ";" << std::endl;

  dilutionIntervalsFile.close();
  dilutionCutoffsFile.close();
} 

void RigidClusterAnalysis::preflightForStripePlot() {
  static bool hasRun = false; 

  if (hasRun) {
    return;
  }

  structure->populateStructureHierarchy();
  generateSetOfCutoffs();
  setupPlotDimensions();

  hasRun = true;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Generate and save the stripe plot portion of a dilution plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveDilutionStripePlot(Plot * plot, PlotInteractivity plotInteractivity) {

  // TODO - removeme - refactor this to joinRigidCluster? (need to assign rigid labels on build, though) 
  if (parameters.publicationStyle) {
    showChainsInStripePlot = false;
  } else {
    showChainsInStripePlot = true;
  }

  parameters.mapFromTaskNameToStatus["Generate dilution plot"] = "Running...";
  outputStatus();

  stringstream height;
  height << stripePlotHeight + 6.;
  stringstream width;
  width << stripePlotWidth; 

  plot->beginDocument((int)stripePlotWidth, (int)stripePlotHeight + 6);

  dilutionPlotScript(plot);

  if (plotInteractivity == NonInteractive) {
    drawYTics(plot);
  }

  drawChainBars(plot);

  identifyVisibleSites();

  drawXTics(plot);
  dilutionPlot(plot);

  if (plotInteractivity == Interactive) {
    dilutionPlotOverlay(plot);
  }

  plot->endDocument();

  parameters.mapFromTaskNameToStatus["Generate dilution plot"] = "Complete";
  outputStatus();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns: the collectionOfRigidClusters of rigid clusters that include more than 
// one site
// 
////////////////////////////////////////////////////////////////////////////////
CollectionOfRigidClusters RigidClusterAnalysis::cleanupCollectionOfRigidClusters(CollectionOfRigidClusters collectionOfRigidClusters) {
  CollectionOfRigidClusters cleanCollectionOfRigidClusters;

  for (CollectionOfRigidClusters::iterator rigidClusterIterator = collectionOfRigidClusters.begin();
      rigidClusterIterator != collectionOfRigidClusters.end();
      rigidClusterIterator++) {
    RigidCluster* rigidCluster = *rigidClusterIterator;

    size_t rigidClusterSize = rigidCluster->size(this);

    if (rigidClusterSize > 1) { // FIXME - hard-coded constant
      cleanCollectionOfRigidClusters.insert(rigidCluster);
    } 
  }

  return cleanCollectionOfRigidClusters;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Identify the set of cutoffs to include in a dilution plot. 
// A cutoff should only be included if passing through it visibly changes the 
// rigid cluster decomposition for visible sites in a dilution plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::generateSetOfCutoffs() {
  static bool hasRun = false;

  if (hasRun) {
    return;
  }

  Cutoff maximumObservedCutoff = *rigidClusterHierarchy->setOfCutoffs.rbegin();

  std::map<float, CollectionOfRigidClusters >::reverse_iterator previousCutoffToCollectionOfLargestRigidClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rbegin();
  std::pair<float, CollectionOfRigidClusters > previousPair = (*previousCutoffToCollectionOfLargestRigidClusters);
  CollectionOfRigidClusters previousCollectionOfRigidClusters = cleanupCollectionOfRigidClusters(previousPair.second);

  Cutoff previousCutoff = maximumObservedCutoff;
  Cutoff minimumSeparationBetweenCutoffs = 0.001; // FIXME - determine this algorithmically 

  for (std::map<float, CollectionOfRigidClusters >::reverse_iterator cutoffToCollectionOfLargestRigidClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rbegin();
      cutoffToCollectionOfLargestRigidClusters != rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rend();
      cutoffToCollectionOfLargestRigidClusters++) {

    std::pair<float, CollectionOfRigidClusters > pair = (*cutoffToCollectionOfLargestRigidClusters);
    Cutoff cutoff = pair.first;
    CollectionOfRigidClusters collectionOfLargestRigidClusters = cleanupCollectionOfRigidClusters(pair.second);

    /*    CollectionOfRigidClusters differenceInRigidClusters;
          set_symmetric_difference(collectionOfLargestRigidClusters.begin(), collectionOfLargestRigidClusters.end(), 
          previousCollectionOfRigidClusters.begin(), previousCollectionOfRigidClusters.end(),
          inserter(differenceInRigidClusters, differenceInRigidClusters.begin()));

          if (differenceInRigidClusters.size() > 0) {
    // TODO - iterate through differenceInRigidClusters and only insert if the
    // difference includes at least one cluster with more than one site
    distinctCutoffs.insert(cutoff);
    }*/

    if (previousCollectionOfRigidClusters.size() != collectionOfLargestRigidClusters.size()) {      
      // FIXME - only base this on the visible sites 
      Cutoff separationBetweenCutoffs = previousCutoff - cutoff;

      if (minimumSeparationBetweenCutoffs < separationBetweenCutoffs) {
        distinctCutoffs.insert(cutoff);
        previousCutoff = cutoff;

      }
    }

    previousCollectionOfRigidClusters = collectionOfLargestRigidClusters;
  }

  hasRun = true;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Compute the dimensions of a dilution plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::setupPlotDimensions() {
  static bool hasRun = false;
  if (hasRun) {
    return;
  }

  verticalOffsetForPlotArea = 0.;
  horizontalOffsetForPlotArea = 0.;

  // TODO - make these user adjustable options
  stripePlotHorizontalPadding = 30.0;
  yLabelWidth = 20;
  yTicsWidth = 100;
  widthPerSite = .5;
  widthPerResidue = 2.5;
  offsetBetweenChains = 10.;
  chartTitleHeight = 50;
  chainSectionHeight = 30;
  xLabelHeight = 15;

  stripeThickness = parameters.stripeThickness;
  separationBetweenStripes = parameters.separationBetweenStripes;

  horizontalOffsetForPlotArea += yLabelWidth + yTicsWidth;

  stripePlotWidth = horizontalOffsetForPlotArea + stripePlotHorizontalPadding;

  switch(xAxisType) {
    case perResidue: 
      for (std::map <unsigned int, Chain*>::iterator chainIterator = structure->mapFromChainIDtoChain.begin();
          chainIterator != structure->mapFromChainIDtoChain.end();
          chainIterator++) {
        Chain *chain = chainIterator->second;

        stripePlotWidth += getWidth(*chain);
        stripePlotWidth += offsetBetweenChains;
      }

      break;

    case perSite:
      stripePlotWidth += widthPerSite * structure->total_sites;
      break;

    case perRigidLabel: 
      stripePlotWidth += widthPerSite * structure->total_sites;
      break;
  }

  verticalOffsetForPlotArea += chartTitleHeight + xLabelHeight;

  if (showChainsInStripePlot) {
    verticalOffsetForPlotArea += chainSectionHeight;
  }

  plotAreaHeight = (stripeThickness + separationBetweenStripes) * (float) distinctCutoffs.size();
  //  plotAreaHeight = (stripeThickness) * (float) distinctCutoffs.size();
  stripePlotHeight = verticalOffsetForPlotArea + plotAreaHeight;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Draw the y-axis tics for a dilution plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::drawYTics(Plot *plot) {
  Point yAxisLabelPoint;
  yAxisLabelPoint.x = 0.0; // FIXME - hard-coded constant
  float yAxisLabelyPosition = (verticalOffsetForPlotArea+(stripePlotHeight-verticalOffsetForPlotArea)/2.);
  yAxisLabelPoint.y = yAxisLabelyPosition + 20.0;

  Style yAxisLabelFontStyle;
  yAxisLabelFontStyle["text-anchor"] = "middle";
  yAxisLabelFontStyle["fill"] = "black";
  yAxisLabelFontStyle["font-size"] = "14pt";
  yAxisLabelFontStyle["font-family"] = "serif";
  yAxisLabelFontStyle["font-weight"] = "normal";

  Attributes attributes;
  stringstream yAxisTransform;
  yAxisTransform << "translate(15) rotate(270 0 " << yAxisLabelyPosition << ")";

  attributes["transform"] = yAxisTransform.str();

  string cutoffUnits = "&#197;"; // FIXME - make this a parameter so we can have Acircle or kcal/mol

  stringstream yAxisLabel;
  yAxisLabel << "Cutoff, &#x03C3;";
  
  if (parameters.timmeScaleFactorA != 1.0) {
    yAxisLabel << " / a";// (" << cutoffUnits << ")";
  }

  plot->drawText(yAxisLabelPoint, yAxisLabelFontStyle, attributes, yAxisLabel.str());

  float minimumCutoff = *distinctCutoffs.begin(); 
  float maximumCutoff = *distinctCutoffs.rbegin();

  float deltaCutoff = (maximumCutoff - minimumCutoff) / 5.0; // FIXME - hard-coded constant

  for (float cutoff = minimumCutoff + deltaCutoff; cutoff <= maximumCutoff; cutoff += deltaCutoff) {
    stringstream yAxisTicLabel;

    yAxisTicLabel << setfill(' ') << std::setw(6) << fixed << setprecision(2) << (cutoff / parameters.timmeScaleFactorA) << " - ";

    Point yTicPoint;
    yTicPoint.x = horizontalOffsetForPlotArea - 10.0;

    yTicPoint.y = getYCoordinateFromCutoff(cutoff) + 14.0 / 2.0; // FIXME - hard-coded constant

    Style yTicStyle;
    yTicStyle["fill"] = "black"; 
    //		yTicStyle["stroke"] = "white";
    //    stringstream fontSize;
    //    fontSize << stripeThickness << "px";
    yTicStyle["font-size"] = "14pt"; // fontSizedt.str();
    yTicStyle["font-family"] = "serif";
    yTicStyle["alignment-baseline"] = "top";
    yTicStyle["font-weight"] = "normal";
    yTicStyle["text-anchor"] = "end";

    plot->drawText(yTicPoint, yTicStyle, yAxisTicLabel.str());
  } 
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Draw the representation of multiple chains in a dilution plot.  
//
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::drawChainBars(Plot *plot){
  unsigned int currentChainID = structure->site_info[1].FIRST_chain_ID;

  float chainVerticalPosition = 10;

  if (showChainsInStripePlot) {

    for (std::map <unsigned int, Chain*>::iterator chainIterator = structure->mapFromChainIDtoChain.begin();
        chainIterator != structure->mapFromChainIDtoChain.end();
        chainIterator++) {
      Chain *chain = chainIterator->second;

      switch(xAxisType) {
        case perSite: 

          break;

        case perRigidLabel: 

          break;

        case perResidue:
          // draw a rectangle from firstSiteIDInThisChain to this siteID and color by chainID
          Chain::iterator residueInChainIterator = chain->begin();

          Residue *firstResidueInChain = residueInChainIterator->second;
          Site_Info *siteInfo = *firstResidueInChain->beginSites();
          currentChainID = siteInfo->FIRST_chain_ID;					

          Point chainBoxPoint;
          chainBoxPoint.x = getXCoordinate(*chain);
          chainBoxPoint.y = chainVerticalPosition;

          Rectangle chainBoxRectangle;
          chainBoxRectangle.point = chainBoxPoint;
          chainBoxRectangle.width = getWidth(*chain);
          chainBoxRectangle.height = (25); // FIXME - hardcoded constant

          Style chainBoxPlotStyle;
          chainBoxPlotStyle["fill"] = plot->getColor(currentChainID);
          chainBoxPlotStyle["stroke"] = "black";
          chainBoxPlotStyle["opacity"] = ".2"; // FIXME - magic number 

          plot->drawRectangle(chainBoxPlotStyle, chainBoxRectangle);

          stringstream chainBoxLabel;

          if (structure->mapFromChainIDtoChain.size()==1) {
            chainBoxLabel << structure->base_name << " "; 

          } else {

            char chainLabel = (char)('A' + currentChainID); // FIXME - assumption about character encoding
            chainBoxLabel << (chainLabel);
          };

          Point chainBoxLabelPoint;
          chainBoxLabelPoint.x = (getXCoordinate(*chain) +
              (getWidth(*chain))/2.);

          chainBoxLabelPoint.y = chainVerticalPosition + 15; // FIXME - hardcoded constant

          Style chainBoxLabelStyle;

          chainBoxLabelStyle["fill"] = "black";
          chainBoxLabelStyle["stroke"] = "none";
          chainBoxLabelStyle["font-size"] = "14pt";
          chainBoxLabelStyle["font-family"] = "roman";
          chainBoxLabelStyle["font-weight"] = "normal";
          chainBoxLabelStyle["text-anchor"] = "middle";

          plot->drawText(chainBoxLabelPoint, chainBoxLabelStyle, chainBoxLabel.str());

          break;
      }			
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Draw the x-axis tics for a dilution plot
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::drawXTics(Plot *plot){
  float xTicVerticalPadding = 2.;
  float xTicFontSize = 12.;

  stringstream xAxisLabel;

  switch (xAxisType) {
    case perResidue:
      xAxisLabel << "Residue #";

      break;

    case perSite:
      xAxisLabel << "Site #";

      break;

    case perRigidLabel: 
      xAxisLabel << "Rigid Label";

      break;
  }

  Point xAxisLabelPoint;
  xAxisLabelPoint.x = (stripePlotWidth)/2. + 50.;// stripePlotHorizontalPadding; // FIXME - hard coded constant
  xAxisLabelPoint.y =  (verticalOffsetForPlotArea - xTicFontSize - xTicVerticalPadding - 1.5 * 25.); // FIXME - hard coded constant

  Style xAxisLabelStyle;
  xAxisLabelStyle["fill"] = "black";
  xAxisLabelStyle["stroke"] = "none";
  xAxisLabelStyle["font-size"] = "14pt";
  xAxisLabelStyle["font-family"] = "serif";
  xAxisLabelStyle["font-weight"] = "normal";
  xAxisLabelStyle["text-anchor"] = "middle";

  plot->drawText(xAxisLabelPoint,  xAxisLabelStyle, xAxisLabel.str());

  bool firstTic = true;

  unsigned int ticFrequency = 1;
  switch(xAxisType) {
    case perRigidLabel: {

                        } break;

    case perSite: {

                  } break;

    case perResidue: {
                       unsigned int totalNumberOfResidues = structure->mapFromResidueToResidueID.size();
                       unsigned int preferredNumberOfXTics = 6;

                       ticFrequency = (unsigned int)getTicSpacing(totalNumberOfResidues, preferredNumberOfXTics);

                     } break;			
  }

  for (std::map <unsigned int, Chain*>::iterator chainIterator = structure->mapFromChainIDtoChain.begin();     
      chainIterator != structure->mapFromChainIDtoChain.end();
      chainIterator++) {
    Chain *chain = chainIterator->second;

    for (Chain::iterator residueIterator = chain->begin();
        residueIterator != chain->end(); 
        residueIterator++) {

      Residue * residue = residueIterator->second;

      switch(xAxisType) {
        case perRigidLabel: {

                            } break;

        case perSite: {

                      } break;

        case perResidue: {

                           unsigned int residueID = structure->mapFromResidueToResidueID[residue];

                           if (residueID % ticFrequency == 0 ) {//|| 
                               //firstTic) { // TODO - parameterize this and algorithmically select which residueIDs to display
                             firstTic = false;
                             Site_Info *firstSiteInResidue = *residue->beginSites();

                             stringstream xAxisTicLabel;
                             xAxisTicLabel << setfill(' ') << std::setw(6) << fixed << setprecision(4) << firstSiteInResidue->seq_number; // TODO - replace ticNumber with the proper label
                             {
                               Point xTicPoint;
                               xTicPoint.x = getXCoordinate(*residue);//+getWidth(*residue) - 14.0/2.;
                               xTicPoint.y = (verticalOffsetForPlotArea - xTicVerticalPadding - 15.);

                               Style xTicStyle;
                               xTicStyle["fill"] = "black";
                               xTicStyle["stroke"] = "none";
                               stringstream xAxisLabelFontSize;
                               xTicStyle["font-size"] = "14pt";
                               xTicStyle["font-family"] = "serif";
                               xTicStyle["font-weight"] = "normal";
                               xTicStyle["text-anchor"] = "middle";

                               Attributes attributes;
                               stringstream transform;
                               //						transform << "translate(0) rotate(270 " << (xTicPoint.x) << " " <<  xTicPoint.y << ")";
                               attributes["transform"] = transform.str();

                               plot->drawText(xTicPoint, xTicStyle, attributes, xAxisTicLabel.str());

                             } {
                              Point xTicPoint;
                              xTicPoint.x = getXCoordinate(*residue)+getWidth(*residue);
                              xTicPoint.y = (verticalOffsetForPlotArea - xTicVerticalPadding - 7);

                              Style xTicStyle;
                              xTicStyle["fill"] = "black";
                              xTicStyle["stroke"] = "none";
                              stringstream xAxisLabelFontSize;xTicStyle["font-size"] = "14pt";
                              xTicStyle["font-family"] = "serif";
                              xTicStyle["text-anchor"] = "middle";
                              Attributes attributes;
                              stringstream transform;
                              transform << "translate(0) rotate(270 " << (xTicPoint.x) << " " <<  xTicPoint.y << ")";
                              attributes["transform"] = transform.str();
                              
                              plot->drawText(xTicPoint, xTicStyle, attributes, " - ");
                              
                             }
                           }

                           break;
                         }
      }
    }
  }
}

float RigidClusterAnalysis::getStripeHeight(float cutoff, float previousCutoff, RigidLabel previousRigidLabel) {
  float stripeHeight = stripeThickness;
  stripeHeight = (getYCoordinateFromCutoff(previousCutoff) - getYCoordinateFromCutoff(cutoff));
  stripeHeight = max(0.0f, stripeHeight - separationBetweenStripes);

  if ((previousRigidLabel == UINT_MAX)) {
    stripeHeight = .1*stripeThickness; // FIXME - add a parameter to allow the user to specify the stripeThickness for single-site clusters
  }

  return stripeHeight;
}

void RigidClusterAnalysis::subclusterBlockDiagonal(RigidCluster * rigidCluster, Plot *plot) {
  SiteID minLabel = rigidCluster->minLabel;
  SiteID maxLabel = rigidCluster->maxLabel;

  Point clusterRectanglePoint;
  clusterRectanglePoint.x = minLabel;
  clusterRectanglePoint.y = minLabel;
  
  Rectangle clusterRectangle;
  clusterRectangle.point = clusterRectanglePoint;
  clusterRectangle.height = maxLabel - minLabel;
  clusterRectangle.width  = maxLabel - minLabel;
  
  RigidLabel rigidLabel = rigidCluster->rigidLabel;
  
  ColorNumber colorNumber(rigidLabel);
  
  Style clusterRectanglePlotStyle;
  clusterRectanglePlotStyle["fill"] = plot->getColor(rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[colorNumber]);
  
  if (rigidLabel > 0) {
    plot->drawRectangle(clusterRectanglePlotStyle, clusterRectangle);
  }

  if (rigidCluster->beginSubClusters() == rigidCluster->endSubClusters()) {
    return;
  }

  for (RigidCluster::RigidClusterIterator subclusterIterator = rigidCluster->beginSubClusters();
      subclusterIterator != rigidCluster->endSubClusters();
      subclusterIterator ++ ) {
  
    RigidCluster * subCluster = *subclusterIterator;

    subclusterBlockDiagonal(subCluster, plot);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
//  Generate and save a block-diagonalized dilution plot 
//
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::blockDiagonalPlot(Plot *plot) {

  for (std::set<RigidCluster *>::iterator clusterIterator = rigidClusterHierarchy->setOfLargestRigidClusters.begin();
    clusterIterator != rigidClusterHierarchy->setOfLargestRigidClusters.end();
    clusterIterator ++) {

    RigidCluster * rigidCluster = *clusterIterator;

    subclusterBlockDiagonal(rigidCluster, plot);
  }
}

void RigidClusterAnalysis::saveBlockDiagonalDilutionPlot(Plot *plot) {
  size_t numberOfLabels = rigidClusterHierarchy->mapFromRigidLabelToCluster.size();
  plot->beginDocument((int)numberOfLabels, (int)numberOfLabels);
  
  blockDiagonalPlot(plot);
  
  plot->endDocument();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Generate and save a dilution plot to plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::dilutionPlot(Plot *plot){
  RigidLabel previousRigidClusterLabel = UINT_MAX; // TODO - replace UINT_MAX with 'unlabeled'

  float previousCutoff = *distinctCutoffs.begin();
  for (std::set<Cutoff>::iterator cutoffIterator = distinctCutoffs.begin();
      cutoffIterator != distinctCutoffs.end();
      cutoffIterator++) {

    float cutoff = *cutoffIterator;
    Sites::iterator siteIterator = sites.begin();
    if (siteIterator == sites.end()) {
      std::cerr << " error: siteIterator == sites.end()" << std::endl;
      continue; 
    }

    SiteID previousSiteID = *siteIterator;
    Site_Info previousSite = structure->site_info[previousSiteID];

    while (siteIterator != sites.end()) {        
      SiteID siteID = *siteIterator;
      Site_Info site = structure->site_info[siteID];

      RigidLabel rigidClusterLabel = getRigidLabel(siteID, cutoff);

      Point clusterRectanglePoint;

      clusterRectanglePoint.x = getXCoordinate(previousSite);

      clusterRectanglePoint.y = getYCoordinateFromCutoff(cutoff);

      float stripeHeight = getStripeHeight(cutoff, previousCutoff, previousRigidClusterLabel);

      if ((previousRigidClusterLabel == UINT_MAX)) {
        clusterRectanglePoint.y += .5*stripeThickness;
      } 
      
      Sites::iterator nextSiteIterator = siteIterator;
      nextSiteIterator++;

      if (nextSiteIterator != sites.end()) {
        SiteID nextSiteID = *nextSiteIterator;
        Site_Info nextSite = structure->site_info[nextSiteID];

        if (nextSite.chain_ID != site.chain_ID) {
          Rectangle clusterRectangle;

          clusterRectangle.point = clusterRectanglePoint;
          clusterRectangle.height = stripeHeight;

          clusterRectangle.width = getXCoordinate(site) - getXCoordinate(previousSite);

          Style clusterRectanglePlotStyle;
          ColorNumber colorNumber(previousRigidClusterLabel);
          clusterRectanglePlotStyle["fill"] = plot->getColor(rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[colorNumber]);

          if ((previousRigidClusterLabel != UINT_MAX)) { // FIXME - && parameters->hideLines
            plot->drawRectangle(clusterRectanglePlotStyle, clusterRectangle);
          }

          previousSiteID = nextSiteID;
          previousSite = nextSite;

          previousRigidClusterLabel = rigidClusterLabel;
        }
      }

      if ((previousRigidClusterLabel != rigidClusterLabel) || 
          (nextSiteIterator == sites.end())) {
        Rectangle clusterRectangle;

        clusterRectangle.point = clusterRectanglePoint;
        clusterRectangle.height = stripeHeight;

        clusterRectangle.width = getXCoordinate(site) - getXCoordinate(previousSite);

        Style clusterRectanglePlotStyle;
        ColorNumber previousColorNumber(previousRigidClusterLabel);
        clusterRectanglePlotStyle["fill"] = plot->getColor(rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[previousColorNumber]);

        if ((previousRigidClusterLabel != UINT_MAX)) { // FIXME - && parameters->hideLines
          plot->drawRectangle(clusterRectanglePlotStyle, clusterRectangle);
        }

        previousSiteID = siteID;
        previousSite = site;
      }

      previousRigidClusterLabel = rigidClusterLabel;

      siteIterator++;
    }

    previousCutoff = cutoff;
  }
}

// borrowed from Brandon's output_RCD_PDB_format method in MolFramework_Output 
// TODO - refactor to MolFramework_Output (DRY!)
void RigidClusterAnalysis::savePDBmodel(MolFramework* model,
    std::map<SiteID, float> &mapFromSiteIDtoOccupancy,
    std::map<SiteID, float> &mapFromSiteIDtoCharge,
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
        << setw(6) << setprecision(2) << (((float)model->site_info[siteID].rigid_label)/100.) // FIXME - magic number and a bit ugly / cryptic
        << "            "
        //			<< setw(4) << "      " // TODO - segID
        //			<< setw(2) << "  " // TODO - something like siteInfo.element
        //			<< setw(2) << setprecision(1) << mapFromSiteIDtoCharge[siteID] // TODO - only integers!
        << endl; 
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Determine which sites should be included in the dilution plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::identifyVisibleSites() {
  static bool hasRun = false;

  if (hasRun) {
    return;
  }

  for (SiteID siteID = 0; // FIXME - poor choice of name as this may be a rigidLabel and edge issue - do we want to start at 0 for siteIDs?

      siteID < structure->total_sites; // FIXME - boundary case issue - it seems that this should be <= but that gives a core dump 
      siteID++) {

    switch (xAxisType) {
      case perSite:

        if (includeSite(siteID)) {
          sites.push_back(siteID);
        }
        break;

      case perRigidLabel: {

                            unsigned int actualSiteId = rigidClusterHierarchy->mapFromRigidLabelToSiteId[siteID];

                            if (includeSite(actualSiteId)) {
                              sites.push_back(actualSiteId);
                            }

                          } break;

      case perResidue: {

                         if (includeSite(siteID)) {
                           sites.push_back(siteID);
                         }

                       } break;
    }

  }

  hasRun = true;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Assign a rigid label to each site.
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::assignRigidClusterLabels() {
  parameters.mapFromTaskNameToStatus["Assign rigid labels"] = "Running...";
  outputStatus();

  size_t rigidLabelOffset = 0;

  for (CollectionOfRigidClusters::iterator largestRigidClusterIterator = rigidClusterHierarchy->setOfLargestRigidClusters.begin();
      largestRigidClusterIterator != rigidClusterHierarchy->setOfLargestRigidClusters.end();
      largestRigidClusterIterator++) {
    RigidCluster * rigidCluster = *largestRigidClusterIterator;

    //    assignRigidClusterLabels(rigidCluster);

    size_t rigidClusterSize = rigidCluster->size();
    assignRigidClusterLabels(rigidCluster, 
        rigidLabelOffset, 
        rigidLabelOffset + rigidClusterSize);

    rigidLabelOffset += rigidClusterSize;
  }

  parameters.mapFromTaskNameToStatus["Assign rigid labels"] = "Complete";
  outputStatus();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// by default, assign an integer rigid cluster label in the interval [0, # of sites) for each siteID
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::assignRigidClusterLabels(RigidCluster * rigidCluster) {

  assignRigidClusterLabels(rigidCluster, 
      0, 
      rigidCluster->size());
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// assign rigid cluster labels on the interval [minLabel, maxLabel)
//
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::assignRigidClusterLabels(RigidCluster * rigidCluster,
    SiteID minLabel, 
    SiteID maxLabel) { 
  size_t rigidClusterSize = rigidCluster->size();

  rigidCluster->minLabel = minLabel;
  rigidClusterHierarchy->setOfSizeAndRigidLabelPairs.insert(std::pair<size_t, unsigned int>(rigidClusterSize, minLabel));

  if (rigidClusterHierarchy->mapFromRigidLabelToCluster.find(minLabel) != rigidClusterHierarchy->mapFromRigidLabelToCluster.end()) {
    RigidCluster * candidateCluster = rigidClusterHierarchy->mapFromRigidLabelToCluster[minLabel];

    if (rigidCluster->size() > candidateCluster->size()) {
      rigidClusterHierarchy->mapFromRigidLabelToCluster[minLabel] = rigidCluster;
      rigidCluster->rigidLabel = minLabel;
      
      //assert(false); // FIXME - we should never reach here ; if we do, something is wrong
    }

  } else {
    rigidClusterHierarchy->mapFromRigidLabelToCluster[minLabel] = rigidCluster;
    rigidCluster->rigidLabel = minLabel;
    
    //RigidCluster * testCluster = rigidClusterHierarchy->mapFromRigidLabelToCluster[minLabel];
    SiteID meanSiteID = getMeanSiteID(rigidCluster);
    SiteID rmsdFromSiteID  = getSiteIDrmsdFromMean(rigidCluster);

    rigidClusterHierarchy->mapFromRigidLabelToMeanSiteID[minLabel] = meanSiteID;
    rigidClusterHierarchy->mapFromRigidLabelToRMSDfromMeanSiteID[minLabel] = rmsdFromSiteID;
  
    rigidClusterHierarchy->mapFromMaxLabelToMeanSiteID[maxLabel] = meanSiteID;
    rigidClusterHierarchy->mapFromMaxLabelToRMSDfromMeanSiteID[maxLabel] = rmsdFromSiteID;
  }

  rigidCluster->maxLabel = maxLabel;

  unsigned int currentMinLabel = minLabel;
  unsigned int deltaLabel = (maxLabel - minLabel)/(rigidCluster->size());

  stable_sort(rigidCluster->beginSubClusters(), rigidCluster->endSubClusters(), sortClustersBySize);

  bool isLargestSubcluster = true;

  for (RigidCluster::RigidClusterIterator subclusterIterator = rigidCluster->beginSubClusters();
      subclusterIterator != rigidCluster->endSubClusters();
      subclusterIterator ++) {

    RigidCluster * subCluster = *subclusterIterator;

    if (isLargestSubcluster) {
      // assign the same mean / stddev as the parent cluster to the largest subcluster

      isLargestSubcluster = false;

    }

    unsigned int subClusterLabelSize = deltaLabel * subCluster->size();

    unsigned int currentMaxLabel = currentMinLabel+subClusterLabelSize;

    // assign labels for all sites in each rigid subcluster
    assignRigidClusterLabels(subCluster, 
        currentMinLabel, 
        currentMaxLabel);

    currentMinLabel = currentMaxLabel;
  } 

  unsigned int siteLabel = currentMinLabel;

  for (RigidCluster::SiteIdIterator siteIDs = rigidCluster->beginSiteIDs();
      siteIDs != rigidCluster->endSiteIDs();
      siteIDs ++) {

    SiteID siteID = *siteIDs;

    // assign a rigid_label for each site in this RigidCluster
    structure->site_info[siteID].rigid_label = (unsigned int)siteLabel; // FIXME - should we just move this from floats to unsigned ints or use floats for everything?

    rigidClusterHierarchy->mapFromSiteIdToRigidLabel[siteID] = (unsigned int)siteLabel;
    rigidClusterHierarchy->mapFromRigidLabelToSiteId[(unsigned int)siteLabel] = siteID;
    rigidClusterHierarchy->mapFromRigidLabelToCluster[(unsigned int)siteLabel] = rigidCluster;

    siteLabel ++;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////
SiteID RigidClusterAnalysis::getSumOfSiteIDs(RigidCluster *rigidCluster) {
  SiteID sumOfSiteIDs = 0;

  for (RigidCluster::SiteIdIterator siteIDiterator = rigidCluster->beginSiteIDs();
      siteIDiterator != rigidCluster->endSiteIDs();
      siteIDiterator ++) {

    SiteID siteID = *siteIDiterator;

    sumOfSiteIDs += siteID;   
  }

  for (RigidCluster::RigidClusterIterator subClusterIterator = rigidCluster->beginSubClusters(); 
      subClusterIterator != rigidCluster->endSubClusters();
      subClusterIterator ++) {
    RigidCluster * subCluster = *subClusterIterator;

    sumOfSiteIDs += getSumOfSiteIDs(subCluster);
  }

  return (sumOfSiteIDs);
}

////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////
SiteID RigidClusterAnalysis::getMeanSiteID(RigidCluster *rigidCluster) {
  SiteID sumOfSiteIDs = getSumOfSiteIDs(rigidCluster);

  SiteID meanSiteID = (SiteID) (sumOfSiteIDs / (float)rigidCluster->size()); 

  return meanSiteID;
}

////////////////////////////////////////////////////////////////////////////////
//
//
////////////////////////////////////////////////////////////////////////////////
SiteID RigidClusterAnalysis::getSumOfSiteIDsquaredDeviations(RigidCluster *rigidCluster, SiteID meanValue) {
  SiteID sumOfSquaredDeviations = 0;

  for (RigidCluster::SiteIdIterator siteIDiterator = rigidCluster->beginSiteIDs();
      siteIDiterator != rigidCluster->endSiteIDs();
      siteIDiterator ++) {

    SiteID siteID = *siteIDiterator;

    sumOfSquaredDeviations += (SiteID) (pow((double)meanValue - (double)siteID, 2));
  }

  for (RigidCluster::RigidClusterIterator subClusterIterator = rigidCluster->beginSubClusters();
      subClusterIterator != rigidCluster->endSubClusters();
      subClusterIterator ++) {
    RigidCluster * subCluster = *subClusterIterator;

    sumOfSquaredDeviations += getSumOfSiteIDsquaredDeviations(subCluster, meanValue);
  }

  return (sumOfSquaredDeviations);
}

////////////////////////////////////////////////////////////////////////////////
//
// 
////////////////////////////////////////////////////////////////////////////////
SiteID RigidClusterAnalysis::getSiteIDrmsdFromMean(RigidCluster *rigidCluster) {
  SiteID meanSiteID = getMeanSiteID(rigidCluster);

  double sumOfSiteIDsquaredDeviations = (double)getSumOfSiteIDsquaredDeviations(rigidCluster, meanSiteID);

  double rootMeanSquaredDeviation = sqrt(sumOfSiteIDsquaredDeviations / (double)rigidCluster->size());

  SiteID siteIDrmsdFromMean = (SiteID) rootMeanSquaredDeviation;

  return siteIDrmsdFromMean;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Compute which colors to assign to each rigid label in a rigid cluster 
// decomposition.
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::computeRigidClusterColoring(){
  Color::setUseSpectrum(parameters.useSpectralColoring);

  std::vector<unsigned int> rigidLabelsBySize;
  std::set<unsigned int> assignedRigidLabels;

  for (std::set<std::pair<size_t, unsigned int> >::reverse_iterator sizeAndLabelPairIterator = rigidClusterHierarchy->setOfSizeAndRigidLabelPairs.rbegin(); 
      sizeAndLabelPairIterator != rigidClusterHierarchy->setOfSizeAndRigidLabelPairs.rend();
      sizeAndLabelPairIterator ++) {

    RigidLabel rigidLabel =  (unsigned int)(*sizeAndLabelPairIterator).second;
    //		size_t size = (*sizeAndLabelPairIterator).first;

    if (assignedRigidLabels.find(rigidLabel) == assignedRigidLabels.end()) {
      rigidLabelsBySize.push_back(rigidLabel);
      assignedRigidLabels.insert(rigidLabel);
    }
  }

  Color::setNumberOfBins(rigidLabelsBySize.size());
  ColorNumber colorNumber = 0;
  for (std::vector<unsigned int>::iterator rigidLabelsBySizeIterator = rigidLabelsBySize.begin();
      rigidLabelsBySizeIterator != rigidLabelsBySize.end();
      rigidLabelsBySizeIterator ++) {
    RigidLabel rigidLabel = (unsigned int )(*rigidLabelsBySizeIterator);	

    if (parameters.useSpectralColoring) {
      RigidCluster * rigidCluster = rigidClusterHierarchy->mapFromRigidLabelToCluster[rigidLabel];

      SiteID minLabel = rigidCluster->minLabel;

      SiteID meanSiteID = rigidClusterHierarchy->mapFromRigidLabelToMeanSiteID[minLabel];
      SiteID rmsdFromMeanSiteID = rigidClusterHierarchy->mapFromRigidLabelToRMSDfromMeanSiteID[minLabel];
      
      //meanSiteID = rigidClusterHierarchy->mapFromMaxLabelToMeanSiteID[maxLabel];
      //rmsdFromMeanSiteID = rigidClusterHierarchy->mapFromMaxLabelToRMSDfromMeanSiteID[maxLabel];

      rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[rigidLabel] = rigidLabel;
      Color::mapFromColorNumberToMeanSiteID[rigidLabel] = meanSiteID;
      Color::mapFromColorNumberOtRMSDfromMeanSiteID[rigidLabel] = rmsdFromMeanSiteID;
    
    } else {
      rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[rigidLabel] = colorNumber;

    }

    colorNumber ++;
  }

  rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[UINT_MAX] = UINT_MAX;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Generate the SVG user interaction overlay for a dilution plot.
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::dilutionPlotOverlay(Plot *plot) {  
  if (typeid(*plot) != typeid(SvgPlot)) {
    return;
  }

  SvgPlot *svgPlot = (SvgPlot*)plot;

  Chain *firstChain = structure->mapFromChainIDtoChain.begin()->second;
  Chain *lastChain = structure->mapFromChainIDtoChain.rbegin()->second;

  Residue * firstResidue = firstChain->begin()->second;
  Residue * lastResidue = lastChain->rbegin()->second;

  float width = getXCoordinate(*lastResidue)-getXCoordinate(*firstResidue);

  Cutoff  previousCutoff = *distinctCutoffs.begin();
  float previousYCoordinate = getYCoordinateFromCutoff(previousCutoff);
  
  for (std::set<Cutoff>::iterator cutoffIterator = distinctCutoffs.begin();
      cutoffIterator != distinctCutoffs.end();
      cutoffIterator++) {

    float cutoff = *cutoffIterator;

    stringstream link;

    link << "javascript:setCutoff(" << cutoff << ", 1)";
    std::map<string, string> linkAttributes;
    linkAttributes["xlink:href"] = link.str();

    //			scalableVectorGraphics.openTag("a", linkAttributes);

    std::map<string, string> cutoffRectangleAttributes;
    cutoffRectangleAttributes["style"] = "fill: black; ";
    cutoffRectangleAttributes["opacity"] = "0";
    cutoffRectangleAttributes["x"] = "1";

    Point overlayPoint;
    overlayPoint.x = getXCoordinate(0);
    overlayPoint.y = getYCoordinateFromCutoff(cutoff);

    Rectangle overlayRectangle;
    overlayRectangle.point = overlayPoint;
    overlayRectangle.width = width;
    overlayRectangle.height = stripeThickness;

    Attributes overlayAttributes;

    // FIXME - use indices for cutoffs rather than the cutoffs themselves (this is causing lots of unnecessary problems)

    float minimumCutoffMagnitude = 0.000001;

    if (pow(cutoff, 2) < minimumCutoffMagnitude) {
      cutoff = 0;

    }

    stringstream onmouseover;
    onmouseover	<<  "highlight(evt, '" << fixed << cutoff << "')"; // asdf 
    overlayAttributes["onmouseover"]=onmouseover.str();

    stringstream onmouseout;
    onmouseout	<<  "unhighlight(evt, '" << fixed << cutoff << "')";
    overlayAttributes["onmouseout"]=onmouseout.str();

    stringstream onclick;
    onclick	<<  "selectCutoff(evt, '" << fixed << cutoff << "')";
    overlayAttributes["onclick"]=onclick.str();      

    stringstream groupID;
    groupID << fixed << cutoff;
    overlayAttributes["id"]=groupID.str();          

    overlayAttributes["opacity"] = "0.0"; // FIXME - magic number 

    Style overlayStyle;
    overlayStyle["fill"] = "black";
    overlayStyle["stroke"] = "black";
    svgPlot->openTag("g",
        overlayAttributes);

    //    svgPlot->drawRectangle(overlayStyle, overlayRectangle);

    Points points;

    // FIXME - hard coded constants (make user selectable)
    float cutoffTriangleVerticalOffset = 10;
    float cutoffTriangleHorizontalOffset = 10;

    /*    float xCoordinate = horizontalOffsetForPlotArea;

          switch (xAxisType) {
          case perSite: {
          xCoordinate +=  widthPerSite * (float)siteNumber;

          }	brmeak;

          case perRigidLabel: {
          unsigned int rigidLabel = rigidClusterHierarchy->mapFromSiteIdToRigidLabel[siteNumber]; // FIXME 

          xCoordinate +=  widthPerSite * (float)rigidLabel;

          } break;


          case perResidue: {
          if (structure->mapFromSiteIDToResidue.find(siteNumber) != structure->mapFromSiteIDToResidue.end()) {
          Residue *residue = structure->mapFromSiteIDToResidue[siteNumber];
          xCoordinate = getXCoordinate(*residue);

          size_t residueSize = residue->size();

          float widthPerSiteInResidue = widthPerResidue / (float)residueSize;

          unsigned int siteNumberInResidue = 0;
          for (Residue::Sites::iterator siteIterator = residue->beginSites();
          siteIterator != residue->endSites();
          siteIterator ++) {
          Site_Info * site = *siteIterator;
          if (site->site_number == siteNumber) {
          break;
          }

          siteNumberInResidue++;
          }
          xCoordinate += (float) siteNumberInResidue * widthPerSiteInResidue; 
          }

          } break;
          }*/

    Point topLeft;
    topLeft.x = overlayPoint.x - cutoffTriangleHorizontalOffset;
    topLeft.y = overlayPoint.y - cutoffTriangleVerticalOffset;
    points.push_back(topLeft);

    Point topMiddleLeft;
    topMiddleLeft.x = overlayPoint.x;
    topMiddleLeft.y = overlayPoint.y;
    points.push_back(topMiddleLeft);

    Point topMiddleRight;
    topMiddleRight.x = stripePlotWidth - stripePlotHorizontalPadding;
    topMiddleRight.y = overlayPoint.y;
    points.push_back(topMiddleRight);

    Point topRight;
    topRight.x = topMiddleRight.x + cutoffTriangleHorizontalOffset;
    topRight.y = topMiddleRight.y - cutoffTriangleVerticalOffset;
    points.push_back(topRight);

    Point bottomMiddleRight;
    bottomMiddleRight.x = stripePlotWidth - stripePlotHorizontalPadding;
    bottomMiddleRight.y = previousYCoordinate;//overlayPoint.y+stripeThickness;    

    Point bottomRight;
    bottomRight.x = bottomMiddleRight.x + cutoffTriangleHorizontalOffset;
    bottomRight.y = bottomMiddleRight.y + cutoffTriangleVerticalOffset;
    points.push_back(bottomRight);

    points.push_back(bottomMiddleRight);

    points.push_back(topMiddleRight);

    points.push_back(topMiddleLeft);

    Point bottomMiddleLeft;
    bottomMiddleLeft.x = overlayPoint.x;
    bottomMiddleLeft.y = previousYCoordinate;//overlayPoint.y+stripeThickness;

    points.push_back(bottomMiddleLeft);

    points.push_back(bottomMiddleRight);

    points.push_back(bottomMiddleLeft);

    Point bottomLeft;
    bottomLeft.x = overlayPoint.x-cutoffTriangleHorizontalOffset;
    bottomLeft.y = previousYCoordinate+cutoffTriangleVerticalOffset;//overlayPoint.y+cutoffTriangleVerticalOffset+stripeThickness;
    points.push_back(bottomLeft);    

    Attributes polygonAttributes;

    polygonAttributes["stroke"] = "#000000";
    polygonAttributes["stroke-width"] = "2";

    svgPlot->drawPolygon(overlayStyle, points);

    Rectangle rectangle;
    rectangle.point = topMiddleLeft;
    rectangle.width = topMiddleRight.x - topMiddleLeft.x;
    rectangle.height = bottomMiddleLeft.y - topMiddleLeft.y;
    Attributes rectangleAttributes;
    rectangleAttributes["opacity"] = "0.0";

    Style rectangleStyle;
    rectangleStyle["fill"] = "#000000";
    svgPlot->drawRectangle(rectangleStyle, rectangleAttributes, rectangle);

    Point labelPoint;

    float labelPadding = 3.;
    labelPoint.x = bottomLeft.x - labelPadding;
    labelPoint.y = 0.5*(overlayPoint.y + previousYCoordinate + 14.); //  (bottomMiddleLeft.y + topMiddleLeft.y) / 2.;

    Style labelStyle; 
    Attributes labelAttributes; 
    labelAttributes["baseline"] = "middle";
    labelAttributes["font-family"] = "roman";
    labelAttributes["font-size"] = "14pt";
    labelAttributes["text-anchor"] = "end";

    stringstream labelText;

    labelText << showpoint << setprecision(2) << cutoff;
    
    svgPlot->drawText(labelPoint, labelStyle, labelAttributes, labelText.str());

    svgPlot->closeTag("g");

    previousYCoordinate = overlayPoint.y;
    previousCutoff = cutoff;
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Save the python script required to adjust the cutoff in PyMol
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::savePymolDilutionScript(string prefix) {
  string pythonDilutionFilename = prefix + "_pymol.py";
  ofstream pythonDilutionFile(pythonDilutionFilename.c_str());

  pythonDilutionFile << "def sortByIntervalLength(firstInterval, secondInterval) : " << std::endl;
  pythonDilutionFile << "  firstIntervalLength = firstInterval[1] - firstInterval[0] " << std::endl;
  pythonDilutionFile << "  secondIntervalLength = secondInterval[1] - secondInterval[0] " << std::endl;
  pythonDilutionFile << "  return secondIntervalLength - firstIntervalLength" << std::endl;

  pythonDilutionFile << "import copy" << std::endl;
  //	pythonDilutionFile << "def setCutoff(desiredCutoff, minimumNumberOfSitesPerCluster = " << parameters.min_output_cluster_size << ", sortByClusterSize = False) : " << std::endl;

  pythonDilutionFile << "mapFromRigidLabelToInterval = {}" << std::endl;	// FIXME - global variable

  pythonDilutionFile << "def intervalSize(interval) :" << std::endl;
  pythonDilutionFile << "  return interval[1]-interval[0]" << std::endl;

  pythonDilutionFile << "def compareClusters(a, b) :" << std::endl;
  pythonDilutionFile << "  global mapFromRigidLabelToInterval" << std::endl;
  pythonDilutionFile << "  clusterA = mapFromRigidLabelToInterval[a]" << std::endl;
  pythonDilutionFile << "  clusterAsize = intervalSize(clusterA)" << std::endl;

  pythonDilutionFile << "  clusterB = mapFromRigidLabelToInterval[b]" << std::endl;
  pythonDilutionFile << "  clusterBsize = intervalSize(clusterB)" << std::endl;

  pythonDilutionFile << "  return clusterBsize - clusterAsize" << std::endl;

  pythonDilutionFile << "def renderDecomposition(representation = \"sticks\", minimumNumberOfSitesPerCluster=" << parameters.minimumRigidClusterSize << ", sortByClusterSize=True) : " << std::endl;
  pythonDilutionFile << "  global mapFromRigidLabelToInterval" << std::endl;
  pythonDilutionFile << "  rigidLabels = mapFromRigidLabelToInterval.keys()" << std::endl;

  pythonDilutionFile << "  if sortByClusterSize :" << std::endl;
  pythonDilutionFile << "    rigidLabels.sort(compareClusters)" << std::endl;

  pythonDilutionFile << "  clusterNumber = 0" << std::endl;
  pythonDilutionFile << "  for rigidLabel in rigidLabels:" << std::endl;
  pythonDilutionFile << "    if (mapFromRigidLabelToColor.has_key(rigidLabel)) :" << std::endl;
  pythonDilutionFile << "      minValue = float(mapFromRigidLabelToInterval[rigidLabel][0])/100-0.001" << std::endl;
  pythonDilutionFile << "      maxValue = float(mapFromRigidLabelToInterval[rigidLabel][1])/100" << std::endl;
  pythonDilutionFile << "      numberOfSitesInInterval = int(100*(maxValue - minValue))" << std::endl;
  pythonDilutionFile << "      if (numberOfSitesInInterval < minimumNumberOfSitesPerCluster):" << std::endl;
  pythonDilutionFile << "        break" << std::endl;
  pythonDilutionFile << "      clusterName = \"RC_%03d\"%clusterNumber" << std::endl;
  pythonDilutionFile << "      clusterColor = mapFromRigidLabelToColor[rigidLabel]" << std::endl;
  pythonDilutionFile << "      cmd.select(clusterName, \"b > %f and b < %f\" % (minValue, maxValue))" << std::endl;
  pythonDilutionFile << "      cmd.as(representation, clusterName)" << std::endl;
  pythonDilutionFile << "    else : " << std::endl;
  pythonDilutionFile << "      clusterColor = \"grey\"" << std::endl;
  pythonDilutionFile << "    cmd.color(clusterColor, clusterName)" << std::endl;
  pythonDilutionFile << "    clusterNumber += 1" << std::endl;
  pythonDilutionFile << "  cmd.deselect()" << std::endl;		

  pythonDilutionFile << "def magicDecomposition(criticalSize = 5, sortByClusterSize = True, representation = \"sticks\") : " << std::endl;
  pythonDilutionFile << "  cmd.delete('RC_*')" << std::endl;
  pythonDilutionFile << "  cmd.refresh()" << std::endl;
  pythonDilutionFile << "  cmd.select('all')" << std::endl;
  pythonDilutionFile << "  cmd.as('lines')" << std::endl;
  pythonDilutionFile << "  cmd.color('grey')" << std::endl;

  pythonDilutionFile << "  global mapFromRigidLabelToInterval" << std::endl;		
  pythonDilutionFile << "  mapFromRigidLabelToInterval = copy.deepcopy(mapFromRigidLabelToLargestInterval)" << std::endl;	// FIXME - start with initial set
  pythonDilutionFile << "  mapFromRigidLabelToPutativeInterval = copy.deepcopy(mapFromRigidLabelToLargestInterval)" << std::endl;	// FIXME - start with initial set
  pythonDilutionFile << "  cutoffs = mapFromCutoffToNewRigidClusterIntervals.keys()" << std::endl;
  pythonDilutionFile << "  cutoffs.sort()" << std::endl; // TODO - the insert order was already reversed - is this necessary? (it seems so :-( )
  pythonDilutionFile << "  cutoffs.reverse()" << std::endl;
  pythonDilutionFile << "  for cutoff in cutoffs :" << std::endl;		
  pythonDilutionFile << "    newRigidClusterIntervals = mapFromCutoffToNewRigidClusterIntervals[cutoff] " << std::endl;
  pythonDilutionFile << "    for addedRigidLabel in newRigidClusterIntervals.keys():" << std::endl;
  pythonDilutionFile << "      newInterval = newRigidClusterIntervals[addedRigidLabel]" << std::endl;
  pythonDilutionFile << "      oldInterval = newInterval" << std::endl;
  pythonDilutionFile << "      if addedRigidLabel in mapFromRigidLabelToInterval :" << std::endl;
  pythonDilutionFile << "        oldInterval = mapFromRigidLabelToPutativeInterval[addedRigidLabel]" << std::endl;
  pythonDilutionFile << "      newInterval[1] = oldInterval[1]" << std::endl;
  pythonDilutionFile << "      mapFromRigidLabelToPutativeInterval[addedRigidLabel] = newInterval" << std::endl;
  pythonDilutionFile << "      if (intervalSize(newInterval) > criticalSize) :" << std::endl;
  pythonDilutionFile << "        mapFromRigidLabelToInterval[addedRigidLabel] = newInterval" << std::endl;
  pythonDilutionFile << "  renderDecomposition(representation, 0, sortByClusterSize)" << std::endl;

  pythonDilutionFile << "def setCutoff(desiredCutoff, minimumNumberOfSitesPerCluster = " << parameters.minimumRigidClusterSize << ", sortByClusterSize = True, representation = \"sticks\") : " << std::endl;
  pythonDilutionFile << "  cmd.delete('RC_*')" << std::endl;
  pythonDilutionFile << "  cmd.refresh()" << std::endl;
  pythonDilutionFile << "  cmd.select('all')" << std::endl;
  pythonDilutionFile << "  cmd.as('lines')" << std::endl;
  pythonDilutionFile << "  cmd.color('grey')" << std::endl;

  pythonDilutionFile << "  global mapFromRigidLabelToInterval" << std::endl;
  pythonDilutionFile << "  mapFromRigidLabelToInterval = copy.deepcopy(mapFromRigidLabelToLargestInterval)" << std::endl;	// FIXME - start with initial set
  pythonDilutionFile << "  cutoffs = mapFromCutoffToNewRigidClusterIntervals.keys()" << std::endl;
  pythonDilutionFile << "  cutoffs.sort()" << std::endl; // TODO - the insert order was already reversed - is this necessary? (it seems so :-( )
  pythonDilutionFile << "  cutoffs.reverse()" << std::endl;
  pythonDilutionFile << "  for cutoff in cutoffs :" << std::endl;
  pythonDilutionFile << "    if cutoff <= desiredCutoff :" << std::endl;  // TODO - verify boundary conditions
  pythonDilutionFile << "      break" << std::endl;

  pythonDilutionFile << "    newRigidClusterIntervals = mapFromCutoffToNewRigidClusterIntervals[cutoff] " << std::endl;
  pythonDilutionFile << "    for addedRigidLabel in newRigidClusterIntervals.keys():" << std::endl;
  pythonDilutionFile << "      mapFromRigidLabelToInterval[addedRigidLabel] = newRigidClusterIntervals[addedRigidLabel]" << std::endl;

  pythonDilutionFile << "  renderDecomposition(representation, minimumNumberOfSitesPerCluster, sortByClusterSize)" << std::endl;

  pythonDilutionFile.close();

  string pymolDilutionFilename = prefix + "_pymol.pml";
  ofstream pymolDilutionFile(pymolDilutionFilename.c_str());

  string dilutionIntervalsFilename = prefix + "_intervals_pymol";
  pymolDilutionFile << "cmd.do('run " << dilutionIntervalsFilename << "');" << std::endl;

  pymolDilutionFile << "mapFromRigidLabelToColor = {}" << std::endl;

  // construct a pymol dictionary from rigid labels to color
  pymolDilutionFile << "mapFromRigidLabelToColor = {}" << std::endl;
  for (std::map<RigidLabel,  ColorNumber>::iterator rigidLabelToColorIterator = rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber.begin();
      rigidLabelToColorIterator != rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber.end();
      rigidLabelToColorIterator++) {
    RigidLabel rigidLabel = (*rigidLabelToColorIterator).first;
    pymolDilutionFile << "mapFromRigidLabelToColor[" << rigidLabel << "] = '" << Color::getPymolColor(rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[rigidLabel]) << "';" << std::endl;
  }

  // Python / PyMol script for showing rigid cluster decomposition 
  // FIXME - refactor this out of here
  string pdbDilutionFilename = prefix + ".pdb";

  pymolDilutionFile << "cmd.load('"+pdbDilutionFilename+"');" << std::endl;
  pymolDilutionFile << "cmd.zoom();" << std::endl;
  pymolDilutionFile << "cmd.do('run " << pythonDilutionFilename << "');" << std::endl;

  pymolDilutionFile.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Save the statistics of a dilution hierarchy
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveDilutionStatistics(string prefix) {
  string dilutionStatisticsFilename = prefix + ".statistics";

  ofstream dilutionStatisticsFile(dilutionStatisticsFilename.c_str());

  dilutionStatisticsFile << "cutoff \t number of rigid clusters per total number of sites\t mean Rigid Cluster Size per total number of sites" << std::endl; // fraction of protein / atoms / etc

  for (std::map<float, CollectionOfRigidClusters >::reverse_iterator cutoffToCollectionOfLargestRigidClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rbegin();
      cutoffToCollectionOfLargestRigidClusters != rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rend();
      cutoffToCollectionOfLargestRigidClusters++) {
    std::pair<float, CollectionOfRigidClusters > pair = (*cutoffToCollectionOfLargestRigidClusters);

    float cutoff = pair.first;

    CollectionOfRigidClusters setOfRigidClusters = pair.second;

    float meanRigidClusterSize = computeMeanRigidClusterSize(setOfRigidClusters);
    float largestRigidClusterSize = computeLargestRigidClusterSize(setOfRigidClusters);

    float totalSites = structure->total_sites;

    dilutionStatisticsFile << cutoff << "\t" << ((float)setOfRigidClusters.size()/totalSites) << "\t" << (meanRigidClusterSize/(totalSites)) << "\t" << (largestRigidClusterSize/(totalSites)) << std::endl; // FIXME - verify that this should be (meanRigidClusterSize/(structure->total_sites+1) and not (meanRigidClusterSize/(structure->total_sites)
  }

  dilutionStatisticsFile.close();

}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Generates the jmol script required to adjust the cutoff in a 
// molecular / 3D visualization tile
//
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveJmolDilutionScript(string prefix) {
  // FIXME - this was copied / pasted and translated from PyMol - would be *VERY* nice to have a generic adaptor so this would be done once and not violate DRY 
  // FIXME - continue translating

  string rcdDilutionFilename = prefix + ".pdb";

  string jmolDilutionFilename = prefix + ".js";
  ofstream jmolDilutionFile(jmolDilutionFilename.c_str());

  jmolDilutionFile << "cutoffNumber = 0;" << std::endl;
  // it would be nice to properly scope cutoffs but it will require a bit of effort to keep working with both javascript and python 
  //  jmolDilutionFile << "viewProperties.listOfCutoffs = cutoffs;" << std::endl;

  jmolDilutionFile << "mapFromCutoffToNewRigidClusterIntervals = mapFromCutoffToNewRigidClusterIntervals;" << std::endl;
  jmolDilutionFile << "mapFromRigidLabelToLargestInterval = mapFromRigidLabelToLargestInterval;" << std::endl;

  jmolDilutionFile << "viewProperties.colorByRigidClusterJmolScript = '';" << std::endl;
  jmolDilutionFile << "viewProperties.showWireframeRigidClusterJmolScript = '';" << std::endl;
  jmolDilutionFile << "viewProperties.showBackbone = true;" << std::endl;
  jmolDilutionFile << "viewProperties.showOverlay = false;" << std::endl;
  jmolDilutionFile << "viewProperties.showNetwork = true;" << std::endl;
  jmolDilutionFile << "viewProperties.showAnimation = false;" << std::endl;
  jmolDilutionFile << "viewProperties.modelURI = '" << rcdDilutionFilename  << "';" << std::endl; // FIXME

  // the following requires dilutionIntervalsFilename to be included 

  // construct a pymol dictionary from rigid labels to color
  jmolDilutionFile << "mapFromRigidLabelToColor = {};" << std::endl;
  for (std::map<unsigned int, unsigned int>::iterator rigidLabelToColorIterator = rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber.begin();
      rigidLabelToColorIterator != rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber.end();
      rigidLabelToColorIterator++) {
    RigidLabel rigidLabel = (*rigidLabelToColorIterator).first;
    jmolDilutionFile << "mapFromRigidLabelToColor[" << rigidLabel << "] = '" << Color::getJmolColor(rigidClusterHierarchy->mapFromRigidClusterLabelToColorNumber[rigidLabel]) << "';" << std::endl;
  }

  // Python / PyMol script for showing rigid cluster decomposition 
  // FIXME - refactor this out of here
  jmolDilutionFile << "sortByIntervalLength = function(firstInterval, secondInterval) { " << std::endl;
  jmolDilutionFile << "  firstIntervalLength = firstInterval[1] - firstInterval[0];" << std::endl;
  jmolDilutionFile << "  secondIntervalLength = secondInterval[1] - secondInterval[0];" << std::endl;
  jmolDilutionFile << "  return secondIntervalLength - firstIntervalLength;" << std::endl;
  jmolDilutionFile << "}" << std::endl;

  jmolDilutionFile << "setCutoff = function(desiredCutoff, minimumNumberOfSitesPerCluster) { " << std::endl;
  //	jmolDilutionFile << "  if (viewProperties.getCutoff() < desiredCutoff) {" << std::endl;
  //	jmolDilutionFile << "    alert(' ' + viewProperties.getCutoff() + ' -> ' +desiredCutoff + ' ' );" << std::endl;
  jmolDilutionFile << "    viewProperties.showWireframeRigidClusterJmolScript = 'select all; wireframe .01; spacefill off;';" << std::endl;
  jmolDilutionFile << "    viewProperties.colorByRigidClusterJmolScript = 'select all; color grey; ';" << std::endl;
  //	jmolDilutionFile << "  }" << std::endl;
  jmolDilutionFile << "  var mapFromRigidLabelToInterval = {};" << std::endl; // FIXME - make a deep copy; this is a bad kludge
  jmolDilutionFile << "  mapFromRigidLabelToInterval[0] = this.mapFromRigidLabelToLargestInterval[0];" << std::endl; // FIXME - make a deep copy; this is a bad kludge

  jmolDilutionFile << "  cutoffNumber = 0;" << std::endl;
  jmolDilutionFile << "  for (var cutoffIndex in cutoffs) {" << std::endl;
  jmolDilutionFile << "    cutoffNumber = cutoffIndex;" << std::endl; // TODO - verify the boundary condition 
  jmolDilutionFile << "    cutoff = cutoffs[cutoffIndex];" << std::endl; // TODO - verify the boundary condition 

  jmolDilutionFile << "    if (cutoff < desiredCutoff) {" << std::endl; // TODO - verify the boundary condition 
  jmolDilutionFile << "        break;" << std::endl;
  jmolDilutionFile << "    }" << std::endl;

  jmolDilutionFile << "    newRigidClusterIntervals = this.mapFromCutoffToNewRigidClusterIntervals[cutoff]; " << std::endl;
  jmolDilutionFile << "    for (addedRigidLabel in newRigidClusterIntervals) {" << std::endl;
  jmolDilutionFile << "      mapFromRigidLabelToInterval[addedRigidLabel] = newRigidClusterIntervals[addedRigidLabel];" << std::endl;
  jmolDilutionFile << "    }" << std::endl;
  jmolDilutionFile << "  }" << std::endl;

  jmolDilutionFile << "  viewProperties.setCutoffNumber(cutoffNumber);" << std::endl;
  jmolDilutionFile << "  clusterNumber = 0" << std::endl;
  jmolDilutionFile << "  for (rigidLabel in mapFromRigidLabelToInterval) {" << std::endl;
  jmolDilutionFile << "    minValue = (mapFromRigidLabelToInterval[rigidLabel][0])/100-0.001;" << std::endl;
  jmolDilutionFile << "    maxValue = (mapFromRigidLabelToInterval[rigidLabel][1])/100;" << std::endl;
  jmolDilutionFile << "    numberOfSitesInInterval = mapFromRigidLabelToInterval[rigidLabel][1] - mapFromRigidLabelToInterval[rigidLabel][0];" << std::endl;
  jmolDilutionFile << "    if (numberOfSitesInInterval < minimumNumberOfSitesPerCluster) {" << std::endl;
  jmolDilutionFile << "      break;" << std::endl;
  jmolDilutionFile << "    }" << std::endl;
  jmolDilutionFile << "    var clusterColor = 'grey';" << std::endl;
  jmolDilutionFile << "    try {" << std::endl;
  jmolDilutionFile << "      clusterColor = this.mapFromRigidLabelToColor[rigidLabel];" << std::endl;
  jmolDilutionFile << "    } catch (exception) {}; // do nothing for now (just use grey)" << std::endl;
  jmolDilutionFile << "    viewProperties.colorByRigidClusterJmolScript += 'select temperature > ' + minValue +' and temperature < ' + maxValue +';';" << std::endl;
  jmolDilutionFile << "    viewProperties.showWireframeRigidClusterJmolScript += 'select temperature > ' + minValue +' and temperature < ' + maxValue +';';" << std::endl;

  jmolDilutionFile << "    viewProperties.colorByRigidClusterJmolScript += 'color '+ this.mapFromRigidLabelToColor[rigidLabel] + ';';" << std::endl;
  jmolDilutionFile << "    viewProperties.showWireframeRigidClusterJmolScript += 'wireframe .2;';" << std::endl;

  //	jmolDilutionFile << "      cmd.as(\"sticks\", clusterName)" << std::endl;
  jmolDilutionFile << "    }" << std::endl;
  jmolDilutionFile << "    clusterNumber += 1;" << std::endl;
  jmolDilutionFile << "  viewProperties.updateView();" << std::endl;

  jmolDilutionFile << "  }" << std::endl;

  jmolDilutionFile.close();
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the rigid label of siteNumber at cutoff
// 
////////////////////////////////////////////////////////////////////////////////
unsigned int RigidClusterAnalysis::getRigidLabel(SiteID siteNumber, 
    Cutoff cutoff) {

  std::pair<SiteID, Cutoff> methodParameters(siteNumber, cutoff);

  static std::map<std::pair<SiteID, Cutoff>, unsigned int> cache;
  bool cached = true;
  if (cached) {
    if (cache.find(methodParameters) != cache.end()) {
      return cache[methodParameters];
    }
  }

  RigidCluster * rigidCluster = rigidClusterHierarchy->mapFromSiteIDToRigidClusterContainingSiteID[siteNumber];

  float workingCutoff = getClosestCutoffBelowCutoff(cutoff);
  RigidLabel rigidLabel = getRigidLabel(rigidCluster,  workingCutoff);

  if (cached) {
    cache[methodParameters] = rigidLabel;
  }

  return rigidLabel;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the rigid label of rigidCluster at cutoff 
//
////////////////////////////////////////////////////////////////////////////////
unsigned int RigidClusterAnalysis::getRigidLabel(RigidCluster * rigidCluster, 
    Cutoff cutoff) {

  unsigned int rigidLabel = UINT_MAX;

  if (rigidCluster == NULL) {

    return UINT_MAX; 
  }

  std::pair<RigidCluster *, Cutoff> methodParameters(rigidCluster, cutoff);

  static std::map<std::pair<RigidCluster *, Cutoff>, unsigned int> cache;
  bool cached = true;
  if (cached) {
    if (cache.find(methodParameters) != cache.end()) {
      return cache[methodParameters];
    }
  }

  rigidLabel = rigidCluster->minLabel;

  // ensure that there is a parent rigid cluster to rigidCluster
  if (rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster.find(rigidCluster) != rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster.end()) {
    RigidCluster * parentRigidCluster = rigidClusterHierarchy->mapFromRigidSubClusterToParentRigidCluster[rigidCluster];

    Cutoff parentRigidClusterCutoff = rigidClusterHierarchy->mapFromRigidClusterToCutoff[parentRigidCluster];

    // continue walking up the rigid cluster hierarchy until we find a rigid cluster for a cutoff that exceeds the parameter cutoff
    if (parentRigidClusterCutoff <= cutoff) {
      rigidLabel = getRigidLabel(parentRigidCluster, cutoff);
    } /*else if (rigidCluster->size() < parameters.minimumRigidClusterSize) { 

      rigidLabel = UINT_MAX; 
      // TODO - would be nice to have a better response ; perhaps a RigidLabel type?
      // also TODO - do we want to label single-site rigid clusters?		
    }		*/
  }

  if (cached) {
    cache[methodParameters] = rigidLabel;
  }

  return rigidLabel;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getTicSpacing(float intervalLength, unsigned int numberOfTics) {
  float ticSpacing = 10.0;

  double idealTicSpacing = (double)intervalLength / (double) numberOfTics;

  double powerOfTen = floor(log(idealTicSpacing) / log(10.0));

  float leadingDigit = idealTicSpacing / pow(10.0, powerOfTen);

  // prefer to use 1, 2, 5, 10, 20, 50, 100, etc... 
  if (leadingDigit < 1.5) {
    leadingDigit = 1.0;
  } else if (leadingDigit < 3.5) {
    leadingDigit = 2.0;
  } else if (leadingDigit < 6.5) {
    leadingDigit = 5.0;
  } else {
    leadingDigit = 10.0;
  }

  ticSpacing = leadingDigit * pow(10.0, powerOfTen);

  return ticSpacing;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the leftmost x coordinate of siteNumber in the rigid cluster dilution plot 
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getXCoordinate(SiteID siteNumber) {
  Site_Info siteInfo = structure->site_info[siteNumber];

  float xCoordinate = horizontalOffsetForPlotArea;

  switch (xAxisType) {
    case perSite: {
                    xCoordinate +=  widthPerSite * (float)siteNumber;

                  }	break;

    case perRigidLabel: {
                          unsigned int rigidLabel = rigidClusterHierarchy->mapFromSiteIdToRigidLabel[siteNumber]; // FIXME 

                          xCoordinate +=  widthPerSite * (float)rigidLabel;

                        } break;

    case perResidue: {
                       if (structure->mapFromSiteIDToResidue.find(siteNumber) != structure->mapFromSiteIDToResidue.end()) {
                         Residue *residue = structure->mapFromSiteIDToResidue[siteNumber];
                         xCoordinate = getXCoordinate(*residue);

                         size_t residueSize = residue->size();

                         float widthPerSiteInResidue = widthPerResidue / (float)residueSize;

                         unsigned int siteNumberInResidue = 0;
                         for (Residue::Sites::iterator siteIterator = residue->beginSites();
                             siteIterator != residue->endSites();
                             siteIterator ++) {
                           Site_Info * site = *siteIterator;
                           if (site->site_number == siteNumber) {
                             break;
                           }

                           siteNumberInResidue++;
                         }

                         xCoordinate += (float) siteNumberInResidue * widthPerSiteInResidue; 
                       }

                     } break;
  }

  return xCoordinate;

}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the leftmost x coordinate or residue in the rigid cluster dilution plot 
//
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getXCoordinate(Residue & residue){
  float xCoordinate = horizontalOffsetForPlotArea;

  switch (xAxisType) {
    case perSite:{
                   Site_Info *firstResidueSiteInfo = *(residue.beginSites());

                   SiteID firstResidueSiteID = firstResidueSiteInfo->site_number;

                   return getXCoordinate(firstResidueSiteID);
                 } break;

    case perRigidLabel: {
                          Site_Info *firstResidueSiteInfo = *(residue.beginSites());

                          SiteID firstResidueSiteID = firstResidueSiteInfo->site_number;

                          return getXCoordinate(firstResidueSiteID);
                        } break;

    case perResidue:
                        Chain *chain = (structure->mapFromResidueToChain)[&residue];
                        Chain::iterator residueIterator = chain->begin();

                        Residue * firstResidueInChain = residueIterator->second;
                        unsigned int firstResidueIdInChain = structure->mapFromResidueToResidueID[firstResidueInChain];
                        xCoordinate = getXCoordinate(*chain);  // FIXME - why does this not properly offset the xCoordinate for the chain?

                        // get residue number in chain
                        unsigned int residueID = structure->mapFromResidueToResidueID[&residue];
                        unsigned int residueIDinChain = residueID - firstResidueIdInChain;
                        xCoordinate += widthPerResidue*(float)(residueIDinChain);

                        break;
  }

  return xCoordinate;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the leftmost x coordinate of site in the rigid cluster dilution plot 
//
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getXCoordinate(Site_Info & site) {
  float xCoordinate = getXCoordinate(site.site_number);

  return xCoordinate;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the leftmost x coordinate of chain in the rigid cluster dilution plot
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getXCoordinate(Chain & chain) {
  float xCoordinate = horizontalOffsetForPlotArea;

  switch (xAxisType) {
    case perSite:

      for (std::map <unsigned int, Chain*>::iterator chainIterator = structure->mapFromChainIDtoChain.begin();
          chainIterator->second != &chain;
          chainIterator++) {
        Chain *previousChain = chainIterator->second;

        for (Chain::iterator residueIterator = previousChain->begin();
            residueIterator != previousChain->end();
            residueIterator++) {
          Residue * residue = residueIterator->second;
          xCoordinate += widthPerSite*(float)residue->size(); 
        }

        xCoordinate += offsetBetweenChains;
      }

      // TODO - implement
      break;

    case perResidue:
      for (std::map <unsigned int, Chain*>::iterator chainIterator = structure->mapFromChainIDtoChain.begin();
          chainIterator->second != &chain;
          chainIterator++) {
        Chain *previousChain = chainIterator->second;
        xCoordinate += widthPerResidue*(float)previousChain->size(); 
        xCoordinate += offsetBetweenChains;
      }

      return xCoordinate;

      break;

    case perRigidLabel: 
      // TODO - implement (does this even make sense for a chain?)
      break;
  }

  return 0.;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the width of siteNumber in the rigid cluster dilution plot 
//
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getWidth(unsigned int siteNumber) {
  Site_Info siteInfo = structure->site_info[siteNumber];

  float width = 0.;

  switch (xAxisType) {
    case perSite:
      width = horizontalOffsetForPlotArea + widthPerSite * (float)siteNumber;

      break;

    case perResidue:
      // get residue for site
      //			Residue *residue = ;

      width = horizontalOffsetForPlotArea + widthPerSite * (float)siteNumber; // WTF? 
      break;

    case perRigidLabel: 
      width = widthPerSite;
      break;
  }

  if ((xAxisType == perSite) || (xAxisType == perSite)){
    width += siteInfo.FIRST_chain_ID;
  }

  return width;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// the width of residue in the rigid cluster dilution plot
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getWidth(Residue & residue){
  float width = 0.;

  switch (xAxisType) {
    case perSite:
      width += widthPerSite*(float)residue.size();

      break;

    case perResidue:
      width += widthPerResidue;

      break;

    case perRigidLabel: 
      width += widthPerSite;
      break;
  }

  return width;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns: 
// the width of site in the rigid cluster dilution plot
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getWidth(Site_Info & site) {
  switch (xAxisType) {
    case perSite:
      return widthPerSite;

      break;

    case perResidue:
      return widthPerResidue;
      break;

    case perRigidLabel: 
      return widthPerSite;
      break;
  }

  return 0.;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns: width of the chain in the plot
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getWidth(Chain & chain) {
  size_t chainSize = chain.size();

  switch (xAxisType) {
    case perSite:
      return widthPerSite;
      break;

    case perRigidLabel:
      return widthPerSite;
      break;

    case perResidue:
      float width = widthPerResidue*(float)chainSize; 
      return width;
      break;
  }

  return 0.;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns: The y coordinate in the dilution plot coorsponding to cutoff
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getYCoordinateFromCutoff(Cutoff cutoff) {	
  float yCoordinate = 0.0;

  if (parameters.useLinearCutoffScale) {
    float minimumCutoff = *distinctCutoffs.begin();
    float maximumCutoff = *distinctCutoffs.rbegin();

    float fractionalYcoordinate = (cutoff - minimumCutoff) / (maximumCutoff - minimumCutoff);
    yCoordinate = verticalOffsetForPlotArea + (stripePlotHeight - verticalOffsetForPlotArea) * (1.0 - fractionalYcoordinate);

  } else {
    std::set<Cutoff>::iterator cutoffIterator = distinctCutoffs.lower_bound(cutoff);
    std::set<Cutoff> setOfCutoffsBelowCurrent;
    setOfCutoffsBelowCurrent.insert(distinctCutoffs.begin(), cutoffIterator); // FIXME - compute the number below cutoffIterator without the set insertion (not the most efficient approach)

    size_t cutoffsBelowCurrent = setOfCutoffsBelowCurrent.size();

    yCoordinate = stripePlotHeight - (stripeThickness + separationBetweenStripes) * (float) cutoffsBelowCurrent;
  }

  return yCoordinate;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Save the set of rigid label intervals corresponding to the rigid cluster decomposition
// at a given cutoff. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveRigidClusterIntervals(ofstream & outputStream,
    float cutoff,
    CollectionOfRigidClusters setOfRigidClusters) {
  outputStream << "{";

  CollectionOfRigidClusters::iterator rigidClusterIterator = setOfRigidClusters.begin();
  while (rigidClusterIterator != setOfRigidClusters.end()) {

    RigidCluster *rigidCluster = *rigidClusterIterator;
    outputStream << getRigidLabel(rigidCluster, cutoff)  << ": ";
    outputStream << "[" << rigidCluster->minLabel << ", " << rigidCluster->maxLabel << "]";
    rigidClusterIterator++;

    if (rigidClusterIterator != setOfRigidClusters.end()) {
      outputStream << ", ";
    }

  }

  outputStream << "}";
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// 
// Returns: 
// The arithmetic mean number of sites included in a rigid cluster 
// at a given cutoff
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::computeMeanRigidClusterSize(float cutoff) {
  float workingCutoff = getClosestCutoffBelowCutoff(cutoff);

  CollectionOfRigidClusters setOfRigidClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters[workingCutoff];

  return computeMeanRigidClusterSize(setOfRigidClusters);
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// 
// Returns: 
// The largest number of sites included in a rigid cluster at a given cutoff
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::computeLargestRigidClusterSize(float cutoff) {
  float workingCutoff = getClosestCutoffBelowCutoff(cutoff);

  CollectionOfRigidClusters setOfRigidClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters[workingCutoff];

  return computeLargestRigidClusterSize(setOfRigidClusters);
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// The largest number of sites included in a rigid cluster from a 
// setOfRigidClusters
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::computeLargestRigidClusterSize(CollectionOfRigidClusters setOfRigidClusters) {
  size_t largestRigidClusterSize = 0;
  for (CollectionOfRigidClusters::iterator rigidCluster = setOfRigidClusters.begin();
      rigidCluster != setOfRigidClusters.end();
      rigidCluster++) {
    size_t currentClusterSize = (*rigidCluster)->size();
    if (largestRigidClusterSize < currentClusterSize) {
      largestRigidClusterSize = currentClusterSize;
    }

  }

  return largestRigidClusterSize;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// The arithmatic mean number of ssites included in a rigid cluster from a 
// setOfRigidClusters
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::computeMeanRigidClusterSize(CollectionOfRigidClusters setOfRigidClusters) {
  unsigned int sumOfRigidClusterSizes = 0;

  for (CollectionOfRigidClusters::iterator rigidClusterIterator = setOfRigidClusters.begin();
      rigidClusterIterator != setOfRigidClusters.end();
      rigidClusterIterator++) {
    RigidCluster * rigidCluster = *rigidClusterIterator;
    sumOfRigidClusterSizes += rigidCluster->size();
  }    

  float meanRigidClusterSize = (float)sumOfRigidClusterSizes/((float)setOfRigidClusters.size());

  return meanRigidClusterSize;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// 
////////////////////////////////////////////////////////////////////////////////
float RigidClusterAnalysis::getClosestCutoffBelowCutoff(Cutoff cutoff) {
  Cutoff closestCutoff = *rigidClusterHierarchy->setOfCutoffs.lower_bound(cutoff);

  return closestCutoff;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns: true if siteID should be included in the rigid cluster dilution plot
// false otherwise
// 
////////////////////////////////////////////////////////////////////////////////
bool RigidClusterAnalysis::includeSite(SiteID siteID) {
  return RigidClusterAnalysis::includeSite(structure, xAxisType, siteID);
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
//
// Returns:
// true if siteID should be included in the rigid cluster dilution plot 
// for structure with xAxisType
//
// false otherwise
// 
////////////////////////////////////////////////////////////////////////////////
bool RigidClusterAnalysis::includeSite(MolFramework *structure, 
    XAxisType xAxisType,
    SiteID siteID) {

  switch(xAxisType) {
    case perSite:
      return true;
      break;

    case perResidue:
      return  ((structure->isBackbone(siteID) == 1) ||
          (structure->is_nucleic_backbone(siteID) == 1) || 
          (structure->isBackbone(siteID) == 3) ||
          (structure->is_nucleic_backbone(siteID) == 3));
      /*      return  ((structure->isBackbone(siteID) == 1) ||
              (structure->is_nucleic_backbone(siteID) == 1));*/

      /*        return ((structure->isBackbone(siteID) != 0) ||
                (structure->is_nucleic_backbone(siteID) != 0));*/
      break;

    case perRigidLabel:
      return true; 
      break;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Description: 
// Generate a the script content for a dilution plot. 
// 
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::dilutionPlotScript(Plot *plot) {
  stringstream scriptContent;
  scriptContent << "var highlightOpacity = .9;" << std::endl;
  scriptContent << "var mouseoverOpacity = .5;" << std::endl;
  scriptContent << "var selectedCutoff = 0;" << std::endl;
  scriptContent << "var selectedTarget = null;" << std::endl;

  scriptContent << " function highlight(evt, cutoff) { " << std::endl;
  scriptContent << "	var target = document.getElementById(cutoff);  " << std::endl;
  scriptContent << "	target.setAttribute('opacity', mouseoverOpacity); " << std::endl;
  scriptContent << "} " << std::endl;

  scriptContent << "function unhighlight(evt, cutoff) { " << std::endl;
  scriptContent << "	var target = document.getElementById(cutoff);  " << std::endl;
  scriptContent << "	if (selectedCutoff != cutoff) { " << std::endl;
  scriptContent << "	  target.setAttribute('opacity', '0'); " << std::endl; // TODO - DRY 
  scriptContent << "  } else {" << std::endl;
  scriptContent << "    target.setAttribute('opacity', highlightOpacity); " << std::endl;
  scriptContent << "  }" << std::endl;
  scriptContent << "} " << std::endl;

  scriptContent << "function selectCutoff(evt, cutoff) { " << std::endl;
  scriptContent << "	if (selectedTarget != null) {" << std::endl;
  scriptContent << "	  selectedTarget.setAttribute('opacity', '0');" << std::endl; // TODO - DRY 
  scriptContent << "  } " << std::endl;
  scriptContent << "	var target = document.getElementById(cutoff);  " << std::endl;
  scriptContent << "	target.setAttribute('opacity', highlightOpacity); " << std::endl;
  scriptContent << "	selectedTarget = document.getElementById(cutoff);  " << std::endl;

  scriptContent << "	selectedCutoff = cutoff; " << std::endl;

  // FIXME - this is a fairly ugly kludge so that we can call the browser's javascript engine 
  // something like this is necessary to work cross-browser at present but it should be cleaned up 
  scriptContent << "	if (typeof browserEval != 'undefined') {" << std::endl;
  scriptContent << "	  browserEval('window.location=\\'javascript:setCutoff(' +cutoff+ ')\\'');" << std::endl;
  scriptContent << "	} else if (typeof top != 'undefined') {" << std::endl;
  scriptContent << "	  top.setCutoff(cutoff)" << std::endl;
  scriptContent << "	} else {" << std::endl;
  scriptContent << "	  window.location = 'javascript:setCutoff(' +cutoff+ ')'" << std::endl;
  scriptContent << "	}" << std::endl;
  // end FIXME

  scriptContent << "} " << std::endl;

  plot->script(scriptContent.str());
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   A function to export several additional hydrogen bond dilution data files.
//   Some of these are then used to define a vector describing the unfolding
//   pathway.
//
//   Written by AJ Rader IUPUI 03.30.06 ajrader@iupui.edu
////////////////////////////////////////////////////////////////////////////////
void RigidClusterAnalysis::saveDilutionData(string prefix) {

  string hbdiluteFilename = prefix + "_hbstats.txt"; 
  ofstream hbdiluteFile(hbdiluteFilename.c_str());
  string lrcDiluteFilename = prefix + "_lrc.txt"; 
  ofstream lrcDiluteFile(lrcDiluteFilename.c_str());
  string pathwayDiluteFilename = prefix + "_pathway.txt"; 
  ofstream pathwayDiluteFile(pathwayDiluteFilename.c_str());

  unsigned int dilutionSteps=0;
  float pathwayAvgRClust[structure->total_residues];
  for(unsigned int a=0;a<structure->total_residues;a++)  
    pathwayAvgRClust[a]=0;

  hbdiluteFile <<"#Step     Ecut     <r>    RCcount   Xflop     Xlrc"<<endl;
  //345678901234567890123456789012345678901234567890123456789012345678901234567890
  std::map<Cutoff, CollectionOfRigidClusters >::reverse_iterator cutoffToRClusters;
  for (cutoffToRClusters = rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rbegin();
      cutoffToRClusters != rigidClusterHierarchy->mapFromCutoffToSetOfLargestRigidClusters.rend();
      cutoffToRClusters++){

    std::pair<Cutoff, CollectionOfRigidClusters > pair = (*cutoffToRClusters);
    Cutoff cutoff = pair.first;

    CollectionOfRigidClusters setOfRigidClusters = pair.second;

    dilutionSteps++;
    unsigned int largestRigidClusterSize = 0;
    RigidLabel largestRigidClusterLabel = 0;


    //---------------- Find the Largest Rigid Cluster from the setOfRigidClusters
    for (CollectionOfRigidClusters::iterator rcIterator = setOfRigidClusters.begin();
        rcIterator != setOfRigidClusters.end(); rcIterator++) {

      RigidCluster *rigid_cluster = *rcIterator;
      RigidLabel rigidClusterLabel = getRigidLabel(rigid_cluster,cutoff);

      if(rigid_cluster->size() > largestRigidClusterSize) {
        largestRigidClusterSize = rigid_cluster->size();
        largestRigidClusterLabel = rigidClusterLabel;
      }
    }

    Sites::iterator siteIterator = sites.begin();
    //SiteID previousSiteID = *siteIterator;
    //Site_Info previousSite = structure->site_info[previousSiteID];
    unsigned int this_residue=0;
    int this_res_id=0;
    int *largestRCCount = new int[structure->total_residues];

    for(this_residue=0;this_residue<structure->total_residues;this_residue++) {
      largestRCCount[this_residue] = 0;
    }

    while(siteIterator !=sites.end()) {
      stringstream residueIDString;
      SiteID siteID = *siteIterator;
      Site_Info site = structure->site_info[siteID];    
      this_residue = site.seq_number; 
      residueIDString << site.seq_number << ";" << site.insert_code << ";"
        << site.chain_ID;
      // resIDString = site.seq_number+";"+site.insert_code+";"+site.chain_ID;
      this_res_id = structure->unique_res_id[residueIDString.str()];
      RigidLabel siteRigidClusterLabel = getRigidLabel(siteID, cutoff);

      if(siteRigidClusterLabel == largestRigidClusterLabel) {
        largestRCCount[this_res_id]++;
      }
      siteIterator++;
    }


    //lrcDiluteFile<<setw(10)<<setprecision(6)<<cutoff;
    lrcDiluteFile<<setw(4)<<dilutionSteps;
    for(this_residue=0;this_residue<structure->total_residues;this_residue++) {
      lrcDiluteFile<<setw(3)<<largestRCCount[this_residue];
      pathwayAvgRClust[this_residue] += (float) largestRCCount[this_residue]/2.0;
    }
    lrcDiluteFile<<endl;
    //float clusterDilutionParameters = calculateRigidClusterParameters(setOfRigidClusters);
    hbdiluteFile<<setw(4)<<dilutionSteps<<setw(11)<<setprecision(5)<<cutoff<<"  "<<structure->mean_coordination;
    hbdiluteFile<<setw(8)<< setOfRigidClusters.size()<<"   "<<endl;
    //		<<"\t"<<largestRCCount[this_res_id]<<endl;	    

  } 
  // std::map<string,int>::iterator resIterator=structure->unique_res_id.begin();
  //   for( resIterator = structure->unique_res_id.begin(); 
  //        resIterator != structure->unique_res_id.end(); resIterator++ ) {
  //     int resid = resIterator->second; 
  //     string residueIdentity = resIterator->first;
  //     pathwayAvgRClust[resid] = pathwayAvgRClust[resid]/((float) dilutionSteps);
  //     pathwayDiluteFile<<setw(10)<<residueIdentity<<"  "<<setw(11)<<setprecision(5)<<pathwayAvgRClust[resid]<<endl;
  //     //pathwayDiluteFile<<setw(6)<<residue<<"  "<<setw(11)<<setprecision(5)
  //     //		     <<pathwayAvgRClust[residue]<<endl;
  for(unsigned int residue=0; residue < structure->total_residues; residue++) {
    pathwayAvgRClust[residue] = pathwayAvgRClust[residue]/((float) dilutionSteps);
    pathwayDiluteFile<<setw(6)<<residue+1<<"  "<<setw(11)<<setprecision(5)
      <<pathwayAvgRClust[residue]<<endl;

  }

  hbdiluteFile.close();
  lrcDiluteFile.close();
  pathwayDiluteFile.close();
}
