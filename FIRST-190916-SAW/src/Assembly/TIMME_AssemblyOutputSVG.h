// Brandon Hespenheide (c) 2007
////////////////////////////////////////////////////////////////////////////////

#ifndef _TIMME_ASSEMBLY_OUTPUT_SVG_
#define _TIMME_ASSEMBLY_OUTPUT_SVG_

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Assembly.h"
#include "AssemblyStepIterator.h"
#include "XML_Writer.h"
#include "MolFramework.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
class TIMME_AssemblyOutputSVG : public AssemblyVisitor<I,T> {
	
private:
	MolFramework *molFramework;
	
	void saveStep(XML_Writer &xml_Writer, 
								I previousCutoff, I cutoff, 
								vector<Piece<I,T>*> &step);
	
	void setupDimensions(Assembly<I,T> *assembly);
	
	void makeTitle(XML_Writer &xml_Writer);
	void makeXtics(XML_Writer &xml_Writer);
	void makeYtics(XML_Writer &xml_Writer);
	
	string getWidthFromElement(T element);
	
	double getXCoordinateFromElement(T element);
	double getYCoordinateFromCutoff(I cutoff);
	
	// FIXME - refactor all of these out of here (bad separation of concerns)
	XAxisType xAxisType;
	I maxCutoff;
	I minCutoff;
	
	string yTicsWidth;
	string titleHeight;
	string dilutionPlotWidth;
	string dilutionPlotHeight;
	
	double cutoffHeightScale;
	double widthPerIndex;
	double separationBetweenChains;
	
	bool showElement(T element);
	unsigned int getXindex(T element);
	
public:
		TIMME_AssemblyOutputSVG(MolFramework *_molFramework ) :
    molFramework( _molFramework )
  {
			// setting defaults 
			// this should probably be factored out of here 
			xAxisType = perResidue; 
			yTicsWidth = "10em";
			titleHeight = "0";
			
			cutoffHeightScale = 1000.0;
			widthPerIndex = 0.4;
			separationBetweenChains = 5.0;
			
			dilutionPlotWidth = "10"; // FIXME - default to zero
			dilutionPlotHeight = "0";
	};
  ~TIMME_AssemblyOutputSVG()
  {};
	
public:
		virtual void visit( Assembly<I,T> *assembly );
	
	XAxisType getXAxisType();
	void setXAxisType(XAxisType);
	
};

template <class I, class T>
XAxisType TIMME_AssemblyOutputSVG<I,T>::getXAxisType() {
	return xAxisType;
}

template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::setXAxisType(XAxisType xAxisType) {
	this->xAxisType = xAxisType;
}

template <class I, class T>
double TIMME_AssemblyOutputSVG<I,T>::getXCoordinateFromElement(T element) {
	unsigned int index = getXindex(element);
	
	double xCoordinate = index * widthPerIndex;
	
	return xCoordinate;
}

template <class I, class T>
bool TIMME_AssemblyOutputSVG<I,T>::showElement(T element){
		switch(xAxisType) {
			case perSite: 
				return true;
				
				break;
				
			case perRigidLabel: 
				return true;
				
				break;
				
			case perResidue:
				return ((molFramework->isBackbone(element) == 1) ||
								 (molFramework->is_nucleic_backbone(element) == 1) || 
								 (molFramework->isBackbone(element) == 3) ||
								 (molFramework->is_nucleic_backbone(element) == 3));
				break;			
		};
	
	return true;
}

template <class I, class T>
unsigned int TIMME_AssemblyOutputSVG<I,T>::getXindex(T element) {
	unsigned int xIndex = 0;
		
	switch(xAxisType) {
		case perSite: 
			xIndex = element;
			
			break;
			
		case perRigidLabel: 
//			xIndex = getRigidLabel(element); // TODO - implement
			
			break;
			
		case perResidue:
			// xIndex = getResidue(element) + chainSeparatorIndexShift * getChain(element); // TODO - implement
			xIndex = element; // FIXME - removeme
			break;			
	};
	
	return xIndex;
}

template <class I, class T>
string TIMME_AssemblyOutputSVG<I,T>::getWidthFromElement(T element) {
	stringstream widthStringstream;
	
	widthStringstream << widthPerIndex;
/*	
	switch(xAxisType) {
		case perSite: 
			widthStringstream << widthPerSite;
			
			break;
			
		case perRigidLabel: 
			widthStringstream << widthPerSite;

			break;
			
		case perResidue:
			if (showElement(element)) {
				widthStringstream << widthPerResidue;

			} else {
				widthStringstream << "0";
				
			}
			
			break;			
		};
	*/
	return widthStringstream.str();
}

template <class I, class T>
double TIMME_AssemblyOutputSVG<I,T>::getYCoordinateFromCutoff(I cutoff) {
	double yCoordinate = (maxCutoff - cutoff)*cutoffHeightScale;
	
	return yCoordinate;
}

template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::makeTitle(XML_Writer &xml_Writer) {
	/*	xml_Writer.newTag("svg");
	xml_Writer.addAttribute("xmlns", "http://www.w3.org/2000/svg");
		
	xml_Writer.addAttribute("x", "0");
	xml_Writer.addAttribute("y", "0");
		
	// width = width of the y-axis label + width of the graph
	xml_Writer.addAttribute("width", "500"); 
	
	// height 
	xml_Writer.addAttribute("height", "10em");	
	xml_Writer.openTag();
	xml_Writer.newLine();
	
	xml_Writer.newTag("text");
	xml_Writer.addAttribute("font-family", "roman");
	xml_Writer.addAttribute("font-size", "24pt");
	xml_Writer.addAttribute("width", "100%");
	xml_Writer.addAttribute("alignment-baseline", "middle");
	
	xml_Writer.openTag();
	xml_Writer.newLine();
	xml_Writer.getStream() << "Title";
	xml_Writer.closeTag();
	xml_Writer.newLine();	
	
	xml_Writer.closeTag();
	xml_Writer.newLine();	*/
}

template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::makeXtics(XML_Writer &xml_Writer) {
	xml_Writer.newTag("svg");
	xml_Writer.addAttribute("xmlns", "http://www.w3.org/2000/svg");
		
	xml_Writer.addAttribute("x", yTicsWidth);
	xml_Writer.addAttribute("y", dilutionPlotWidth);
		
	xml_Writer.addAttribute("width", "dilutionPlotWidth");
	xml_Writer.addAttribute("height", "10em");	
	xml_Writer.openTag();
	xml_Writer.newLine();
	
	// TODO - iterate through appropriate labels and mark the tics accordingly
	switch(xAxisType) {
		case perSite: 
			/*xml_Writer.newTag("text");
			 xml_Writer.addAttribute("font-family", "roman");
			 xml_Writer.addAttribute("font-size", "24pt");
			 
			 xml_Writer.openTag();
			 xml_Writer.newLine();
			 xml_Writer.getStream() << "X Tics";
			 xml_Writer.closeTag();
			 xml_Writer.newLine();	*/
			break;
			
		case perRigidLabel: 
			/*
			 xml_Writer.newTag("text");
			 xml_Writer.addAttribute("font-family", "roman");
			 xml_Writer.addAttribute("font-size", "24pt");
			 
			 xml_Writer.openTag();
			 xml_Writer.newLine();
			 xml_Writer.getStream() << "X Tics";
			 xml_Writer.closeTag();
			 xml_Writer.newLine();	
			 */
			break;
			
		case perResidue:
			/*
			 xml_Writer.newTag("text");
			 xml_Writer.addAttribute("font-family", "roman");
			 xml_Writer.addAttribute("font-size", "24pt");
			 
			 xml_Writer.openTag();
			 xml_Writer.newLine();
			 xml_Writer.getStream() << "X Tics";
			 xml_Writer.closeTag();
			 xml_Writer.newLine();	
			 */
			break;			
		};
	
	xml_Writer.closeTag();
	xml_Writer.newLine();	
}

template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::makeYtics(XML_Writer &xml_Writer) {
	xml_Writer.newTag("svg");
	xml_Writer.addAttribute("xmlns", "http://www.w3.org/2000/svg");
		
	xml_Writer.addAttribute("x", "1");
	xml_Writer.addAttribute("y", titleHeight);
		
	xml_Writer.addAttribute("width", yTicsWidth);
	xml_Writer.addAttribute("height", dilutionPlotHeight);	
	xml_Writer.openTag();
	xml_Writer.newLine();
	
	for (I cutoff = minCutoff;
			 cutoff <= maxCutoff;
			 cutoff += (maxCutoff - minCutoff) / 5.0) { // FIXME - prefer to step by 
																									// an even spacing like .1 
																									// instead of .101437
		xml_Writer.newTag("text");
		xml_Writer.addAttribute("font-family", "roman");
		xml_Writer.addAttribute("font-size", "24pt");
		xml_Writer.addAttribute("text-anchor", "end");
		
		xml_Writer.addAttribute("x", "4em"); // FIXME - does this just break in 
																				 // Safari when set to yTicsWidth?
		stringstream yCoordinateStringstream;
		yCoordinateStringstream << getYCoordinateFromCutoff(cutoff);
		xml_Writer.addAttribute("y", yCoordinateStringstream.str());
		
		xml_Writer.openTag();
		xml_Writer.newLine();
		xml_Writer.getStream() << setfill(' ') << std::setw(6) << fixed << 
			setprecision(2) << cutoff;
		xml_Writer.closeTag();
		xml_Writer.newLine();	
		
	}
		
	xml_Writer.closeTag();
	xml_Writer.newLine();	
}

template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::saveStep(XML_Writer &xml_Writer, 
																						I previousCutoff, I cutoff,
																						vector<Piece<I,T>*> &step) {		
	
	for (size_t pieceNumber=0;
			 pieceNumber<step.size();
			 pieceNumber++) {
		
		Piece<I,T>* piece = step[pieceNumber];
		
		set<T> elements = piece->getElementsAndSubElements();
		T previousElementInSegment = *elements.begin();
		T firstElementInSegment = *elements.begin();
		
		// NOTE - remember to use typename with an iterator over a generic collection
		for (typename set<T>::iterator elementIterator = elements.begin(); 
				 elementIterator != elements.end();
				 elementIterator ++) {
			
			T element = *elementIterator;
			
			if (!showElement(element)) {
				continue;
			}
			
			// if the Xindex of the current 
			if (getXindex(element) - getXindex(previousElementInSegment) == 1) {
				std::cout << previousElementInSegment << " -> " << element << std::endl;
				previousElementInSegment = element;

			} else {
				xml_Writer.newTag("rect");
				stringstream xCoordinateStringstream;
				xCoordinateStringstream << getXCoordinateFromElement(firstElementInSegment);
				xml_Writer.addAttribute("x", xCoordinateStringstream.str());
				
				stringstream yCoordinateStringstream;
				yCoordinateStringstream << getYCoordinateFromCutoff(cutoff);
				xml_Writer.addAttribute("y", yCoordinateStringstream.str());
				
				double width = getXCoordinateFromElement(element) -
					getXCoordinateFromElement(firstElementInSegment) + 1.; // FIXME 
					
				std::cout << firstElementInSegment << " - " << previousElementInSegment << " -> " << element << std::endl;

				
				stringstream widthStringstream;
				widthStringstream << width;
				xml_Writer.addAttribute("width", widthStringstream.str());
				
				double height = abs(getYCoordinateFromCutoff(cutoff) -
														getYCoordinateFromCutoff(previousCutoff));
				stringstream heightStringstream;
				heightStringstream << height;
				xml_Writer.addAttribute("height", heightStringstream.str());
				
				xml_Writer.addAttribute("fill", "blue");		
				
				xml_Writer.openTag();		
				xml_Writer.closeTag();
				xml_Writer.newLine();		
				
				firstElementInSegment = element;
				previousElementInSegment = element;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::setupDimensions(Assembly<I,T> *assembly) {
	maxCutoff = assembly->getMaxIndex();
	minCutoff = assembly->getMinIndex();
	
	// rescale maxCutoff and minCutoff a bit to get nicer numbers for the plot
	// FIXME - automagically determine the correct scaling
	maxCutoff = ceil(maxCutoff*5.)/5.;
	minCutoff = floor(minCutoff*50.)/50.;
	
	stringstream heightStringstream;
	heightStringstream << getYCoordinateFromCutoff(minCutoff);
	dilutionPlotHeight = heightStringstream.str();
	
	stringstream widthStringstream;
	widthStringstream << molFramework->total_sites; // TODO - get the width of the dilutionPlot
	dilutionPlotWidth = widthStringstream.str();	
}

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void TIMME_AssemblyOutputSVG<I,T>::visit( Assembly<I,T> *assembly ){
	
	setupDimensions(assembly);
	
  cout << "TIMME_AssemblyOutputSVG currently being implemented." << endl;
  cout << "molFramework pointer has been set. Total sites: " << molFramework->total_sites << endl;
  
	XML_Writer xml_Writer;
	
	string outputFilename =	molFramework->base_name + "_dilution_new.svg";
  ofstream outputFile(outputFilename.c_str());
  xml_Writer.setFile( outputFile.rdbuf() );
	
	xml_Writer.newTag("xml");
	xml_Writer.addAttribute("encoding", "UTF-8");
  xml_Writer.addAttribute("version", "1.0");
  xml_Writer.openEmptyTag("<?");
  xml_Writer.closeEmptyTag("?>");
  xml_Writer.newLine();
	
	xml_Writer.newTag("svg");
  xml_Writer.addAttribute("xmlns", "http://www.w3.org/2000/svg");
  xml_Writer.openTag();
  xml_Writer.newLine();
	
	//	makeTitle(xml_Writer);
	makeXtics(xml_Writer);
	makeYtics(xml_Writer);
	
	xml_Writer.newTag("svg");
	xml_Writer.addAttribute("xmlns", "http://www.w3.org/2000/svg");
	
	xml_Writer.addAttribute("x", yTicsWidth);
	xml_Writer.addAttribute("y", titleHeight);
	
	xml_Writer.addAttribute("width", dilutionPlotWidth);
	xml_Writer.addAttribute("height", dilutionPlotHeight);
	
	xml_Writer.openTag();
	xml_Writer.newLine();
	
	I previousCutoff = 0.0; // FIXME - replace this with first cutoff
		
	AssemblyStepIterator<I,T,Piece<I,T>*> *stepIterator = assembly->getStepIterator();
  while( stepIterator->hasNext() ){
		
    // BMH This is where the meat of the code needs to go. The 'currentStep' variable
    //     is a vector of Piece<I,T>* types. These are the pieces that exist at the 
    //     current or given index (cutoff). To obtain all the 'elements' in a given
    //     piece, call something like: 
    //     set<T> elements = currentStep->getAllElements();
    //     Then do something with the elements. 
		
    vector<Piece<I,T>*> currentStep = stepIterator->next();
		I cutoff = stepIterator->getCurrentStepIndex();
		
		saveStep(xml_Writer, previousCutoff, cutoff, currentStep);
		
		previousCutoff = cutoff;
	}
	
	xml_Writer.closeTag();
  xml_Writer.newLine();
	
	xml_Writer.closeTag();
  xml_Writer.newLine();
	
  xml_Writer.closeTag();
  xml_Writer.newLine();	
	
};

#endif
