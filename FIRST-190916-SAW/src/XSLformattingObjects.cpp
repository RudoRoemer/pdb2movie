#include "XSLformattingObjects.h"
 
// FIXME - need to find a better way to do the XML here (this is ugly and error-prone)

XSLformattingObjects::XSLformattingObjects() {
	
}

void XSLformattingObjects::beginDocument(std::ofstream & document) {
	document << "<?xml version='1.0' encoding='utf-8'?>" << std::endl;
	document << "<fo:root xmlns:fo='http://www.w3.org/1999/XSL/Format'>" << std::endl;
}

void XSLformattingObjects::endDocument(std::ofstream & document) {
	document << "</fo:root>" << std::endl;
}

void XSLformattingObjects::setupSimplePageLayout(std::ofstream & document,
												 float widthInCM,
												 float heightInCM) {
	document << "<fo:layout-master-set>" << std::endl;
	
	// FIXME - hard coded constants (make them parameters)
	document << "<fo:simple-page-master master-name='simple'" << std::endl;
	document << " page-height='"<< heightInCM << "cm'" << std::endl;
	document << " page-width='"<< widthInCM << "cm'" << std::endl;
	document << " margin-top='0cm'" << std::endl;
	document << " margin-bottom='0cm'" << std::endl;
	document << " margin-left='0cm'" << std::endl;
	document << " margin-right='0cm'" << std::endl;
	document << ">" << std::endl;
	
	document << "<fo:region-body margin-top='0cm'/>" << std::endl;
	document << "<fo:region-before extent='0cm'/>" << std::endl;
	document << "<fo:region-after extent='0cm'/>" << std::endl;

	document << "</fo:simple-page-master>" << std::endl;
	document << "</fo:layout-master-set>" << std::endl;
}

void XSLformattingObjects::beginSimplePageSequence(std::ofstream & document) {
	document << "<fo:page-sequence master-reference='simple'>" << std::endl;
}

void XSLformattingObjects::endSimplePageSequence(std::ofstream & document) {
	document << "</fo:page-sequence>" << std::endl;
}

void XSLformattingObjects::beginFlow(std::ofstream & document) {
	document << "<fo:flow flow-name='xsl-region-body'>" << std::endl;
}

void XSLformattingObjects::endFlow(std::ofstream & document) {
	document << "</fo:flow>" << std::endl;
}

void XSLformattingObjects::beginBlock(std::ofstream & document) {
	// FIXME - hard coded constants (make them parameters)
	document << "<fo:block";
	document << " font-size='18pt'" << std::endl;
	document << " font-family='sans-serif'" << std::endl;
	document << " line-height='24pt'" << std::endl;
	document << " space-after.optimum='15pt'" << std::endl;
	document << " background-color='white'" << std::endl;
	document << " color='black'" << std::endl;
	document << " text-align='center'" << std::endl;
	document << " padding-top='3pt'" << std::endl;
	document << ">" << std::endl;
}

void XSLformattingObjects::endBlock(std::ofstream & document) {
	document << "</fo:block>" << std::endl;
}


/*
 dilutionStripePlotFile << "<>" << std::endl;
 dilutionStripePlotFile << "</>" << std::endl;
 */



