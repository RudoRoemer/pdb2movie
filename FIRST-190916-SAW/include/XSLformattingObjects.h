#ifndef XSLFORMATTINGOBJECTS_H_
#define XSLFORMATTINGOBJECTS_H_

#include <fstream>
#include <iostream>

class XSLformattingObjects {
public:
	
	XSLformattingObjects();
	
	void beginDocument(std::ofstream &document);
	void endDocument(std::ofstream &document);
	
	void setupSimplePageLayout(std::ofstream &document,
							   float widthInCM,
							   float heightInCM);
	
	void beginSimplePageSequence(std::ofstream &document);
	void endSimplePageSequence(std::ofstream &document);

	void beginFlow(std::ofstream &document);
	void endFlow(std::ofstream &document);
	
	void beginBlock(std::ofstream &document);
	void endBlock(std::ofstream &document);
};

#endif /*XSLFORMATTINGOBJECTS_H_*/
