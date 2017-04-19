// Brandon Hespenheide (c) 2007
////////////////////////////////////////////////////////////////////////////////

#ifndef _OUTPUT_GRAPHML_
#define _OUTPUT_GRAPHML_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>

#include "Assembly.h"
#include "XML_Writer.h"
#include "AssemblyPieceIterator.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
class OutputGraphML : public AssemblyVisitor<I,T> {
	
public:
	MolFramework *molFramework;
	
	OutputGraphML(MolFramework *_molFramework) :
    molFramework( _molFramework )
	{	};
  ~OutputGraphML()
  {};
	
public:
		virtual void visit( Assembly<I,T> *assembly );
  
};

////////////////////////////////////////////////////////////////////////////////
template <class I, class T>
void OutputGraphML<I,T>::visit( Assembly<I,T> *assembly ){
	
  AssemblyPieceIterator<I,T,Piece<I,T>*> *pieceIterator = assembly->getPieceIterator();
	
  XML_Writer xml_Writer;
	
  ofstream outputFile("timmeCutoffs.xml");
  xml_Writer.setFile( outputFile.rdbuf() );
  
  // Print XML Header stuff
  //////////////////////////////////////////////////////////////////////
  xml_Writer.newTag("xml");
  xml_Writer.addAttribute("version", "1.0");
  xml_Writer.addAttribute("encoding", "UTF-8");
  xml_Writer.openEmptyTag("<?");
  xml_Writer.closeEmptyTag("?>");
  xml_Writer.newLine();
	
  xml_Writer.newTag("graphml");
  xml_Writer.addAttribute("xmlns", "http://graphml.graphdrawing.org/xmlns");
  xml_Writer.addAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
  xml_Writer.addAttribute("xsi:schemaLocation", "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd");
  xml_Writer.openTag();
  xml_Writer.newLine();
	
  xml_Writer.openComment();
  xml_Writer.getStream() << "Due to the nature of the TIMME dilution, there will be a cluster(s) that" << endl;
  xml_Writer.getStream() << "is not a subset of any other cluster. However, the GraphML specification" << endl;
  xml_Writer.getStream() << "requires that every edge have a source and a target. Therefore, these clusters" << endl;
  xml_Writer.getStream() << "in the TIMME dilution (there will always be >= 1 of these) have the source" << endl;
  xml_Writer.getStream() << "and target values equal to eachother." << endl;
  xml_Writer.closeComment();
	
  xml_Writer.newTag("graph");
  xml_Writer.addAttribute("id", "test");
  xml_Writer.addAttribute("edgedefault", "undirected");
  xml_Writer.openTag();
  xml_Writer.newLine();
	
  while( pieceIterator->hasNext() ){
		
    Piece<Cutoff, SiteID> *currentPiece = pieceIterator->next();
    
    xml_Writer.newTag("node");
    xml_Writer.addAttribute("id", convertToString( currentPiece->getID()) );
    
    if( currentPiece->totalElements() ){
			
      xml_Writer.openTag();
      xml_Writer.newLine();
			
      set<SiteID> pieceElements = currentPiece->getElements();
      
			xml_Writer.newTag("data");
			xml_Writer.openTag();
			xml_Writer.newLine();				

      for(set<SiteID>::iterator pieceElementsIter = pieceElements.begin();
					pieceElementsIter != pieceElements.end();
					pieceElementsIter++){
				
				xml_Writer.newTag("element");
			
				xml_Writer.openTag();
				xml_Writer.getStream() << *pieceElementsIter;				
				xml_Writer.closeTag();	
				xml_Writer.newLine();
      }
			
			xml_Writer.closeTag();	
			xml_Writer.newLine();			
			
      xml_Writer.closeTag();
      xml_Writer.newLine();
    }
    else{
      xml_Writer.openEmptyTag();    
      xml_Writer.closeEmptyTag();
      xml_Writer.newLine();
    }
    
    Piece<Cutoff, SiteID> *superPiece = currentPiece->getSuperPiece();
		
		string source;
		string target = convertToString( currentPiece->getID());
		
    if( superPiece != NULL ){
			source = convertToString( superPiece->getID());
			
		} else {
			source = target;
		}
		
		stringstream cutoffStringstream;
		cutoffStringstream << currentPiece->getIndex();
		string cutoff = cutoffStringstream.str();
		
		xml_Writer.newTag("edge");
		xml_Writer.addAttribute("source", source );
		xml_Writer.addAttribute("target", target );
		xml_Writer.addAttribute("cutoff", cutoff );
		xml_Writer.openTag(); 
		xml_Writer.newLine();
		
/*		xml_Writer.newTag("data");
		xml_Writer.openTag();
		xml_Writer.newLine();
		
		xml_Writer.newTag("index");
		xml_Writer.openTag();
		xml_Writer.getStream() << currentPiece->getIndex();
		xml_Writer.closeTag();
		xml_Writer.newLine();
		
		xml_Writer.closeTag();	
		xml_Writer.newLine();*/
		
		xml_Writer.closeTag();
		xml_Writer.newLine();
		
  }
	
  xml_Writer.closeTag();
  xml_Writer.newLine();
  xml_Writer.closeTag();
  xml_Writer.newLine();
};

#endif
