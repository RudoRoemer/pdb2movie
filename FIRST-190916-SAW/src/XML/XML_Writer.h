#ifndef _XML_WRITER_
#define _XML_WRITER_

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <assert.h>
#include "XML_Tag.h"

////////////////////////////////////////////////////////////////////////////////
// Description:
// FILO implementation of an XML writer. 
////////////////////////////////////////////////////////////////////////////////
class XML_Writer {

private:
  vector< XML_Tag* > xml_Tags;
  ostream file;
  int tagSanityCheck;

public:
  XML_Writer() :
    file( cout.rdbuf() ), 
    tagSanityCheck(0)
  {};
  ~XML_Writer()
  {
    assert( !tagSanityCheck );
  };

public:
  
  void setFile( streambuf *strmBuf );
  ostream& getStream();
  void newTag( string name );
  void addAttribute( string name, string value );
  void openTag( string openChars = "<", string closeChars = ">" );
  void closeTag( string openChars = "</", string closeChars = ">" );
  void openEmptyTag( string openChars = "<" );
  void closeEmptyTag( string closeChars = "/>" );
  void openComment( string openChars = "<!--" );
  void closeComment( string closeChars = "-->" );
  void indent( unsigned int indent );
  void newLine();
};

template <class T>
string convertToString( T t ){
  ostringstream o;
  if( !( o << t ) ){
    exit(12);
  }
  return o.str();
};

#endif
