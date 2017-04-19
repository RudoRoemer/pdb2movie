#include "XML_Writer.h"

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::setFile( streambuf *strmBuf ){
    file.rdbuf( strmBuf );
}

////////////////////////////////////////////////////////////////////////////////
ostream& XML_Writer::getStream(){ 
  return file; 
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::newTag( string name ){

  xml_Tags.push_back( new XML_Tag( name ) );
};

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::addAttribute( string name, string value ){

  if( !xml_Tags.empty() ){
    XML_Tag *currentTag = xml_Tags.back();
    currentTag->addAttribute( name, value );
  }
  else{
    cout << "Error: Adding attribute to non-existant tag" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::openTag( string openChars , string closeChars ){
  
  if( xml_Tags.empty() )
    return;
  
  XML_Tag *currentTag = xml_Tags.back();
  
  file << openChars
       << currentTag->getName();
  
  // Print any attributes
  file << currentTag->getAttributes();
  
  // Close the tag      
  file << closeChars;
  
  ++tagSanityCheck;
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::closeTag( string openChars, string closeChars ){
  
  if( xml_Tags.empty() )
    return;
  
  --tagSanityCheck;
  
  XML_Tag *currentTag = xml_Tags.back();
  file << openChars
       << currentTag->getName()
       << closeChars;
  
  xml_Tags.pop_back();
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::openEmptyTag( string openChars ){
  
  if( xml_Tags.empty() )
    return;
  
  XML_Tag *currentTag = xml_Tags.back();
  
  // Open the element
  file << openChars
       << currentTag->getName();
  
  // Print any attributes
  file << currentTag->getAttributes();
  
  ++tagSanityCheck;
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::closeEmptyTag( string closeChars ){
  
  // Close the tag      
  file << closeChars;
  
  xml_Tags.pop_back();
  
  --tagSanityCheck;
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::openComment( string openChars ){

  file << openChars << endl;
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::closeComment( string closeChars ){

  file << closeChars << endl;
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::indent( unsigned int indent ){ 

  file << setw(indent) << ""; 
}

////////////////////////////////////////////////////////////////////////////////
void XML_Writer::newLine(){ 

  file << endl; 
}

