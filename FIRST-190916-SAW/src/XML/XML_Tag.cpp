#include "XML_Tag.h"

////////////////////////////////////////////////////////////////////////////////
void XML_Tag::addAttribute( string name, string value ){

  attributes.insert( pair<string,string> (name, value) );  
}

////////////////////////////////////////////////////////////////////////////////
string XML_Tag::getAttributes(){
  
  ostringstream attributeList;

  if( attributes.size() > 0 ){
    map<string,string>::iterator attributesIter = attributes.begin();
    while( attributesIter != attributes.end() ){
      attributeList << " " 
		    << attributesIter->first
		    << "=\"" 
		    << attributesIter->second
		    << "\""; 
      attributesIter++;
    }
  }
  
  return attributeList.str();
}

////////////////////////////////////////////////////////////////////////////////
string XML_Tag::getName(){

  return name;
}
  

