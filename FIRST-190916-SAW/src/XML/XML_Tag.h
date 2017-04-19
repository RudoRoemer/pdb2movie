#ifndef _XML_TAG_
#define _XML_TAG_

#include <map>
#include <string>
#include <sstream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
class XML_Tag {

private:
  string name;
  map<string, string> attributes;

public:
  XML_Tag( string name_, unsigned short indent = 0, bool lineBreak_ = true ) :
    name( name_ )
  {};
  ~XML_Tag()
  {};
  
public:
  void addAttribute( string name, string value );
  string getAttributes();
  string getName();
};

#endif
