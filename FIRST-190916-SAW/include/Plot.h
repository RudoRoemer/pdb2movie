#ifndef PLOT_H_
#define PLOT_H_

#include <string>
#include <map>
#include <vector>

#include "Color.h"

using namespace std;

struct Point{ // TODO - replace with a Point class
	float x;
	float y;
};

struct Rectangle{ // TODO - replace with a Rectangle class
	Point point;
	float width;
	float height;
};

typedef std::map<string, string> Style; // TODO - replace with a better abstraction 
typedef std::map<string, string> Attributes; // TODO - replace with a better abstraction
typedef std::vector<Point> Points;

class Plot {
public:
	virtual ~Plot();
	
	virtual void beginDocument(int width, int height);
	virtual void endDocument();

	virtual void drawRectangle(Style style, Rectangle rectangle);
  virtual void drawRectangle(Style style, Attributes attributes, Rectangle rectangle);
  
  virtual void drawPolyline(Style style, Points points);
  virtual void drawPolyline(Style style, Attributes attributes, Points points);
  
  virtual void drawPolygon(Style style, Points points);
  virtual void drawPolygon(Style style, Attributes attributes, Points points);
  
	virtual void drawText(Point point, Style style, string text);
	virtual void drawText(Point point, Style style, Attributes attributes, string text);
	
  virtual string getColor(unsigned int colorNumber);
  
	virtual void script(string script);
};

#endif /*PLOT_H_*/
