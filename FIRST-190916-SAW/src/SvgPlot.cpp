#include <sstream>
#include "global_defs.h"
#include "SvgPlot.h"

using namespace std;

SvgPlot::SvgPlot(string filename) {
	this->document = new std::ofstream(filename.c_str());
	
}

SvgPlot::~SvgPlot() {
	if (this->document != NULL) {
		this->document->close();
	}
}

void SvgPlot::drawRectangle(Style style, Rectangle rectangle) {
	Attributes attributes;
  
  drawRectangle(style, attributes, rectangle);
}
  
void SvgPlot::drawRectangle(Style style, Attributes attributes, Rectangle rectangle) {
	Point point = rectangle.point;
	
	stringstream x;
	x << point.x << "";
	attributes["x"] = x.str();
	
	stringstream y;
	y << point.y << "";
	attributes["y"] = y.str();
	
	stringstream width;
	width << rectangle.width << "";
	attributes["width"] = width.str();
	
	stringstream height;
	height << rectangle.height  << "";
	attributes["height"] = height.str();

	stringstream styleStringStream;
	for (Style::iterator styleIterator = style.begin();
		 styleIterator != style.end();
		 styleIterator++) {
		styleStringStream << styleIterator->first << ": " << styleIterator->second << "; ";
	}
	
	attributes["style"] = styleStringStream.str();
	
	openTag("rect", attributes);
	closeTag("rect");
}

void SvgPlot::drawText(Point point, Style style, string text) {
	Attributes attributes;

	drawText(point, style, attributes, text);
}

void SvgPlot::drawText(Point point, Style style, Attributes attributes, string text) {
	stringstream x;
	x << point.x << "";
	attributes["x"] = x.str();
	
	stringstream y;
	y << point.y << "";
	attributes["y"] = y.str();
	
	stringstream styleStringStream;
	for (Style::iterator styleIterator = style.begin();
		 styleIterator != style.end();
		 styleIterator++) {
		styleStringStream << styleIterator->first << ": " << styleIterator->second << "; ";
	}
	
	attributes["style"] = styleStringStream.str();
	
	openTag("text",  attributes);
	
	content(text);
	
	closeTag("text");
}


// TODO - refactor these into an XML superclass?
void SvgPlot::openTag(string tagName, std::map<string, string> mapFromAttributeNameToValue) {
	beginOpenTag(tagName);
	
	for (std::map<string, string>::iterator attributeNameToValueIterator = mapFromAttributeNameToValue.begin();
		 attributeNameToValueIterator != mapFromAttributeNameToValue.end();
		 attributeNameToValueIterator++) {
		string attributeName = (*attributeNameToValueIterator).first;
		string attributeValue = (*attributeNameToValueIterator).second;
		
		addAttribute(attributeName, 
					 attributeValue);
	}
	
	endOpenTag();
}

void SvgPlot::beginOpenTag(string tagName) {
	*document << "<" << tagName << " ";
}

void SvgPlot::addAttribute(string attributeName, string attributeValue) {
	*document << " " << attributeName  << "=\"" << attributeValue << "\" ";
}

void SvgPlot::endOpenTag() {
	*document << ">" << std::endl;;
}

void SvgPlot::content(string content) {
	*document << content;
}


void SvgPlot::closeTag(string name) {
	*document << "</" << name << ">" << std::endl;
}

void SvgPlot::beginScript() {
	std::map<string, string> attributeNameToValue;
	attributeNameToValue["type"] = "text/ecmascript";
	
	openTag("script", attributeNameToValue);
	*document << "<![CDATA[" << std::endl;
}

void SvgPlot::endScript() {
	*document << "]]>" << std::endl;
	closeTag("script");
}

string SvgPlot::getColor(ColorNumber colorNumber) {
  
  return Color::getSVGColor(colorNumber);
}

void SvgPlot::script(string scriptContent) {
	beginScript();
	
	content(scriptContent);
	
	endScript();
}

void SvgPlot::beginDocument(int width, int height) {
  *document << "<?xml version=\"1.0\" standalone=\"no\"?>" << std::endl;
	*document << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"" << std::endl;
	*document << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << std::endl;

  
/*  *document << " <?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
  *document << " <!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG Tiny 1.1//EN\"" << std::endl;
  *document << " \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11-tiny.dtd\">" << std::endl;*/
    
	std::map<string, string> svgAttributes;
	svgAttributes["x"] = "0";
	svgAttributes["y"] = "0";

  stringstream widthStringstream;
  widthStringstream << width  << "";
	svgAttributes["width"] = widthStringstream.str();
  
  stringstream heightStringstream;
  heightStringstream << height << "";
	svgAttributes["height"] = heightStringstream.str();
	
	svgAttributes["version"] = "1.1";
	svgAttributes["xmlns"] = "http://www.w3.org/2000/svg";
	svgAttributes["xmlns:xlink"] = "http://www.w3.org/1999/xlink";
	
	openTag("svg", svgAttributes);
}

void SvgPlot::endDocument() {
	closeTag("svg");
}

void SvgPlot::drawPolyline(Style style, Points points) {	  
  Attributes attributes;
  
  drawPolyline(style, attributes, points);
}

void SvgPlot::drawPolyline(Style style, Attributes attributes, Points points) {
  
	stringstream styleStringStream;
	for (Style::iterator styleIterator = style.begin();
       styleIterator != style.end();
       styleIterator++) {
		styleStringStream << styleIterator->first << ": " << styleIterator->second << "; ";
	}
	
	attributes["style"] = styleStringStream.str();
	
  stringstream pointsStringStream;
  for (Points::iterator point = points.begin();
       point != points.end();
       point ++) {
    
    pointsStringStream << point->x << "," << point->y << " ";
  }
  
  attributes["points"] = pointsStringStream.str();
  
	openTag("polyline", attributes);
	closeTag("polyline");
  
}

void SvgPlot::drawPolygon(Style style, Points points) {	  
  Attributes attributes;

  drawPolygon(style, attributes, points);
}

void SvgPlot::drawPolygon(Style style, Attributes attributes, Points points) {

	stringstream styleStringStream;
	for (Style::iterator styleIterator = style.begin();
       styleIterator != style.end();
       styleIterator++) {
		styleStringStream << styleIterator->first << ": " << styleIterator->second << "; ";
	}
	
	attributes["style"] = styleStringStream.str();
	
  stringstream pointsStringStream;
  for (Points::iterator point = points.begin();
       point != points.end();
       point ++) {
    
    pointsStringStream << point->x << "," << point->y << " ";
  }
  
  attributes["points"] = pointsStringStream.str();
    
	openTag("polygon", attributes);
	closeTag("polygon");
  
}

