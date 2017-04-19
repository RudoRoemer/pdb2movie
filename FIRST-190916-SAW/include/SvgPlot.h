#ifndef SVG_PLOT_H_
#define SVG_PLOT_H_

#include <fstream>
#include <map>

#include "Plot.h"

using namespace std;

class SvgPlot : public Plot {
public:
	SvgPlot(string filename);
	~SvgPlot();
	
	void beginDocument(int width, int height);
	void endDocument();

	void drawRectangle(Style style, Rectangle rectangle);
	void drawRectangle(Style style, Attributes attributes, Rectangle rectangle);
  
  void drawPolyline(Style style, Points points);
  void drawPolyline(Style style, Attributes attributes, Points points);
  
  void drawPolygon(Style style, Points points);
  void drawPolygon(Style style, Attributes attributes, Points points);

	void drawText(Point point, Style style, string text);	
	void drawText(Point point, Style style, Attributes attributes, string text);

	string getColor(unsigned int colorNumber);
  
	void script(string scriptContent);

public:
	void openTag(string tagName, 
               std::map<string, string> mapFromAttributeNameToValue);
	
	void beginScript();
	void endScript();
		
	void content(string content);
	
	void closeTag(string name);
  
private:
  void beginOpenTag(string tagName);
	void addAttribute(string name, string value);
	void endOpenTag();
	
	std::ofstream * document;
};

#endif /*SVG_PLOT_H_*/
