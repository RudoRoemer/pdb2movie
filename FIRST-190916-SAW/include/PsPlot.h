#ifndef PS_PLOT_H_
#define PS_PLOT_H_

#include "Plot.h"

class PsPlot : public Plot {
  
 public:
  PsPlot(string filename);
  ~PsPlot();

  void beginDocument( int width = 0, int height = 0);
  void endDocument();
  
  void drawRectangle(Style style, Rectangle rectangle);
  void drawRectangle(Style style, Attributes attributes, Rectangle rectangle);
  
  void drawPolygon(Style style, Points points);
  void drawPolygon(Style style, Attributes attributes, Points points);

  void drawPolyline(Style style, Points points);
  void drawPolyline(Style style, Attributes attributes, Points points);
  
  void drawText(Point point, Style style, string text);	
  void drawText(Point point, Style style, Attributes attributes, string text);
  
  string getColor(unsigned int colorNumber);
  
 private:
  ofstream *document;
  int width;
  int height;
};

#endif /*PS_PLOT_H_*/
