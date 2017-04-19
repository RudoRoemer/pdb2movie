#include "global_defs.h"
#include "Parameters.h"
#include "PsPlot.h"

extern const Parameters parameters;

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
PsPlot::PsPlot( string filename ) {
  document = new ofstream( filename.c_str() );
  *document << showpoint << setiosflags(ios::fixed);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
PsPlot::~PsPlot(){
  if( document != NULL ) {
    document->close();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void PsPlot::beginDocument( int width, int height ){

  this->width = width;
  this->height = height;
  *document << "%!PS-Adobe-2.0 EPSF-2.0" << endl;
  *document << "%%Title: FIRST dilution plot" << endl;
  *document << "%%Creator: FIRST" << endl;
  *document << "%%CreationDate: 2006" << endl;
  //*document << "%%BoundingBox:0 0 620 800" << endl;
  *document << "%%BoundingBox:0 0 " << width << " " << height << endl;
  *document << "%%DocumentFonts: Times-Roman" << endl;
  *document << "%%EndComments" << endl;
  *document << "%%BeginProlog" << endl;
  *document << "/Col { sethsbcolor } bind def" << endl;
  *document << "%% These are the colors for the flexibility scale and the" << endl;
  *document << "%% geometric primitives that are used to create the plot." << endl;

  *document << "/BLACK {0.0000  0.0000  0.0000 setrgbcolor } def" << endl; 
  *document << "/BLUE  {0.0000  0.0000  1.0000 setrgbcolor } def" << endl;
  *document << "/RED   {1.0000  0.0000  0.0000 setrgbcolor } def" << endl;

  *document << "/poly { moveto lineto lineto lineto fill } bind def" << endl;
  *document << "/rect { 8 copy poly moveto moveto moveto moveto closepath stroke } bind def" << endl;
  *document << "/line { moveto lineto stroke } bind def" << endl; 
  *document << "%%EndProlog" << endl;
  *document << "%%Page:    1   1" << endl;
} 

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void PsPlot::endDocument(){

  *document << "showpage" << endl;
  *document << "%%EOF" << endl;

}

////////////////////////////////////////////////////////////////////////////////
// Description:
////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawRectangle( Style style, Rectangle rectangle ){
  
  string spacer = " ";
  
  *document << setprecision(2);

  // If the height of the rectangle is 0.0, then we're drawing a line. 
  //////////////////////////////////////////////////////////////////////
  if( rectangle.height == 0.5 ){
    *document << "BLACK " << endl;
    *document << rectangle.point.x << spacer
	      << height -rectangle.point.y << spacer
	      << rectangle.point.x + rectangle.width << spacer
	      << height -rectangle.point.y << spacer
	      << "line" << endl;
  }
  else{
    //TODO get color from attributes
    *document << style["fill"];
    *document << rectangle.point.x << spacer
	      << height -rectangle.point.y << spacer
	      << rectangle.point.x + rectangle.width << spacer
	      << height -rectangle.point.y << spacer
	      << rectangle.point.x + rectangle.width << spacer
	      << height -rectangle.point.y - rectangle.height << spacer
	      << rectangle.point.x << spacer 
	      << height -rectangle.point.y - rectangle.height << spacer
	      << "rect" << endl;
  }
  
}

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawRectangle(Style style, Attributes attributes, Rectangle rectangle){};

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawPolyline(Style style, Points points) {};

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawPolyline(Style style, Attributes attributes, Points points) {};

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawPolygon(Style style, Points points) {};

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawPolygon(Style style, Attributes attributes, Points points) {};

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawText(Point point, Style style, string text){};	

////////////////////////////////////////////////////////////////////////////////
void PsPlot::drawText(Point point, Style style, Attributes attributes, string text){};

////////////////////////////////////////////////////////////////////////////////
string PsPlot::getColor( unsigned int colorNumber ){
  
  return( Color::getPSColor( colorNumber ) );
}

