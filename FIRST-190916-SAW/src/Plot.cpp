#include "Plot.h"

Plot::~Plot() {
}

void Plot::beginDocument(int width, int height) {
}

void Plot::endDocument() {
}

void Plot::drawRectangle(Style style, Rectangle rectangle) {
}

void Plot::drawRectangle(Style style, Attributes attributes, Rectangle rectangle) {
}

void Plot::drawPolyline(Style style, Points points) {
  
}

void Plot::drawPolyline(Style style, Attributes attributes, Points points) {
  
}

void Plot::drawPolygon(Style style, Points points) {
  
}

void Plot::drawPolygon(Style style, Attributes attributes, Points points) {
  
}

void Plot::drawText(Point point, Style style, string text) {
}

void Plot::drawText(Point point, Style style, Attributes attributes, string text) {
}

void Plot::script(string script) {
} 

string Plot::getColor(unsigned int colorNumber) {
  return "";
}

