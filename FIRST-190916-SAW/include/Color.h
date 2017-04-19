#ifndef _COLOR_
#define _COLOR_

#include <string>
#include <map>

#include "SiteID.h"

using namespace std;

struct RGBfractionComponents{float red; float green; float blue;};
struct RGBintegerComponents{int red; int green; int blue;};

struct HSBfractionComponents {
  float hue; 
  float saturation; 
  float brightness;
};

typedef unsigned int ColorNumber;

class Color {
public:		
	static string getHexString(unsigned int n);

	static string getPymolColor(unsigned int n);
	static string getJmolColor(unsigned int n);
	static string getSVGColor(unsigned int n);
	static string getPSColor(unsigned int n);
	
	static RGBfractionComponents getRGBfractionComponents(ColorNumber colorNumber);
	static RGBintegerComponents getRGBintegerComponents(ColorNumber colorNumber);
	
	static HSBfractionComponents getHSBfractionComponents(ColorNumber colorNumber);
	
	static unsigned int convertFractionToInt(float fractionValue, unsigned int maxValue);
	
	static RGBfractionComponents convertFromHSBtoRGB(HSBfractionComponents hsbFractionComponents);
	
	static float HUEtoRGB(float m1, float m2, float hue);

        static void setUseSpectrum(bool);
        static void setNumberOfBins(unsigned int);

        static std::map<ColorNumber, SiteID> mapFromColorNumberToMeanSiteID;
        static std::map<ColorNumber, SiteID> mapFromColorNumberOtRMSDfromMeanSiteID;

private:
        static bool useSpectrum;
        static unsigned int numberOfBins;
};

#endif 
