#include "GradientHistogram.h"
#include "GradientMagnitude.h"
#include "ColorSpace.h"

class Channel
{
public:
	int shrink;
    ColorSpace pColor;
	GradientMagnitude pGradMag;
	GradientHistogram pGradHist;
	//pCustom
	int complete;
};
