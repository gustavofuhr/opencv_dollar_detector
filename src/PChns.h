#include "GradientHistogram.h"
#include "GradientMagnitude.h"
#include "ColorChannels.h"

class PChns
{
public:
	int shrink;
	ColorChannels pColor;
	GradientMagnitude pGradMag;
	GradientHistogram pGradHist;
	//pCustom
	int complete;
};
