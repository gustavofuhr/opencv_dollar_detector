#include "QuantizedGradientChannel.h"
#include "GradientMagnitudeChannel.h"
#include "ColorChannel.h"

class ChannelFeatures
{
public:
	int shrink;
   	ColorChannel pColor;
	GradientMagnitudeChannel pGradMag;
	QuantizedGradientChannel pGradHist;
	//pCustom //custom channels (i dont think it is used)
	int complete;
};
