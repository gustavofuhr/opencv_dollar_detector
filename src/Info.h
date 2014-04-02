#ifndef INFO_H
#define INFO_H

#include "ColorChannel.h"
#include "GradientMagnitudeChannel.h"
#include "QuantizedGradientChannel.h"


//this whole class might be useless, it is the same as ChannelFeatures
class Info
{
public:
    	ColorChannel colorCh;
	GradientMagnitudeChannel gradMag;
	GradientHistogramChannel gradHist;
};
#endif
