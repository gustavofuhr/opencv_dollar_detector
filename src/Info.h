#ifndef INFO_H
#define INFO_H

#include "ColorChannels.h"
#include "GradientMagnitude.h"
#include "GradientHistogram.h"

class Info
{
public:
	ColorChannels colorCh;
	GradientMagnitude gradMag;
	GradientHistogram gradHist;
};
#endif