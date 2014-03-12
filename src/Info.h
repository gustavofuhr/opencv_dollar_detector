#ifndef INFO_H
#define INFO_H

#include "ColorSpace.h"
#include "GradientMagnitude.h"
#include "GradientHistogram.h"

class Info
{
public:
    ColorSpace colorCh;
	GradientMagnitude gradMag;
	GradientHistogram gradHist;
};
#endif
