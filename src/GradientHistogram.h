#ifndef GRADHIST_H
#define GRADHIST_H

#include "opencv/cv.h"

using namespace cv;

class GradientHistogram
{
public:
    int enabled;
	int* binSize;
	int nOrients;
	int softBin;
	int useHog;
	double clipHog;
	int nChannels;
	String padWith;
};

#endif
