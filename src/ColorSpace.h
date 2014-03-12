#ifndef COLORCH_H
#define COLORCH_H

#include "opencv/cv.h"

using namespace cv;

class ColorSpace
{
public:
    int enabled;
    int smooth;
    String colorSpaceType;
	int nChannels;
	String padWith;
};
#endif
