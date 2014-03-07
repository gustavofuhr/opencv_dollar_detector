#ifndef COLORCH_H
#define COLORCH_H

#include "opencv/cv.h"

using namespace cv;

class ColorChannels
{
public:
    int enabled;
    int smooth;
	String colorSpace;
	int nChannels;
	String padWith;
};
#endif
