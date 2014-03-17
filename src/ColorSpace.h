#ifndef COLORCH_H
#define COLORCH_H

#include "opencv\cv.h"

using namespace cv;

class ColorSpace
{
public:
    int enabled;
    int smooth;
    String colorSpaceType;
	int nChannels;
	String padWith;

    Mat rgbConvert(Mat);
    Mat rgb2luv(Mat);
    Mat rgb2hsv(Mat);
    Mat rgb2gray(Mat);
};
#endif
