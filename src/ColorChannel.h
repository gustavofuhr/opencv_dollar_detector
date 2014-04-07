#ifndef COLORCH_H
#define COLORCH_H

#include "cv.h"

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX};

using namespace cv;

class ColorChannel
{
public:
    int enabled;
    int smooth;
    String colorSpaceType;
	int nChannels;
	String padWith;

	Mat convConst(Mat, int, int);
	Mat convolutionWithTriangleFilter(Mat);
	Mat rgbConvert(Mat);
	Mat rgb2luv(Mat);
	Mat rgb2hsv(Mat);
	Mat rgb2gray(Mat);
};
#endif
