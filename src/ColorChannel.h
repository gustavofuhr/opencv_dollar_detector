#ifndef COLORCH_H
#define COLORCH_H

#include "cv.h"
#include "sse.hpp"
//#include "wrappers.hpp"

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

	void readColorChannel(FileNode);
	Mat convolution(Mat, int, int, int);
	void triangleFilterConvolution(float*,float*,int,int,int,int, int);
	void convTri1Y(float*, float*, int, float, int);
	void convTri1(float*,float*,int,int,int,float,int);
	Mat rgbConvert(Mat);
	Mat rgb2luv(Mat);
	Mat rgb2hsv(Mat);
	Mat rgb2gray(Mat);
};
#endif
