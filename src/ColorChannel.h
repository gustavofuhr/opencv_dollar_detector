#ifndef COLORCH_H
#define COLORCH_H

#include "opencv.hpp"
#include "sse.hpp"
#include <stdio.h>
#include "utils.h"

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX};

class ColorChannel
{
public:
	int enabled;
	int smooth;
	cv::String colorSpaceType;
	int nChannels;
	cv::String padWith;

	void readColorChannel(cv::FileNode);
	cv::Mat convolution(cv::Mat, int, int, int);
	void triangleFilterConvolution(float*, float*, int, int, int, int, int);
	void convTriY(float*, float*, int, int, int);
	void convTri1Y(float*, float*, int, float, int);
	void convTri1(float*, float*, int, int, int, float, int);
	cv::Mat rgbConvert(cv::Mat);
	cv::Mat rgb2luv(cv::Mat);
	cv::Mat rgb2hsv(cv::Mat);
	cv::Mat rgb2gray(cv::Mat);
};
#endif
