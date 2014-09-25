#ifndef INFO_H
#define INFO_H

#include "ColorChannel.h"
#include "GradientMagnitudeChannel.h"
#include "QuantizedGradientChannel.h"

//return class for the chnsCompute function
//it contains the channels and the result of the computations
//i might think of another way of doing this (or another name, at least), but i need to know how it is used first.
class Info
{
public:
	cv::Mat image;
	cv::Mat gradientMagnitude;
	std::vector<cv::Mat> gradientHistogram;
};
#endif
