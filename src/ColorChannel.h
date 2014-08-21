#ifndef COLORCH_H
#define COLORCH_H

#include "utils.h"

class ColorChannel
{
public:
	int enabled;
	int smoothingRadius;
	cv::String colorSpaceType;
	int nChannels;
	cv::String padWith;

	void readColorChannel(cv::FileNode);
	cv::Mat rgbConvert(cv::Mat);
};
#endif
