#ifndef CHFEAT_H
#define CHFEAT_H

#include "QuantizedGradientChannel.h"
#include "GradientMagnitudeChannel.h"
#include "ColorChannel.h"

class ChannelFeatures
{
public:
	int shrink; //amount to subsample the computed channels
  	ColorChannel pColor;
	GradientMagnitudeChannel pGradMag;
	QuantizedGradientChannel pGradHist;
	//pCustom //to use custom channels (not implemented)
	int complete;

	void readChannelFeatures(cv::FileNode);
};

#endif
