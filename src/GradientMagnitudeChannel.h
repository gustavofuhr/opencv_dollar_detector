#ifndef GRADMAG_H
#define GRADMAG_H

#include "utils.h"

class GradientMagnitudeChannel
{
public:
	int enabled; //if true enable gradient magnitude channel
	int colorChannelIndex; //if>0 color channel to use for grad computation
	int normalizationRadius; //normalization radius for gradient
	float normalizationConstant; //normalization constant for gradient
	int full; //if true compute angles in [0,2*pi) else in [0,pi)
	int nChannels;
	cv::String padWith;

	void readGradientMagnitude(cv::FileNode);
	std::vector<float*> mGradMag(float* source, int rows, int cols, int channel);
	void gradMagNorm(float *M, float *S, int h, int w); 
};

#endif
