#ifndef GRADHIST_H
#define GRADHIST_H

#include "utils.h"

#define PI 3.14159265f

class QuantizedGradientChannel
{
public:
	int enabled; //if true enable gradient histogram channels
	int binSize; //spatial bin size (defaults to shrink)
	int orientationChannels; //number of orientation channels
	int useSoftBinning; //if true use "soft" bilinear spatial binning
	int useHogNormalization; //if true perform 4-way hog normalization/clipping
	double clipHog; //value at which to clip hog histogram bins
	int nChannels;
	cv::String padWith;

	void readGradientHistogram(cv::FileNode);
	//comes from gradientMex.cpp, the return still needs a type
	std::vector<cv::Mat> mGradHist(cv::Mat, cv::Mat, int);
};

#endif
