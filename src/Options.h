#ifndef OPTIONS_H
#define OPTIONS_H

#include "Pyramid.h"

class Options
{
public:
	Pyramid pPyramid;
	int modelDs[2];
	int modelDsPad[2];
	int stride;
	int cascadeThreshold;

	int suppressionType;
	double overlapArea; // area of overlap for bbs
	cv::String overlapDenominator; // area of overlap denominator ('union' or 'min')
	float suppressionThreshold;

	void readOptions(cv::FileNode);
};

#endif
