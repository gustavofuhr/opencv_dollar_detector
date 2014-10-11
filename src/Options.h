#ifndef OPTIONS_H
#define OPTIONS_H

#include "Pyramid.h"
#include "PNms.h"
#include "PBoost.h"

class Options
{
public:
	Pyramid pPyramid;
	int modelDs[2];
	int modelDsPad[2];
	PNms pNms;
	int stride;
	int cascadeThreshold;
	double cascadeCalibration;
	int nWeak[4];
	PBoost pBoost;
	int seed;
	cv::String name;
	cv::String posGtDir;
	cv::String posImgDir;
	cv::String negImgDir;
	cv::String posWinDir;
	cv::String negWinDir;
	//function imreadf;
	cv::Mat imreadp;
	double pLoadSquarify[2];
	int nPos;
	int nNeg;
	int nPerNeg;
	int nAccNeg;
	int pJitterFlip;
	int winsSave;

	void readOptions(cv::FileNode);
};

#endif
