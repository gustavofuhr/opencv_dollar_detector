#ifndef DETECTOR_H
#define DETECTOR_H

#include "Options.h"
#include "Info.h"
#include "BoundingBox.h"
#include "utils.h"

class Detector
{
public:
	Options opts; //opts contains the Pyramid

	//Clf clf;
	cv::Mat fids;
	cv::Mat thrs;
	cv::Mat child;
	cv::Mat hs;
	cv::Mat weights;
	cv::Mat depth;
	cv::Mat errs;
	cv::Mat losses;
	int treeDepth;

	BB_Array detections;

	void exportDetectorModel(cv::String);
	void importDetectorModel(cv::String);
	void acfDetect(cv::Mat);
	BB_Array nonMaximalSuppression(BB_Array bbs);
};

#endif
