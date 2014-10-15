#ifndef DETECTOR_H
#define DETECTOR_H

#include "Options.h"
#include "Clf.h"
#include "Info.h"
#include "BoundingBox.h"
#include "utils.h"

class Detector
{
public:
	Options opts; //opts contains the Pyramid
	Clf clf;

	BB_Array detections;
	// debug
	int totalDetections;

	void exportDetectorModel(cv::String);
	void importDetectorModel(cv::String);
	void acfDetect(cv::Mat);
	BB_Array bbNms(BB_Array bbs);
	BB_Array nmsMax(BB_Array source, bool greedy);
};

#endif
