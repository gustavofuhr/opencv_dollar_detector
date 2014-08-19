#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
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

	BB_Array_Array detections;
	// debug
	int totalDetections;

	void exportDetectorModel(cv::String);
	void importDetectorModel(cv::String);
	void getChild(float*, uint32_t*, uint32_t*,	float*, uint32_t, uint32_t&, uint32_t&);
	void acfDetect(cv::Mat);
	BB_Array bbNms(BB_Array bbs, int size);
	BB_Array nmsMax (BB_Array source, int size, bool greedy);
};

#endif
