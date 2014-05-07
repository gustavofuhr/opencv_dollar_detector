#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "Options.h"
#include "Clf.h"
#include "Info.h"
#include "BoundingBox.h"

class Detector
{
public:
	Options opts; //opts contains the Pyramid
	Clf clf;

	BB_Array_Array detections;

	void exportDetectorModel(cv::String);
	void importDetectorModel(cv::String);
	void getChild(float*, uint32_t*, uint32_t*,	float*, uint32_t, uint32_t&, uint32_t&);
	void acfDetect(cv::Mat);
	BB_Array bbNms(BoundingBox* bbs, int size);
	BB_Array nmsMax (BB_Array source, int size, bool greedy);
};

#endif
