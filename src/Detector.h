#ifndef DETECTOR_H
#define DETECTOR_H

#include <vector>
#include "Options.h"
#include "Clf.h"
#include "Info.h"
#include "BoundingBox.h"

using namespace cv;
using namespace std;

typedef std::vector<BoundingBox> BB_Array;

class Detector
{
public:
	Options opts; //opts contains the Pyramid
	Clf clf;

	void readColorChannel(FileNode);	
	void readChannelFeatures(FileNode);
	void readPyramid(FileNode);
	void readOptions(FileNode);
	void readPNms(FileNode);
	void readInfo(FileNode);
	void readGradientMagnitude(FileNode);
	void readGradientHistogram(FileNode);
	void exportDetectorModel(String);
	void importDetectorModel(String);
	void getChild(float*, uint32_t*, uint32_t*,	float*, uint32_t, uint32_t&, uint32_t&);
	BB_Array* acfDetect(Mat);
	BB_Array bbNms(BoundingBox* bbs, int size);
	BB_Array nmsMax (BB_Array source, int size, bool greedy);
};

#endif
