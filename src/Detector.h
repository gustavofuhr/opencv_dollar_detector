#ifndef DETECTOR_H
#define DETECTOR_H

#include "Options.h"
#include "Clf.h"
#include "Info.h"
#include "BoundingBox.h"

using namespace cv;
using namespace std;

class Detector
{
public:
	Options opts; //opts contains the Pyramid
	Clf clf;

	void readColorChannel(FileNode);	
	void readChannelFeatures(FileNode);
	void readPyramid(FileNode);
	void readOptions(FileNode);
	void readInfo(FileNode);
	void readGradientMagnitude(FileNode);
	void readGradientHistogram(FileNode);
	void importDetectorModel(String);
	void getChild(float*, uint32_t*, uint32_t*,	float*, uint32_t, uint32_t&, uint32_t&);
	BoundingBox* acfDetect(Mat);
	
};

#endif
