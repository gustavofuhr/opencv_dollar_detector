#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	Mat image;
	Detector d;

	//call to read the detector file
	//d.readDetectorModel(DETECTOR_FILE);
	//d.readDetectorModel("../opencv_dollar_detector/detector.xml");

	image = imread("../opencv_dollar_detector/frame0254.png");
	//BoundingBox* bbs = d.acfDetect(image);

	return 0;
}
