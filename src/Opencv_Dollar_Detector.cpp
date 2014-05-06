#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	Mat image;
	Detector d;

	//call to read the detector file
	//d.importDetectorModel(DETECTOR_FILE);
	String fileString = "../opencv_dollar_detector/detector.xml";
	d.importDetectorModel(fileString);
	
	image = imread("../opencv_dollar_detector/frame0254.png");
	
	//imshow("main",image);

	//still have a segmentation fault here
	//but first i need to work a bit more on acfDetect's return
	//BoundingBox** bbs;
	//bbs = d.acfDetect(image);
	d.acfDetect(image);

	d.opts.pPyramid.debugWindow("done",0);

	return 0;
}
