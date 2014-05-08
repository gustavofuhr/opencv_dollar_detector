#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	cv::Mat image;
	Detector d;

	cv::String fileString = "../opencv_dollar_detector/detector.xml";
	d.importDetectorModel(fileString);
	
	image = cv::imread("../opencv_dollar_detector/frame0254.png");
	
	//imshow("main",image);

	//still have a segmentation fault here
	//but first i need to work a bit more on acfDetect's return
	//BB_Array* bbs = (BB_Array*)malloc(100 * sizeof(BB_Array));

	//bbs = d.acfDetect(image);
	d.acfDetect(image);

	BB_Array test = d.detections[0];
	// this next line generates a floating point error
	//BoundingBox test2 = test[1];
	
	// this will print 0, so nothing is being detected
	d.opts.pPyramid.debugWindow("done",d.totalDetections);

	return 0;
}
