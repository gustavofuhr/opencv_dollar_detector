#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	cv::Mat image;
	Detector d;
	Info result;
	cv::FileStorage file("some_name.ext", cv::FileStorage::WRITE);

	cv::String fileString = "../opencv_dollar_detector/detector.xml";
	d.importDetectorModel(fileString);
	
	image = cv::imread("../opencv_dollar_detector/frame0254.png");
	
	//imshow("main",image);

	//still have a segmentation fault here
	//but first i need to work a bit more on acfDetect's return
	//BB_Array* bbs = (BB_Array*)malloc(100 * sizeof(BB_Array));

	//bbs = d.acfDetect(image);
	//d.acfDetect(image);
	
	d.opts.pPyramid.debugWindow("blah",0);

	result = d.opts.pPyramid.computeSingleScaleChannelFeatures(image);

	d.opts.pPyramid.debugWindow("blah",1);

	file << result.image;
	//file << result.gradientMagnitude;
	//file << result.gradientHistogram;

	//BB_Array test = d.detections[0];
	
	// this will print 0, so nothing is being detected
	//d.opts.pPyramid.debugWindow("number of detections",d.totalDetections);

	return 0;
}
