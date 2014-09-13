#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	cv::Mat image;
	Detector d;
	Info result;
	//cv::FileStorage file("some_name.ext", cv::FileStorage::WRITE);

	cv::String fileString = "../opencv_dollar_detector/detector.xml";
	//cv::String fileString = "../opencv_dollar_detector/calculate_all_scales_detector.xml";
	d.importDetectorModel(fileString);
	
	image = cv::imread("../opencv_dollar_detector/frame0254.png");

	d.acfDetect(image);

	// this will print 0, so nothing is being detected
	//d.opts.pPyramid.debugWindow("number of detections",d.totalDetections);

	return 0;
}
