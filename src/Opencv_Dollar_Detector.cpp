#include "Opencv_Dollar_Detector.h"

int main() 
{
	cv::Mat image;
	Detector d;
	Info result;

	cv::String fileString = "../opencv_dollar_detector/detector.xml";
	//cv::String fileString = "../opencv_dollar_detector/calculate_all_scales_detector.xml";
	d.importDetectorModel(fileString);
	
	image = cv::imread("../opencv_dollar_detector/frame0254.png");
	// image = cv::imread("../opencv_dollar_detector/person_011.bmp");

	d.acfDetect(image);

	return 0;
}
