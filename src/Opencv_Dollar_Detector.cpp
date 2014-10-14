#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	cv::Mat image;
	Detector d;
	Info result;

	cv::String fileString = "../opencv_dollar_detector/detector.xml";
	//cv::String fileString = "../opencv_dollar_detector/calculate_all_scales_detector.xml";
	d.importDetectorModel(fileString);
	
	image = cv::imread("../opencv_dollar_detector/frame0254.png");

	d.acfDetect(image);

	return 0;
}
