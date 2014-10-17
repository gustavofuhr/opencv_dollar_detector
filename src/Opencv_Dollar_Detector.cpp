#include "Opencv_Dollar_Detector.h"

int main() 
{
	Detector d;
	Info result;

	cv::String fileString = "../opencv_dollar_detector/detector.xml";
	//cv::String fileString = "../opencv_dollar_detector/calculate_all_scales_detector.xml";
	d.importDetectorModel(fileString);
	
	cv::Mat image = cv::imread("../opencv_dollar_detector/frame0254.png");
	cv::Mat image2 = cv::imread("../opencv_dollar_detector/person_011.bmp");
	cv::Mat image3 = cv::imread("../opencv_dollar_detector/per00001.ppm");
	
	std::vector<cv::Mat> dataSet;
	dataSet.push_back(image);
	dataSet.push_back(image2);
	dataSet.push_back(image3);	

	d.acfDetect(dataSet);

	return 0;
}
