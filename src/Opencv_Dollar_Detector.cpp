#include <sstream>
#include "Opencv_Dollar_Detector.h"

// valgrind --tool=memcheck --leak-check=yes --log-file=valgrind.log ./opencv_dollar_detector ../opencv_dollar_detector/detector.xml ../datasets/small

// call: ./opencv_dollar_detector ../opencv_dollar_detector/odd.conf
int main(int argc, char *argv[]) 
{
	if (argc < 2)
	{
		std::cout << " # Argument Error: this program requires a conf file." << std::endl;
		return 1;
	}
	else
	{
		clock_t start = clock();
		int firstFrame=0, lastFrame=666666666;
		
		OddConfig odd_config(argv[1]);
		Detector d(odd_config);

		// loads all detector settings from the provided xml file
		d.importDetectorModel(odd_config.detectorFileName);

		cv::Mat H(3, 3, CV_32F);
		H.at<float>(0,0) = odd_config.calibrationP->at<float>(0,0);
		H.at<float>(0,1) = odd_config.calibrationP->at<float>(0,1);
		H.at<float>(0,2) = odd_config.calibrationP->at<float>(0,3);
		H.at<float>(1,0) = odd_config.calibrationP->at<float>(1,0);
		H.at<float>(1,1) = odd_config.calibrationP->at<float>(1,1);
		H.at<float>(1,2) = odd_config.calibrationP->at<float>(1,3);
		H.at<float>(2,0) = odd_config.calibrationP->at<float>(2,0);
		H.at<float>(2,1) = odd_config.calibrationP->at<float>(2,1);
		H.at<float>(2,2) = odd_config.calibrationP->at<float>(2,3);

		cv::Point groundPoint = imagePoint2groundPlanePoint(0.0, 50.0, 1.0, H);

		std::cout << "before find, x=" << groundPoint.x << ", y=" << groundPoint.y << std::endl;

		cv::Point imagePoint = worldPoint2imagePoint(groundPoint.x, groundPoint.y, 1.0, H);

		std::cout << "before find, u=" << imagePoint.x << ", v=" << imagePoint.y << std::endl;

		float height = findWorldHeight(*(odd_config.calibrationP), 0.0, 0.0, groundPoint.x, groundPoint.y);
		std::cout << "worldHeight=" << height << std::endl;
		std::cin.get();

		// gets names for all the files inside the data set folder
		std::vector<std::string> imageNames = getDataSetFileNames(odd_config.dataSetDirectory);

		// apply the detection on all images
		d.acfDetect(imageNames, odd_config.dataSetDirectory, odd_config.firstFrame, odd_config.lastFrame);

		clock_t end = clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;

		std::cout << "\nTotal processing time was " << elapsed_secs << " seconds.\n";
		std::cout << "Time elapsed calculating features: " << d.opts.pPyramid.totalTimeForRealScales << std::endl;
		std::cout << "Time elapsed approximating features: " << d.opts.pPyramid.totalTimeForApproxScales << std::endl;
		std::cout << "Time elapsed during detections: " << d.timeSpentInDetection << std::endl;

		return 0;
	}
}
