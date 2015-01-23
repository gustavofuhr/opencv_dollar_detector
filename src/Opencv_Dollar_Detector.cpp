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
		
		OddConfig odd_config(argv[1]);
		Detector d(odd_config);

		// loads all detector settings from the provided xml file
		d.importDetectorModel(odd_config.detectorFileName);

		/*
		// debug: tests new functions
		std::vector<cv::Point2f> points = findGroundPlaneAndImageIntersectionPoints(768*1.5, 576*1.5, 41, 100, 2500.0, *(odd_config.projectionMatrix), *(odd_config.homographyMatrix));
		//std::vector<cv::Point2f> trimmedPoints = trimGroundPlanesBottomPoints(768*1.5, 576*1.5, 41, 100, 2000, points, *(odd_config.projectionMatrix), *(odd_config.homographyMatrix));
		int octaves = findNecessaryNumberOfOctaves(768*1.5, 576*1.5, 41, 100, 1500.0, 2100.0, *(odd_config.projectionMatrix), *(odd_config.homographyMatrix));
		std::cout << "octaves=" << octaves << std::endl;
		// bottom left corner
		double scaleBL = findLastNecessaryScaleInAPoint(0, 576*1.5, 576*1.5, 100, 2100.0, *(odd_config.projectionMatrix), *(odd_config.homographyMatrix));
		// bottom center
		double scaleBC = findLastNecessaryScaleInAPoint(768*1.5/2, 576*1.5, 576*1.5, 100, 2100.0, *(odd_config.projectionMatrix), *(odd_config.homographyMatrix));
		// bottom right corner
		double scaleBR = findLastNecessaryScaleInAPoint(768*1.5, 576*1.5, 576*1.5, 100, 2100.0, *(odd_config.projectionMatrix), *(odd_config.homographyMatrix));
		std::cin.get();
		// debug */


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
