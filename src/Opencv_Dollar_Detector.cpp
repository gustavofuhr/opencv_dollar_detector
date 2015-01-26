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
		
		OddConfig settings(argv[1]);
		Detector d(settings);

		// loads all detector settings from the provided xml file
		d.importDetectorModel(settings.detectorFileName);

		/*
		// debug: tests new functions
		std::vector<cv::Point2f> points = findGroundPlaneAndImageIntersectionPoints(768*1.5, 576*1.5, 41, 100, 2500.0, *(settings.projectionMatrix), *(settings.homographyMatrix));
		//std::vector<cv::Point2f> trimmedPoints = trimGroundPlanesBottomPoints(768*1.5, 576*1.5, 41, 100, 2000, points, *(settings.projectionMatrix), *(settings.homographyMatrix));
		int octaves = findNecessaryNumberOfOctaves(768*1.5, 576*1.5, 41, 100, 1500.0, 2100.0, *(settings.projectionMatrix), *(settings.homographyMatrix));
		std::cout << "octaves=" << octaves << std::endl;
		// bottom left corner
		double scaleBL = findLastNecessaryScaleInAPoint(0, 576*1.5, 576*1.5, 100, 2100.0, *(settings.projectionMatrix), *(settings.homographyMatrix));
		// bottom center
		double scaleBC = findLastNecessaryScaleInAPoint(768*1.5/2, 576*1.5, 576*1.5, 100, 2100.0, *(settings.projectionMatrix), *(settings.homographyMatrix));
		// bottom right corner
		double scaleBR = findLastNecessaryScaleInAPoint(768*1.5, 576*1.5, 576*1.5, 100, 2100.0, *(settings.projectionMatrix), *(settings.homographyMatrix));
		std::cin.get();
		// debug */

		// gets names for all the files inside the data set folder
		std::vector<std::string> imageNames = getDataSetFileNames(settings.dataSetDirectory);

		if (settings.lastFrame <= settings.firstFrame)
			settings.lastFrame = imageNames.size();
		
		// apply the detection on all images
		d.acfDetect(imageNames, settings.dataSetDirectory, settings.firstFrame, settings.lastFrame);

		clock_t end = clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;

		std::cout << "\nTotal processing time was " << elapsed_secs << " seconds.\n";
		std::cout << "Time elapsed calculating features: " << d.opts.pPyramid.totalTimeForRealScales << std::endl;
		std::cout << "Time elapsed approximating features: " << d.opts.pPyramid.totalTimeForApproxScales << std::endl;
		std::cout << "Time elapsed during detections: " << d.timeSpentInDetection << std::endl;

		return 0;
	}
}
