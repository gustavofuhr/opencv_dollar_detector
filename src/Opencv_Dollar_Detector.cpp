#include "Opencv_Dollar_Detector.h"

// valgrind --tool=memcheck --leak-check=yes --log-file=valgrind.log ./opencv_dollar_detector ../opencv_dollar_detector/detector.xml ../datasets/small


// call: ./opencv_dollar_detector ../opencv_dollar_detector/detector.xml ../datasets/small
int main(int argc, char *argv[]) 
{
	if (argc < 2)
	{
		std::cout << " # Argument Error: this program expects at least two arguments (detector file name and data set directory)." << std::endl;
		return 1;
	}
	else
	{
		clock_t start = clock();
		int firstFrame=0, lastFrame=666666666;

		// reads index of the first and last frames 
		if (argc > 3)
		{
			firstFrame = atoi(argv[3]);
		}
		if (argc > 4)
		{
			lastFrame = atoi(argv[4]);
		}

		Detector d;

		/*
		// experimental
		// reads homography matrix from xml file
		std::vector<cv::Mat> P_and_H = readProjectionAndHomographyFromCalibrationFile("../opencv_dollar_detector/towncentre_calib.xml");
		cv::Mat projection = P_and_H[0];
		cv::Mat homography = P_and_H[1];

		cv::Point groundPoint = imagePoint2groundPlanePoint(0, 50, 1, homography);
		float boundingBoxWorldHeight = findWorldHeight(projection, 0, 0, groundPoint.x, groundPoint.y); 
		// experimental */

		// loads all detector settings from the provided xml file
		cv::String detectorFileName = argv[1];
		d.importDetectorModel(detectorFileName);

		// gets names for all the files inside the data set folder
		std::string dataSetDirectory = argv[2];
		std::vector<std::string> imageNames = getDataSetFileNames(dataSetDirectory);

		// apply the detection on all images
		d.acfDetect(imageNames, dataSetDirectory, firstFrame, lastFrame);

		clock_t end = clock();
		double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;

		std::cout << "\nTotal processing time was " << elapsed_secs << " seconds.\n";
		std::cout << "Time elapsed calculating features: " << d.opts.pPyramid.totalTimeForRealScales << std::endl;
		std::cout << "Time elapsed approximating features: " << d.opts.pPyramid.totalTimeForApproxScales << std::endl;
		std::cout << "Time elapsed during detections: " << d.timeSpentInDetection << std::endl;

		return 0;
	}
}
