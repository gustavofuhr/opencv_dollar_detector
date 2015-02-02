#include <sstream>
#include "Opencv_Dollar_Detector.h"

// memory check: 	valgrind --tool=memcheck --leak-check=full --log-file=valgrind.log --error-limit=no ./opencv_dollar_detector ../opencv_dollar_detector/odd.conf
// profiler:		valgrind --tool=callgrind ./opencv_dollar_detector ../opencv_dollar_detector/odd.conf
// call: 			./opencv_dollar_detector ../opencv_dollar_detector/odd.conf
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

		// gets names for all the files inside the data set folder
		std::vector<std::string> imageNames = getDataSetFileNames(settings.dataSetDirectory);

		/*
		// debug: tests new functions
		cv::Mat I = cv::imread(settings.dataSetDirectory + '/' + imageNames[0]);
		cv::Mat image;
		cv::normalize(I, image, 0.0, 1.0, cv::NORM_MINMAX, CV_32FC3);
		cv::imshow("original", image);
		float* floatImage = (float*) malloc(image.rows*image.cols*3*sizeof(float));
		cvMat2floatArray(image, floatImage, 3);
		cv::Mat convertedImage = floatArray2cvMat(floatImage, image.rows, image.cols, 3);
		cv::imshow("converted", convertedImage);

		cv::Mat gray_image;
 		cv::cvtColor(image, gray_image, CV_BGR2GRAY );
 		float* floatGray = (float*) malloc(image.rows*image.cols*sizeof(float));
		cvMat2floatArray(gray_image, floatGray, 1);
		cv::Mat convertedGray = floatArray2cvMat(floatGray, image.rows, image.cols, 1);
		cv::imshow("convertedGray", convertedGray);
		
		cv::waitKey();
		// debug */

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

		delete settings.projectionMatrix;
		delete settings.homographyMatrix;

		return 0;
	}
}
