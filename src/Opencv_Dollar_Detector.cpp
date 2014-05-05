#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	Mat image;
	Detector d;

	//call to read the detector file
	//d.importDetectorModel(DETECTOR_FILE);
	String fileString = "../opencv_dollar_detector/detector.xml";
	d.importDetectorModel(fileString);
	
	image = imread("../opencv_dollar_detector/frame0254.png");
	
	//imshow("main",image);

	BoundingBox** bbs = d.acfDetect(image);

	waitKey();

	return 0;
}
