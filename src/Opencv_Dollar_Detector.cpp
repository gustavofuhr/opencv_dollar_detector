#include "Opencv_Dollar_Detector.h"

int main(int argc, char** argv) 
{
	Mat image;
	Detector d;

	//call to read the detector file
	//d.importDetectorModel(DETECTOR_FILE);
	String fileString = "../opencv_dollar_detector/detector.xml";
	d.importDetectorModel(fileString);
	//these next two matices seem to be okay
	cout << "fids00: " << d.clf.fids.at<int>(0,0) << "\n";
	cout << "thrs00: " << d.clf.thrs.at<double>(0,0) << "\n";

	//problem reading the child matrix
	cout << "child00: " << d.clf.child.at<int>(1,0) << "\n";
	
	//this one seems to be okay also
	cout << "hs00:" << d.clf.hs.at<double>(0,0) << "\n";

	cout << "weights00:" << d.clf.weights.at<double>(0,0) << "\n";

	image = imread("../opencv_dollar_detector/frame0254.png");
	//BoundingBox* bbs = d.acfDetect(image);
	
	//namedWindow("main");
	imshow("main",image);

	waitKey();

	return 0;
}
