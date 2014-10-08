#include <cv.h>

// this class belongs to opts in the detector model but was empty
class PNms
{
public:
	cv::String type;
	double overlap; // area of overlap for bbs
	cv::String ovrDnm; // area of overlap denominator ('union' or 'min')
	int threshold;

	void readPNms(cv::FileNode);
};
