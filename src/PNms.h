#include "cv.h"

using namespace cv;

// this class belongs to opts in the detector model but was empty
class PNms
{
public:
	String type;
	double overlap; // area of overlap for bbs
	String ovrDnm; // area of overlap denominator ('union' or 'min')
	int threshold;

	void readPNms(FileNode);
};
