#include "opencv/cv.h"

using namespace cv;

class BoundingBox
{
public:
	Point firstPoint;
	int height;
	int width;
};