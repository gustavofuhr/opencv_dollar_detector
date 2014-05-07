//#include "cv.h"

class BoundingBox
{
public:
	cv::Point firstPoint;
	int height;
	int width;
	int score;
};

typedef std::vector<BoundingBox> BB_Array;
typedef std::vector<BB_Array> BB_Array_Array;
