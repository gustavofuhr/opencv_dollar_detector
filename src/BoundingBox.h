//#include "cv.h"

class BoundingBox
{
public:
	cv::Point firstPoint;
	int height;
	int width;
	int score;

	// this was added to be able to sort BB_Array objects
	bool operator< (const BoundingBox &other) const 
	{
		//return score < other.score;
  	return score > other.score;
  }
};

typedef std::vector<BoundingBox> BB_Array;
typedef std::vector<BB_Array> BB_Array_Array;
