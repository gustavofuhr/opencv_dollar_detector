#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H

#include <cv.h>

class BoundingBox
{
public:
	cv::Point    topLeftPoint;
	cv::Point2i  get_bottomRightPoint();


	int height;
	int width;
	float score;
	int scale;
	float world_height;

	void plot (cv::Mat &frame, cv::Scalar color);
	std::string toString (int frameIndex);

	// this was added to be able to sort BB_Array objects
	bool operator< (const BoundingBox &other) const 
	{
		//return score < other.score;
  		return score > other.score;
 	}

};

typedef std::vector<BoundingBox> BB_Array;
typedef std::vector<BB_Array> BB_Array_Array;

#endif