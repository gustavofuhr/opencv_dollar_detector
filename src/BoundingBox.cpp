#include "BoundingBox.h"

void BoundingBox::plot (cv::Mat &frame, cv::Scalar color)
{
	cv::Point br;
	br.x = topLeftPoint.x + width;
	br.y = topLeftPoint.y + height;
	cv::rectangle(frame, topLeftPoint, br, color, 2.0);
}

std::string BoundingBox::toString (int frameIndex)
{
	std::ostringstream result;

	result << frameIndex << " " << topLeftPoint.x << " " << topLeftPoint.y << " " << height << " " << width << " " << score;

	return result.str();
}



cv::Point2i  BoundingBox::get_bottomRightPoint() {

	cv::Point2i br;
	br.x = topLeftPoint.x + width;
	br.y = topLeftPoint.y + height;

	return br;

}