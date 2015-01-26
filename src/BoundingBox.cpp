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

	result << "Frame: " << std::setw(4) << frameIndex << ", TopLeftPoint: (" << std::setw(4) << topLeftPoint.x << "," << std::setw(4) << topLeftPoint.y << 
			"), Height: " << std::setw(4) << height << ", Width: " << std::setw(3) << width << ", Score: " << score << std::endl;

	return result.str();
}

cv::Point2i  BoundingBox::get_bottomRightPoint() {

	cv::Point2i br;
	br.x = topLeftPoint.x + width;
	br.y = topLeftPoint.y + height;

	return br;

}