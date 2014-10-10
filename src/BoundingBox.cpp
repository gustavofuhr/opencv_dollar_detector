#include "BoundingBox.h"

void BoundingBox::plot (cv::Mat &frame, cv::Scalar color)
{
	cv::Point br;
	br.x = firstPoint.x + width;
	br.y = firstPoint.y + height;
	cv::rectangle(frame, firstPoint, br, color, 2.0);
}