#ifndef CLF_H
#define CLF_H

#include <cv.h>

class Clf
{
public:
	cv::Mat fids;
	cv::Mat thrs;
	cv::Mat child;
	cv::Mat hs;
	cv::Mat weights;
	cv::Mat depth;
	cv::Mat errs;
	cv::Mat losses;
	int treeDepth;

	void readClf(cv::FileNode);
};

#endif
