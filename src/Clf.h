#ifndef CLF_H
#define CLF_H

#include "cv.h"

using namespace cv;

class Clf
{
public:
	Mat fids;
	Mat thrs;
	Mat child;
	Mat hs;
	Mat weights;
	Mat depth;
	Mat errs;
	Mat losses;
	int treeDepth;

	void readClf(FileNode);
};

#endif
