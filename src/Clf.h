#ifndef CLF_H
#define CLF_H

class Clf
{
public:
   Mat fids;
   Mat thrs;
	Mat child;
	Mat hs;
	Mat weights;
	unsigned int depth[7][2048];
    Mat errs;
	double losses[2048];
	int treeDepth;
};

#endif
