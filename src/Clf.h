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
	Mat depth;
   Mat errs;
	Mat losses;
	int treeDepth;
};

#endif
