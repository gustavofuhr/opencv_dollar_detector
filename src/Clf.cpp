#include "Clf.h"

void Clf::readClf(cv::FileNode clfNode)
{
	clfNode["fids"] >> fids;
	clfNode["thrs"] >> thrs;
	clfNode["child"] >> child;
	clfNode["hs"] >> hs;
	clfNode["weights"] >> weights;
	clfNode["depth"] >> depth;
	clfNode["errs"] >> errs;
	clfNode["losses"] >> losses;		
	clfNode["treeDepth"] >> treeDepth;
}
