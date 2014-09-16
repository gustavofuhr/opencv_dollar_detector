#include "Clf.h"

void Clf::readClf(cv::FileNode clfNode)
{
	clfNode["fids"] >> fids;
	clfNode["child"] >> child;
	clfNode["thrs"] >> thrs;
	clfNode["hs"] >> hs;
	clfNode["weights"] >> weights;
	clfNode["depth"] >> depth;
	clfNode["errs"] >> errs;
	clfNode["losses"] >> losses;		
	clfNode["treeDepth"] >> treeDepth;
}
