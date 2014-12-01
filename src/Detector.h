#ifndef DETECTOR_H
#define DETECTOR_H

#include "Options.h"
#include "Info.h"
#include "BoundingBox.h"
#include "utils.h"

class Detector
{
public:
	Options opts; //opts contains the Pyramid

	//Clf clf;
	cv::Mat fids;
	cv::Mat thrs;
	cv::Mat child;
	cv::Mat hs;
	cv::Mat weights;
	cv::Mat depth;
	cv::Mat errs;
	cv::Mat losses;
	int treeDepth;

	BB_Array_Array detections;

	double timeSpentInDetection;

	void exportDetectorModel(cv::String);
	void importDetectorModel(cv::String);
	BB_Array applyDetectorToFrame(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, 
									float *hs, uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns);
	void acfDetect(std::vector<std::string> imageNames, std::string dataSetDirectoryName, int firstFrame, int lastFrame);
	BB_Array nonMaximalSuppression(BB_Array bbs);


private:
	BB_Array generateCandidates(int imageHeight, int imageWidth, cv::Mat_<float> &P, 
							float meanHeight = 1700, float stdHeight = 100, float factorStdHeight = 2.0);

	int findClosestScaleFromBbox(std::vector<Info> &pyramid, BoundingBox &bb,
												int modelHeight, int imageHeight);


	BB_Array applyDetectorToFrameSmart(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, float *hs, 
										uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns, int imageWidth, int imageHeight, 
										cv::Mat_<float> &P, cv::Mat &debug_image);

	void bbTopLeft2PyramidRowColumn(int *r, int *c, BoundingBox &bb, int modelHt, int modelWd, int ith_scale, int stride);

};

#endif
