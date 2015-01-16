#ifndef DETECTOR_H
#define DETECTOR_H

#include "Options.h"
#include "Info.h"
#include "BoundingBox.h"
#include "utils.h"
#include <fstream>
#include <string>


struct OddConfig {
	float resizeImage;
	std::string detectorFileName;
	std::string dataSetDirectory;
	int firstFrame, lastFrame;

	bool displayDetections, saveFrames, saveLog; 
	std::string outputFolder;
	std::string logFilename;

	bool useCalibration;

	cv::Mat_<float> *projectionMatrix;
	cv::Mat_<float> *homographyMatrix;

	float supressionThreshold;

	OddConfig(std::string config_file);
};


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
	OddConfig config;

	BB_Array_Array detections;

	double timeSpentInDetection;

	void exportDetectorModel(cv::String);
	void importDetectorModel(cv::String);
	BB_Array applyDetectorToFrame(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, 
									float *hs, uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns);
	void acfDetect(std::vector<std::string> imageNames, std::string dataSetDirectoryName, int firstFrame, int lastFrame);
	BB_Array nonMaximalSuppression(BB_Array bbs);

	Detector(OddConfig _config): config(_config) { };

private:
	BoundingBox pyramidRowColumn2BoundingBox(int r, int c,  int modelHt, int modelWd, int ith_scale, int stride);

	//BB_Array *generateCandidates(int imageHeight, int imageWidth, cv::Mat_<float> &P, double *maxHeight, float meanHeight = 1800, float stdHeight= 100, 
	//float factorStdHeight = 2.0, float groundPlaneMinX, float groundPlaneMaxX, float groundPlaneMinY, float groundPlaneMaxY) ;

	BB_Array *generateCandidates(int imageHeight, int imageWidth, float groundPlaneMinX, float groundPlaneMaxX, float groundPlaneMinY, float groundPlaneMaxY, 
		cv::Mat_<float> &P, double *maxHeight, float meanHeight = 1800, float stdHeight = 100, float factorStdHeight = 2.0);

	BB_Array* generateCandidatesFaster(int imageHeight, int imageWidth, int shrink, cv::Mat_<float> &P, double *maxHeight,
							cv::Mat &im_debug, float meanHeight = 1800, float stdHeight = 100, float factorStdHeight = 2.0);

	int findClosestScaleFromBbox(std::vector<Info> &pyramid, BoundingBox &bb,
												int modelHeight, int imageHeight);

	int findClosestScaleFromBbox2(std::vector<Info> &pyramid, BoundingBox &bb,
												int modelHeight, double shrink);


	BB_Array applyCalibratedDetectorToFrame(std::vector<Info> pyramid, BB_Array* bbox_candidates, int shrink, int modelHt, int modelWd, int stride, float cascThr, 
											float *thrs, float *hs, uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns, int imageWidth, 
											int imageHeight, cv::Mat_<float> &P, cv::Mat &debug_image);

	void bbTopLeft2PyramidRowColumn(int *r, int *c, BoundingBox &bb, int modelHt, int modelWd, int ith_scale, int stride);
	BB_Array nonMaximalSuppressionSmart(BB_Array bbs, double meanHeight, double stdHeight);

};

#endif
