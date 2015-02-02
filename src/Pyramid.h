#ifndef PYRAMID_H
#define PYRAMID_H

#include "ChannelFeatures.h"
#include "Info.h"
#include "utils.h"

class Pyramid
{
public:
	ChannelFeatures pChns; //parameters for creating channels
	int scalesPerOctave; //number of scales per octave, substitutes nPerOct
	int upsampledOctaves; //number of upsampled octaves to compute, substitutes nOctUp
	int approximatedScales; //number of approx. scales
	int providedLambdas;
	double lambdas[3]; //for power law scaling
	int padSize;
	std::vector<int> pad; //amount to pad channels (along T/B and L/R)
	int minImgSize[2]; //minimum image size for channel computation, substitutes minDs
	float smoothRadius; //radius for channel smoothing (using convTri)
	int concatenateChannels; //if true, concatenate channels
	int completeInput; //if true does not check/set default vals in pPyramid

	// output attributes
	int computedScales; //number of scales computed
	std::vector<double> scales;
	std::vector<double> scales_w;
	std::vector<double> scales_h;

	// time attributes
	double totalTimeForRealScales;
	double totalTimeForApproxScales;

	void readPyramid(cv::FileNode);
	std::vector<Info> computeFeaturePyramid(cv::Mat I, bool useCalibration);
	Info computeSingleScaleChannelFeatures(float* I, int rows, int cols);
	cv::Mat TriangleFilterConvolution(cv::Mat, int, int, int);
	void getScales(int, int , int);
	void calibratedGetScales(int h, int w, int shrink, int boundingBoxImageWidth, int boundingBoxImageHeight, double maxPedestrianWorldHeight, 
	cv::Mat_<float> &projection, cv::Mat_<float> &homography);
};

#endif