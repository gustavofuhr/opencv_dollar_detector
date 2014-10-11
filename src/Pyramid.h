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
	int *pad; //amount to pad channels (along T/B and L/R)
	int minImgSize[2]; //minimum image size for channel computation, substitutes minDs
	int smoothRadius; //radius for channel smoothing (using convTri)
	int concatenateChannels; //if true, concatenate channels
	int completeInput; //if true does not check/set default vals in pPyramid

	//output parameters
	int channelTypes; //number of channel types
	int computedScales; //number of scales computed
	Info* computedChannels; //[nScales x nTypes] cell array of computed channels
	double* scales;
	double* scales_w;
	double* scales_h;

	void readPyramid(cv::FileNode);
	void computeMultiScaleChannelFeaturePyramid(cv::Mat);
	Info computeSingleScaleChannelFeatures(cv::Mat);
	cv::Mat TriangleFilterConvolution(cv::Mat, int, int, int);
	void getScales(int, int , int);

	// debug
	void readScalesFromXML(cv::FileNode pyramid);
};
