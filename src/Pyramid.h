#include "ChannelFeatures.h"
#include "Info.h"

class Pyramid
{
	public:
		ChannelFeatures pChns; //parameters for creating channels
		int scalesPerOctave; //number of scales per octave
		int upsampledOctaves; //number of upsampled octaves to compute
		int approximatedScales; //number of approx. scales
		double lambdas[3]; //for power law scaling
		int pad[2]; //amount to pad channels (along T/B and L/R)
		int minImgSize[2]; //minimum image size for channel computation
		int smoothRadius; //radius for channel smoothing (using convTri)
		int concatenateChannels; //if true, concatenate channels
		int completeInput; //if true does not check/set default vals in pPyramid

		//output parameters
		int channelTypes; //number of channel types
		int computedScales; //number of scales computed
		Info computedChannels[]; //[nScales x nTypes] cell array of computed channels
		double* scales;
		Point* scaleshw;

		void readPyramid(FileNode);
		Pyramid computeMultiScaleChannelFeaturePyramid(Mat);
		Info computeSingleScaleChannelFeatures(Mat);
		Mat TriangleFilterConvolution(Mat I, int r, int s, int nomex);
		void getScales(int h, int w, int shrink);
};
