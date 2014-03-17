#include "Channel.h"

class Pyramid
{
	public:
        Channel pChns; //parameters for creating channels
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
        int computedChannels[]; //[nScales x nTypes] cell array of computed channels


        Pyramid computeMultiScaleChannelFeaturePyramid(Mat);
        Pyramid computeSingleScaleChannelFeaturePyramid(Mat, Channel);
};
