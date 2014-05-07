#include "ChannelFeatures.h"

void ChannelFeatures::readChannelFeatures(cv::FileNode chFeatNode)
{
	pColor.readColorChannel(chFeatNode["pColor"]);
	
	pGradMag.readGradientMagnitude(chFeatNode["pGradMag"]);

	pGradHist.readGradientHistogram(chFeatNode["pGradHist"]);	

	shrink = chFeatNode["shrink"];
	complete = chFeatNode["complete"];

}
