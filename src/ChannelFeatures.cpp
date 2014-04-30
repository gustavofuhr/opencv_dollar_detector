#include "ChannelFeatures.h"

void ChannelFeatures::readChannelFeatures(FileNode chFeatNode)
{
	pColor.readColorChannel(chFeatNode["pColor"]);
	
	pGradMag.readGradientMagnitude(chFeatNode["pGradMag"]);

	pGradHist.readGradientHistogram(chFeatNode["pGradHist"]);	

	complete = chFeatNode["complete"];
}
