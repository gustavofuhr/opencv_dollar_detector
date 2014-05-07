#include "PNms.h"

void PNms::readPNms(cv::FileNode pNmsNode)
{
	type = (cv::String)pNmsNode["type"];
	overlap = pNmsNode["overlap"];
	ovrDnm = (cv::String)pNmsNode["ovrDnm"];
}
