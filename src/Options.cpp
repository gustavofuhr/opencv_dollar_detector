#include "Options.h"

void Options::readOptions(cv::FileNode optionsNode)
{
	pPyramid.readPyramid(optionsNode["pPyramid"]);

	modelDs[0] = optionsNode["modelDs"][0];
	modelDs[1] = optionsNode["modelDs"][1];
	modelDsPad[0] = optionsNode["modelDsPad"][0];
	modelDsPad[1] = optionsNode["modelDsPad"][1];
	stride = optionsNode["stride"];
	cascadeThreshold = optionsNode["cascThr"];

	cv::String temp = (cv::String)optionsNode["pNms"]["type"];
	if (temp == "max")
		suppressionType = MAX;
	else if (temp == "maxg")
		suppressionType = MAXG;
	else if (temp == "ms")
		suppressionType = MS;
	else if (temp == "cover")
		suppressionType = COVER;
	else 
		suppressionType = NONE;

	overlapArea = optionsNode["pNms"]["overlap"];
	overlapDenominator = (cv::String)optionsNode["pNms"]["ovrDnm"];
	suppressionThreshold = 0.0;
	suppressionThreshold = optionsNode["pNms"]["threshold"];
}
