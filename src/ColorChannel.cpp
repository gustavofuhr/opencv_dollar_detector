#include "ColorChannel.h"

void ColorChannel::readColorChannel(cv::FileNode colorNode)
{
	enabled = colorNode["enabled"];
	smoothingRadius = colorNode["smooth"];
	nChannels = colorNode["nChns"];
	padWith = (cv::String)colorNode["padWith"];
    if ((cv::String)colorNode["colorSpace"] == "gray")
        colorSpaceType = GRAY;
    if ((cv::String)colorNode["colorSpace"] == "rgb")
        colorSpaceType = RGB;
    if ((cv::String)colorNode["colorSpace"] == "luv")
        colorSpaceType = LUV;
    if ((cv::String)colorNode["colorSpace"] == "hsv")
        colorSpaceType = HSV;
    if ((cv::String)colorNode["colorSpace"] == "orig")
        colorSpaceType = ORIG;
}

/******************************************************/
