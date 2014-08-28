#include "ColorChannel.h"

void ColorChannel::readColorChannel(cv::FileNode colorNode)
{
	enabled = colorNode["enabled"];
	smoothingRadius = colorNode["smooth"];
	colorSpaceType = (cv::String)colorNode["colorSpace"];
	nChannels = colorNode["nChns"];
	padWith = (cv::String)colorNode["padWith"];
}

/******************************************************/
//colorspace conversions
cv::Mat ColorChannel::rgbConvert(cv::Mat I)
{
    cv::Mat result;
    cv::Mat uChar;

    if (this->colorSpaceType == "luv")
    {
        I.convertTo(uChar, CV_8UC3, 255.0);
        cvtColor(uChar, result, CV_BGR2Luv);
        result.convertTo(result, CV_32FC3, 1.0/255.0);
    }
    else
        if (this->colorSpaceType == "hsv")
           cvtColor(I, result, CV_BGR2HSV);
        else
            if (this->colorSpaceType == "gray")
                cvtColor(I, result, CV_BGR2GRAY);
            else
                return I;

    return result;
}