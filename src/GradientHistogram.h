#ifndef GRADHIST_H
#define GRADHIST_H

#include "opencv/cv.h"

using namespace cv;

class GradientHistogram
{
public:
    int enable; //if true enable gradient histogram channels
    int* binSize; //spatial bin size (defaults to shrink)
    int orientationChannels; //number of orientation channels
    int useSoftBinning; //if true use "soft" bilinear spatial binning
    int useHogNormalization; //if true perform 4-way hog normalization/clipping
    double clipHog; //value at which to clip hog histogram bins
    int nChannels;
    String padWith;
};

#endif
