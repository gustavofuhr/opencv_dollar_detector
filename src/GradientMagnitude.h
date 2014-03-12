#ifndef GRADMAG_H
#define GRADMAG_H

class GradientMagnitude
{
public:
    int enabled; //if true enable gradient magnitude channel
    int colorChannelIndex; //if>0 color channel to use for grad computation
    int normalizationRadius; //normalization radius for gradient
    double normalizationConstant; //normalization constant for gradient
    int full; //if true compute angles in [0,2*pi) else in [0,pi)
    int nChannels;
    String padWith;
};

#endif
