#include "Pyramid.h"

Pyramid computeMultiScaleChannelFeaturePyramid(Mat)
{
    ;
}

Pyramid computeSingleScaleChannelFeaturePyramid(Mat I, Channel pChns)
{
    //crop I so it becomes divisible by shrink
    int height = I.rows - (I.rows % pChns.shrink);
    int width =  I.cols - (I.cols % pChns.shrink);


    //compute color channels

}
