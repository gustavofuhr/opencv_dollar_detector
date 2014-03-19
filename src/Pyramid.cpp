#include "Pyramid.h"

Pyramid computeMultiScaleChannelFeaturePyramid(Mat)
{
    ;
}

Pyramid Pyramid::computeSingleScaleChannelFeaturePyramid(Mat I, Channel pChns)
{
    //crop I so it becomes divisible by shrink
    int height = I.rows - (I.rows % pChns.shrink);
    int width =  I.cols - (I.cols % pChns.shrink);


    //compute color channels
    I = this->pChns.pColor.rgbConvert(I);
}

Mat TriangleFilterConvolution(Mat I, int r, int s, int nomex)
{
    Mat result = I;
    if(!I.empty() && !(r==0 && s==1))
    {
        int m = min(I.rows, I.cols);
        if(m>=4 && m>2*r+1)
        {
            if(nomex==0)
            {
                if(r>0 && r<=1 && s<=2)
                    result = convConst("convTri1",I,12/r/(r+2)-2,s);
                else
                    result = convConst("convTri", I, r, s);
            }
            else
            {
                if (r <= 1)
                {
                    double p=12/r/(r+2)-2;
                    f = [1 p 1]/(2+p);
                    r=1;
                }
                else
                {
                    f=[1:r r+1 r:-1:1]/(r+1)^2;
                }
                result = padarray(I,[r r],"symmetric","both");
                //result = convn(convn(result,f,"valid"),f',"valid");
                if (s>1)
                {
                    t=floor(s/2)+1;
                    result = result(t:s:end-s+t,t:s:end-s+t,:);
                }
            }
        }
    }
    return result;
}
