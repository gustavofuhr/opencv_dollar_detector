#include "ColorSpace.h"

Mat ColorSpace::rgbConvert(Mat I)
{
    Mat result;
    if (this->colorSpaceType == "luv")
        result = rgb2luv(I);
    else
        if (this->colorSpaceType == "hsv")
           result = rgb2hsv(I);
        else
            if (this->colorSpaceType == "gray")
                result = rgb2gray(I);
    return result;
}

Mat rgb2luv(Mat I)
{
    //image will be represented as h,s,v
    Mat result;
    uchar h, s, v, maxValue=0, minValue=0;
    for (int i = 0; i < I.rows*I.cols*I.channels(); i = i+3)
    {
        uchar b = I.data[i*I.step];
        uchar g = I.data[i*I.step+1];
        uchar r = I.data[i*I.step+2];
        if (r == g && g == b)
        {
            result.data[i*I.step]   = 0; //hue
            result.data[i*I.step+1] = 0; //saturation
            result.data[i*I.step+2] = r; //value
        }
        else
        {
          if (r>=g && g>=b)
          {
              maxValue = r;
              if (g<b)
                  minValue = g;
              else
                  minValue = b;
              h = (g-b)/(maxValue-minValue)+6;
              if(h>=6)
                  h-=6;
          }
          else
            if(g>=r && g>=b)
            {
                maxValue = g;
                if(r<b)
                    minValue = r;
                else
                    minValue = b;
                h = (b-r)/(maxValue-minValue)+2;
            }
            else
            {
                maxValue = b;
                if (r<g)
                    minValue = r;
                else
                    minValue = g;
                h = (r-g)/(maxValue-minValue)+4;
            }
          result.data[i*I.step] = h*(1/6.0);
          result.data[i*I.step+1] = 1-minValue/maxValue;
          result.data[i*I.step+2] = maxValue;
        }
    }
    return result;
}

Mat rgb2hsv(Mat I)
{
    for (int i = 0; i < I.rows*I.cols*I.channels(); i = i+3)
    {
        ;
    }
}

Mat rgb2gray (Mat I)
{
    Mat grayImage;
    double redMultiplier   = 0.2989360213;
    double greenMultiplier = 0.5870430745;
    double blueMultiplier  = 0.1140209043;
    for (int i = 0; i < I.rows*I.cols*I.channels(); i = i+3)
    {
        grayImage.data[i] = I.data[i*I.step]*blueMultiplier +
        I.data[i*I.step+1]*greenMultiplier + I.data[i*I.step]*redMultiplier;
    }
    return grayImage;
}
