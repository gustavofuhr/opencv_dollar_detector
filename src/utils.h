#ifndef UTILS_H
#define UTILS_H

#include "opencv.hpp"
#include "sse.hpp"
//#include <typeinfo> // for sse rgbConvert functions

#include "Info.h"

#define COLOR_CHANNEL 0
#define PI 3.14159265f

typedef unsigned int uint32;

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX};
enum colorSpaceType {GRAY=0, RGB, LUV, HSV, ORIG};
enum paddingType { REPLICATE=1, SYMMETRIC, CIRCULAR };

float* cvImage2floatArray(cv::Mat source, int channels);
cv::Mat floatArray2cvImage(float* source, int rows, int cols, int channels);
float* cvMat2floatArray(cv::Mat source, int length1, int length2, int length3);
cv::Mat floatArray2cvMat(float* source, int length1, int length2, int length3);
uint32_t* cvMat2charArray(cv::Mat source, int channels);
cv::Mat charArray2cvMat(uint32_t* source, int rows, int cols, int channels);
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);
cv::Mat convolution(cv::Mat source, int channels, int radius, int s, int flag);
cv::Mat resample(cv::Mat source, int ori_h, int ori_w, int new_h, int new_w, float nrm, int channels);
cv::Mat rgbConvert(cv::Mat, int);
cv::Mat padImage(cv::Mat source, int channels, int *pad, int padSize, int type);
void printElements(float* source, cv::String name);
void testFeatures(Info features, cv::String name);

#endif