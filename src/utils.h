#ifndef UTILS_H
#define UTILS_H

#include "opencv.hpp"
#include "sse.hpp"
#include <iomanip>
#include <typeinfo> // for sse rgbConvert functions

#include "Info.h"

#define PI 3.14159265f

typedef unsigned int uint32;

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX};
enum colorSpaceType {GRAY=0, RGB, LUV, HSV, ORIG};
enum paddingType { REPLICATE=1, SYMMETRIC, CIRCULAR };

// matrix conversions
float* cvImage2floatArray(cv::Mat source, int channels);
cv::Mat floatArray2cvImage(float* source, int rows, int cols, int channels);

// Dóllar's allocation functions
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);

// image operations
cv::Mat convolution(cv::Mat source, int channels, int radius, int s, int flag);
cv::Mat resample(cv::Mat source, int ori_h, int ori_w, int new_h, int new_w, float nrm, int channels);
cv::Mat rgbConvert(cv::Mat, int);
cv::Mat padImage(cv::Mat source, int channels, int *pad, int padSize, int type);

// debug prints
void printElements(float* source, int rows, cv::String name);
void testFeatures(Info features, cv::String name);

#endif