#ifndef UTILS_H
#define UTILS_H

#include <cv.h>
#include <highgui.h> //to use imread
#include "sse.hpp"
#include <iomanip>
#include <typeinfo> // for sse rgbConvert functions
#include <dirent.h>

#include "Info.h"

#define PI 3.14159265f

typedef unsigned int uint32;

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX };
enum colorSpaceType { GRAY=0, RGB, LUV, HSV, ORIG };
enum paddingType { REPLICATE=1, SYMMETRIC, CIRCULAR };
enum suppressionTypeID { NONE, MAX, MAXG, MS, COVER };

// matrix conversions
float* cvImage2floatArray(cv::Mat source, int channels);
cv::Mat floatArray2cvImage(float* source, int rows, int cols, int channels);
void features2floatArray (Info features, float* result, int rows, int cols, int colorChannels, int magChannels, int histChannels);

// Dóllar's allocation functions
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);

// image operations
float* convolution(float* source, int rows, int cols, int channels, int radius, int s, int flag);
float* resample(float* source, int ori_h, int ori_w, int new_h, int new_w, float nrm, int channels);
float* rgbConvert(float *I, int n, int d, int flag, float nrm);
cv::Mat padImage(cv::Mat source, int channels, int *pad, int padSize, int type);

std::vector<std::string> getDataSetFileNames(std::string directory);

#endif