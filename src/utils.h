#ifndef UTILS_H
#define UTILS_H

#include "opencv.hpp"
#include "sse.hpp"
#include <stdio.h>
#include <typeinfo>

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX};

float* cvMat2floatArray(cv::Mat source);
cv::Mat floatArray2cvMat(float* source, int rows, int cols, int type);
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);
cv::Mat convolution(cv::Mat source, int radius, int s, int flag);
cv::Mat resample(cv::Mat source, int h, int w, float nrm);

#endif