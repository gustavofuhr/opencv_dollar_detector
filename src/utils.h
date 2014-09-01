#ifndef UTILS_H
#define UTILS_H

#include "opencv.hpp"
#include "sse.hpp"
#include <stdio.h>
#include <typeinfo>

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX};

float* cvImage2floatArray(cv::Mat source, int channels);
cv::Mat floatArray2cvImage(float* source, int rows, int cols, int channels);
float* cvMat2floatArray(cv::Mat source, int length1, int length2, int length3);
cv::Mat floatArray2cvMat(float* source, int length1, int length2, int length3);
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);
cv::Mat convolution(cv::Mat source, int channels, int radius, int s, int flag);
cv::Mat resample(cv::Mat source, int ori_h, int ori_w, int new_h, int new_w, float nrm, int channels);

#endif