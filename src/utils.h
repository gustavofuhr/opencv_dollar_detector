#ifndef UTILS_H
#define UTILS_H
#include "opencv.hpp"

float* cvMat2floatArray(cv::Mat source);
cv::Mat floatArray2cvMat(float* source, int rows, int cols, int type);
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);

#endif