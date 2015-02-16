#ifndef UTILS_H
#define UTILS_H

#include <cv.h>
#include <highgui.h> //to use imread
#include "sse.hpp"
#include <iomanip>
#include <typeinfo> // for sse rgbConvert functions
#include <dirent.h>
#include "BoundingBox.h"

#include "Info.h"

#define PI 3.14159265f

typedef unsigned int uint32;

enum method { CONV_BOX, CONV_TRI, CONV_11, CONV_TRI1, CONV_MAX };
enum colorSpaceType { GRAY=0, RGB, LUV, HSV, ORIG };
enum paddingType { REPLICATE=1, SYMMETRIC, CIRCULAR };
enum suppressionTypeID { NONE, MAX, MAXG, MS, COVER };

// matrix conversions
void cvMat2floatArray(cv::Mat source, float* result, int channels);
cv::Mat floatArray2cvMat(float* source, int rows, int cols, int channels);
//void floatArray2cvData(float* source, float* result, int rows, int cols, int channels);
void features2floatArray (Info features, float* result, int rows, int cols, int colorChannels, int magChannels, int histChannels);

// Dóllar's allocation functions
void* alMalloc( size_t size, int alignment );
void alFree(void* aligned);

// image operations
void convolution(float* source, float* result, int rows, int cols, int channels, int radius, int s);
void resample(float *A, float *B, int ha, int hb, int wa, int wb, int d, float r );
//float* rgbConvert(float *I, int n, int d, int flag, float nrm);
void rgbConvert(float *I, float* J, int n, int d, int flag, float nrm);
cv::Mat padImage(cv::Mat source, int channels, int *pad, int padSize, int type);

std::vector<std::string> getDataSetFileNames(std::string directory);

cv::Point2f imagePoint2groundPlanePoint(float imageU, float imageV, float imageZ, cv::Mat homography);
cv::Point worldPoint2imagePoint(float worldX, float worldY, float worldZ, cv::Mat homography);
std::vector<cv::Point2f> findGroundPlaneAndImageIntersectionPoints(int imageWidth, int imageHeight, int boundingBoxImageWidth, int boundingBoxImageHeight, 
  float maxPedestrianHeight, cv::Mat_<float> &projection, cv::Mat_<float> &homography);
double findWorldHeight(int u, int bottom_v, int top_v, cv::Mat_<float> &projection, cv::Mat_<float> &homography);
std::vector<cv::Point2f> trimGroundPlanesBottomPoints(int imageWidth, int imageHeight, int boundingBoxImageWidth, int boundingBoxImageHeight, 
  double targetPedestrianWorldHeight, std::vector<cv::Point2f> groundPlaneLimits, cv::Mat_<float> &projection, cv::Mat_<float> &homography);
int findNecessaryNumberOfOctaves(int imageWidth, int imageHeight, int boundingBoxImageWidth, int boundingBoxImageHeight, 
  double minPedestrianWorldHeight, double maxPedestrianWorldHeight, cv::Mat_<float> &projection, cv::Mat_<float> &homography);
double findLastNecessaryScaleInAPoint(int u, int v, int imageHeight, int boundingBoxImageHeight, double maxPedestrianWorldHeight, 
  cv::Mat_<float> &projection, cv::Mat_<float> &homography);

cv::Mat world2image(cv::Mat &w_point, cv::Mat_<float> &P);
BoundingBox wcoord2bbox(cv::Point2f w_point, cv::Mat_<float> &P, float w_height, float aspect_ratio);

void print_fmatrix(const std::string &title, const cv::Mat &m);
void print_dmatrix(const std::string &title, const cv::Mat &m);

double log_base_n(double y, double base);


#endif