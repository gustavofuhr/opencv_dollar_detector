#pragma once

#include "opencv/cv.h"
#include "opencv/highgui.h"
#include "Detector.h"
#include <QApplication>
#include "mainwindow.h"

#define PREFIX "D:\\Workspace\\CppAcfDetect\\CppAcfDetect\\"
#define EXTENSION ".jpg"

using namespace cv;
using namespace std;

BoundingBox* acfDetectImg(Mat,Detector);
