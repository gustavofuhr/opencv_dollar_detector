opencv_dollar_detector
======================

This is a work-in-progress port of the Dóllar pedestrian detector to OpenCV in C++





To Do List:
======================

Current Total: 7

Opencv_Dollar_Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;enable full use of data sets as parameters  
&nbsp;&nbsp;&nbsp;&nbsp;enable detection and feature calculation to be used in multiple images  

Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;acfDetect:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;separate detection in a scale to a single function  
&nbsp;&nbsp;&nbsp;&nbsp;bbNms:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly relocate suppression to other class  

ColorChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

GradientMagnitudeChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

QuantizedGradientChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;computeMultiScaleChannelFeaturePyramid:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;change literal number of channels to nChannels whenever needed  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add calculation of lambdas  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;investigate segmentation fault in scale approximation for images of other dimensions from previously tested 
&nbsp;&nbsp;&nbsp;&nbsp;computeSingleScaleChannelFeatures:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly add computation of custom channels  

utils.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  
