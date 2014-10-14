opencv_dollar_detector
======================

This is a work-in-progress port of the Dóllar pedestrian detector to OpenCV in C++





To Do List:
======================

Current Total: 14

Opencv_Dollar_Detector.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;enable full use of datasets as parameters
&nbsp;&nbsp;&nbsp;&nbsp;enable detection and feature calculation to be used in multiple images

Detector.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;acfDetect:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;separate detection in a scale to a single function
&nbsp;&nbsp;&nbsp;&nbsp;nmsMax:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 248: sort bounding boxes
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 253: add test for greediness
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 256: test if result[j] was discarded
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 291: discard the bounding box
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 296: adjust result

ColorChannel.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;OK!

GradientMagnitudeChannel.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;OK!
	
QuantizedGradientChannel.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;OK!

Pyramid.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;computeMultiScaleChannelFeaturePyramid:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;change literal number of channels to nChannels whenever needed
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 113: the whole for statement needs work, write a substitution for is=is(2:3) 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 140: compute lambdas
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 189: fix approximated (resampled) scales of H
&nbsp;&nbsp;&nbsp;&nbsp;computeSingleScaleChannelFeatures:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;line 247: possibly add computation of custom channels

pNms.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;check if it is necessary!

utils.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;OK!
