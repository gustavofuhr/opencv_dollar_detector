opencv_dollar_detector
======================

&nbsp;&nbsp;&nbsp;&nbsp;This is a work-in-progress port of the Dóllar pedestrian detector to OpenCV in C++ by Charles Arnoud, oriented by Cláudio Rosito Jüng and with the help of Gustavo Führ.  


Current Status  
======================  

&nbsp;&nbsp;&nbsp;&nbsp;Detection working for multiple images, but data set reading is not yet implemented. Still has a lot of prints used for debug.  


To Do List:  
======================  

Current Total: 7  

Opencv_Dollar_Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;enable use of data sets as parameters  
&nbsp;&nbsp;&nbsp;&nbsp;add some way to calculate time elapsed in detection  

Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;nonMaximalSuppression:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly relocate suppression to other class  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add the other two types of suppression    

ColorChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

GradientMagnitudeChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

QuantizedGradientChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;computeMultiScaleChannelFeaturePyramid:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add calculation of lambdas  
&nbsp;&nbsp;&nbsp;&nbsp;computeSingleScaleChannelFeatures:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly add computation of custom channels  

Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;decide what to do with the channelTypes variable    

utils.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  
