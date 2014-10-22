opencv_dollar_detector
======================

This is a work-in-progress port of the Dóllar pedestrian detector to OpenCV in C++ by Charles Arnoud, under the mentorship of Cláudio Rosito Jüng and with the help of Gustavo Führ.  


Current Status  
======================  

Detection working for multiple images, but process is being killed for lack of memory when applied to a big data set. 


To Do List:  
======================  

Current Total: 10  

Opencv_Dollar_Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;maybe change the way time is being calculated  

Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;nonMaximalSuppression:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly relocate suppression to other class  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add the other two types of suppression    

ColorChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

GradientMagnitudeChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;check if variable d is correct inside mGradMag  
&nbsp;&nbsp;&nbsp;&nbsp;check if mGradMag is not modifying the input  

QuantizedGradientChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  

Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;use matrices I, M, O, S, H as float* inside this class, transform them afterwards  
&nbsp;&nbsp;&nbsp;&nbsp;computeMultiScaleChannelFeaturePyramid:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add calculation of lambdas  
&nbsp;&nbsp;&nbsp;&nbsp;computeSingleScaleChannelFeatures:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly add computation of custom channels  

Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;decide what to do with the channelTypes variable   
&nbsp;&nbsp;&nbsp;&nbsp;fix crop inside chnsCompute 

utils.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;OK!  
