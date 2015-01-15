opencv_dollar_detector
======================

This is an attemp to improve the Dóllar pedestrian detector's performance using camera calibration by Charles Arnoud, under the mentorship of Cláudio Rosito Jüng and with the help of Gustavo Führ.  


Current Status  
======================  

Detection working, but has two memory leaks and is too slow. Added calibrated detector.


To Do List:  
======================  

Current Total: 13  

Opencv_Dollar_Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;find out why the program is so slow!  
&nbsp;&nbsp;&nbsp;&nbsp;maybe change the way time is being calculated  
&nbsp;&nbsp;&nbsp;&nbsp;save detections to .txt file  
&nbsp;&nbsp;&nbsp;&nbsp;add more comments  

Detector.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;nonMaximalSuppression:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly relocate suppression to other class  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add the other two types of suppression    
  
GradientMagnitudeChannel.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;check if variable d is correct inside mGradMag  
&nbsp;&nbsp;&nbsp;&nbsp;check if mGradMag is not modifying the input  
  
Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;computeMultiScaleChannelFeaturePyramid:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;add calculation of lambdas  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;find sources of two memory leaks  
&nbsp;&nbsp;&nbsp;&nbsp;computeSingleScaleChannelFeatures:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;possibly add computation of custom channels  

Pyramid.cpp:  
&nbsp;&nbsp;&nbsp;&nbsp;decide what to do with the channelTypes variable   

utils.cpp:
&nbsp;&nbsp;&nbsp;&nbsp;find out why convTri function is so much slower than in the original mex file