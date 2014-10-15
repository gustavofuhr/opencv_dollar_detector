cmake_minimum_required (VERSION 2.6)
SET(CMAKE_CXX_COMPILER g++)
project (opencv_dollar_detector)
find_package( OpenCV REQUIRED )
add_executable(opencv_dollar_detector 
./src/Opencv_Dollar_Detector.cpp 
./src/Detector.cpp
./src/ColorChannel.cpp
./src/GradientMagnitudeChannel.cpp
./src/Pyramid.cpp
./src/QuantizedGradientChannel.cpp
./src/ChannelFeatures.cpp
./src/Options.cpp
./src/Clf.cpp
./src/BoundingBox.cpp
./src/utils.cpp)
include_directories("/usr/include/opencv2")
include_directories("/usr/include/opencv2/highgui")
target_link_libraries( opencv_dollar_detector ${OpenCV_LIBS} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3")
