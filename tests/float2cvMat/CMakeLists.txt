cmake_minimum_required (VERSION 2.6)
SET(CMAKE_CXX_COMPILER g++)
project (mat_conversion)
find_package( OpenCV REQUIRED )
add_executable(mat_conversion 
				../../src/ColorChannel.h
				main.cpp) 

#include_directories("/usr/include/opencv2")
#include_directories("/usr/include/opencv2/highgui")

target_link_libraries( mat_conversion ${OpenCV_LIBS} )
