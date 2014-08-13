#include "opencv.hpp"
#include "highgui.hpp"

float* cvMatToFloatArray(cv::Mat source)
{
	float* result;
	float* tempFloat;
	cv::Mat tempMat;
	int resultIndex=0;

	result = (float*)malloc(source.rows*source.cols*3*sizeof(float));

	//first, we need the type conversion
	source.convertTo(tempMat, CV_32FC3, 1.0/255.0);
	tempFloat = (float*)tempMat.data;

	//cvMat is BGR, row then column
	//i need result to be column then row, one channel after the other
	//the next step is changing the way the rows, columns and channels are arranged
	for (int column=0; column < source.cols; column++)
		for (int row=0; row < source.rows; row++)
			for (int channel=0; channel < 3; channel++)
			{
				result[resultIndex++] = tempFloat[column*3+row*source.cols+channel];
			}

	return result;
}

cv::Mat floatArrayToCvMat(float* source, int rows, int cols, int type)
{
	cv::Mat result(rows, cols, type);
	float* tempFloat;
	tempFloat = (float*)malloc(rows*cols*3*sizeof(float));
	int tempIndex=0;


	for (int channel=0; channel < 3; channel++) {
		for (int row=0; row < rows; row++) {
			for (int column=0; column < cols; column++) {
				tempFloat[tempIndex++] = source[channel*rows*cols+column*rows+row];
			}
		}
	}

	result.data = (uchar*)tempFloat;

	return result;
}

void print_array(float *v, int n) {
	for (int i=0; i<n; ++i)
		std::cout << v[i] << " ";
		//printf("%f ", v[i]);
	//printf("\n");
	std::cout << std::endl;
}


int main() {

	cv::Mat image = cv::imread("./simple.png");
	
	cv::namedWindow("before");
	cv::waitKey();
	cv::imshow("before", image);
	cv::waitKey();
	
	cv::Mat fimage;
	image.convertTo(fimage, CV_32FC3, 1/255.0);
	print_array((float*)fimage.data, fimage.cols*fimage.rows*3);

	float* I;
	std::cout << "before cvMatToFloatArray" << std::endl;
	I = cvMatToFloatArray(image);
	std::cout << "after cvMatToFloatArray" << std::endl;

	print_array(I, image.cols*image.rows*3);

	cv::Mat testMat;
	testMat = floatArrayToCvMat(I, image.rows, image.cols, CV_32FC3);

	/*cv::namedWindow("testing conversion");
	cv::waitKey();
	cv::imshow("testing conversion", image);
	cv::waitKey();*/

	return 0;

}