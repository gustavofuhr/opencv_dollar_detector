#include <cv.h>
#include <highgui.h>


int main () {

	cv::Mat image1, image_temp;

	// um exemplo dessa conversão de cv::Mat para float* 
	//e de volta para cv::Mat que não altere a imagem original?
	image1 = cv::imread("./frame0254.png");//,0);

	cv::imshow("Original image cv::Mat", image1);
	cv::waitKey(0);

	// convertendo a imagem para float
	// nota que eu estou escalando todos os pixels para ficarem
	// entre 0 e 1, pois assim a visualizacao funciona legal.
	image1.convertTo(image_temp, CV_32FC1, 1.0/255.0);
	cv::imshow("Converted to float (still cv::Mat)", image_temp);
	cv::waitKey(0);

	float *float_mat;
	float_mat = (float*)image_temp.data;
	// o ponteiro agora referencia os dados de image_temp.

	// eu posso criar uma nova imagem do mesmo tamanho e apontar
	// para os dados apontados pelo ponteiro de float.
	cv::Mat image2(image_temp.rows, image_temp.cols, image_temp.type()); // nao eh preciso alocar espaco aqui
	image2.data = (uchar*)float_mat;

	// exibimos a image2 para ver se esta tudo certo
	cv::imshow("cv::Mat created from float pointer", image2);
	cv::waitKey(0);
	

	return 0;
}