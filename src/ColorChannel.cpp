#include "ColorChannel.h"

void ColorChannel::readColorChannel(cv::FileNode colorNode)
{
	enabled = colorNode["enabled"];
	smooth = colorNode["smooth"];
	colorSpaceType = (cv::String)colorNode["colorSpace"];
	nChannels = colorNode["nChns"];
	padWith = (cv::String)colorNode["padWith"];
}

float* cvMatToFloatArray(cv::Mat img)
{
    assert(img.type() == CV_8UC3);
    int h = img.rows; //height
    int w = img.cols; //width
    int d = 3; //nChannels
    float multiplier = 1 / 255.0f; //rescale pixels to range [0 to 1]

    uchar* img_data = &img.data[0];
    float* I = (float*)malloc(h * w * d * sizeof(float)); //img transposed to Matlab data layout. 
    for(int y=0; y<h; y++){
        for(int x=0; x<w; x++){
            for(int ch=0; ch<d; ch++){
                I[x*h + y + ch*w*h] = img_data[x*d + y*w*d + ch] * multiplier; 
            }
        }
    }
    return I;
}

cv::Mat floatArrayToCvMat(float* I, int h, int w)
{
	cv::Mat img(h, w, CV_32F);
	int indexForI = 0;

	for (int i=0; i < h; i++)
		for (int j=0; j < w; j++)
			img.data[i,j] = I[indexForI++];

	return img;
}

//convolutions taken from the convConst.cpp file
//I'm going to add just the ones that are actually called

//Probably need to change all of this or refactor some things

//this is the wrapper function that call the appropriate one
//this one could go away, if we refactor the other ones to
//operate on cvMat structures rather than just float*
cv::Mat ColorChannel::convolution(cv::Mat source, int radius, int s, int flag)
{
	int indexForI = 0;
	float* O;
	cv::Mat result;
<<<<<<< Updated upstream
	
=======
	int h = source.rows;
	int w = source.cols;
	int d = source.dims;

	cv::imshow("source_1",source);
	cv::waitKey();
	//cv::destroyAllWindows();

	cv::imshow("floatMat_1",floatMat);
	cv::waitKey();
	//cv::destroyAllWindows();

	// float* I = (float*)malloc(source.rows * source.step * sizeof(float));
	// for (int k=0; k < source.rows; k++)
	// 	for (int l=0; l < source.cols; l++)
	// 		I[k*source.step+l] = floatMat.at<float>(k,l);

	cv::Mat tempMat(source.rows, source.cols, CV_32F);
	tempMat.data = (uchar*)floatMat.data;
	
	cv::imshow("tempMat",tempMat);
	cv::waitKey();
	//cv::destroyAllWindows();
>>>>>>> Stashed changes
	
	switch(flag)
	{
		case CONV_TRI: 	
<<<<<<< Updated upstream
					triangleFilterConvolution(source, O, radius, s);
=======
					triangleFilterConvolution((float*)floatMat.data, O, h, w, d, radius, s);
>>>>>>> Stashed changes
					cv::imshow("conv2",source);
					cv::waitKey();				
					//cv::destroyAllWindows();	
					break;
		case CONV_TRI1: 
					int p = 12/radius/(radius+2)-2;
<<<<<<< Updated upstream
					//convTri1(I, O, h, w, d, p, s);
=======
					convTri1((float*)floatMat.data, O, h, w, d, p, s);
>>>>>>> Stashed changes
					break;
	}
	floatMat.data = (uchar*)O;

	floatMat.convertTo(result, CV_8UC3);
	cv::imshow("result", result);
	cv::waitKey();

	//result.data = (uchar*)O;
	return result;
}

// convolve I by a 2rx1 triangle filter (uses SSE)
//aka convTri
void ColorChannel::triangleFilterConvolution(cv::Mat source, float *O, int r, int s)
{
	int h = source.rows;
	int w = source.cols;
	int d = source.dims;
	r++; float nrm = 1.0f/(r*r*r*r); int i, j, k=(s-1)/2, h0, h1, w0;
	if(h%4==0) 
		h0=h1=h; 
	else 
	{ 
		h0=h-(h%4); 
		h1=h0+4; 
	} 
	w0=(w/s)*s;
	
	float *T=(float*) malloc(2*h1*sizeof(float)), *U=T+h1;

<<<<<<< Updated upstream
	//start of the debug section
	cv::Mat floatMat(h, w, CV_32F);
	source.convertTo(floatMat, CV_32F, 1.0/255.0);
	cv::imshow("testing Mat after convertTo applied to source",floatMat);
	cv::waitKey();
	cv::destroyAllWindows();

	float* I = cvMatToFloatArray(source);
	cv::Mat test1(source.rows, source.cols, CV_32F, I);
	cv::imshow("testing Mat after assignment from I",test1);

	float* I2 = (float*) floatMat.data;
	cv::Mat test2(source.rows, source.cols, CV_32F, I2);
	cv::imshow("testing Mat after assignment from I2 without step",test2);

	float* I3 = (float*) floatMat.data;
	cv::Mat test3(source.rows, source.cols, CV_32F, I3, floatMat.step);
	cv::imshow("testing Mat after assignment from I3 with step",test3);

=======
	cv::Mat input(h, w, CV_32F);
	//input.convertTo(input, CV_32F);
	input.data = (uchar*)I;
	cv::imshow("triangle 0",input);
>>>>>>> Stashed changes
	cv::waitKey();
	cv::destroyAllWindows();


	// in here, we go back to normal processing
	while(d-- > 0) 
	{
		// initialize T and U
		for(j=0; j<h0; j+=4) 
			STR(U[j], STR(T[j], LDu(I[j])));
		for(i=1; i<r; i++) 
			for(j=0; j<h0; j+=4) 
				INC(U[j],INC(T[j],LDu(I[j+i*h])));
		for(j=0; j<h0; j+=4) 
			STR(U[j],MUL(nrm,(SUB(MUL(2,LD(U[j])),LD(T[j])))));
		for(j=0; j<h0; j+=4) 
			STR(T[j],0);
		for(j=h0; j<h; j++ ) 
			U[j]=T[j]=I[j];
		for(i=1; i<r; i++) 
			for(j=h0; j<h; j++ ) 
				U[j]+=T[j]+=I[j+i*h];
		for(j=h0; j<h; j++ ) 
		{ 
			U[j] = nrm * (2*U[j]-T[j]); 
			T[j]=0; 
		}

		// prepare and convolve each column in turn
		for( i=0; i<w0; i++ ) 
		{
			float *Il, *Ir, *Im; 
			Il=Ir=Im=I; 
			Im+=(i-1)*h;
			if( i<=r ) 
			{ 
				Il+=(r-i)*h; 
				Ir+=(r-1+i)*h; 
			}
			else 
				if( i<=w-r ) 
				{ 
					Il-=(r+1-i)*h; 
					Ir+=(r-1+i)*h; 
				}
				else 
				{ 
					Il-=(r+1-i)*h; 
					Ir+=(2*w-r-i)*h; 
				}
			if(i) 
				for( j=0; j<h0; j+=4 ) 
				{
					__m128 del = SUB(ADD(LDu(Il[j]),LDu(Ir[j])),MUL(2,LDu(Im[j])));
					INC(U[j], MUL(nrm,(INC(T[j],del))));
				}
			if(i) 
				for( j=h0; j<h; j++ ) 
					U[j]+=nrm*(T[j]+=Il[j]+Ir[j]-2*Im[j]);
			cv::imshow("triangle 1",floatMat);
			cv::waitKey();
			cv::destroyAllWindows();
			k++; 
			if(k==s) 
			{
				k=0; 
				//the segmentation fault is caused in convTri1Y
				convTri1Y(U,O,h,r-1,s); 
				O+=h/s; 
			}
			cv::imshow("triangle 2",floatMat);
			cv::waitKey();
			cv::destroyAllWindows();
		}
		I+=w*h;
	}
	free(T);
}

// convolve one column of I by [1 p 1] filter (uses SSE)
void ColorChannel::convTri1Y( float *I, float *O, int h, float p, int s ) {
	cv::Mat testMat = cv::Mat::zeros( 300, 300, CV_8UC3 );
	cv::imshow("inside convTri1Y",testMat);
	cv::waitKey();
	cv::destroyAllWindows();

	#define C4(m,o) ADD(ADD(LDu(I[m*j-1+o]),MUL(p,LDu(I[m*j+o]))),LDu(I[m*j+1+o]))
  int j=0, k=((~((size_t) O) + 1) & 15)/4, h2=(h-1)/2;
  if( s==2 ) 
  {
  	cv::imshow("inside convTri1Y's if",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
    for( ; j<k; j++ ) 
    	O[j]=I[2*j]+p*I[2*j+1]+I[2*j+2];
    for( ; j<h2-4; j+=4 ) 
    	STR(O[j],_mm_shuffle_ps(C4(2,1),C4(2,5),136));
    for( ; j<h2; j++ ) 
    	O[j]=I[2*j]+p*I[2*j+1]+I[2*j+2];
    if( h%2==0 ) 
    	O[j]=I[2*j]+(1+p)*I[2*j+1];
  } 
  else 
  {
  	cv::imshow("inside convTri1Y's else 1",testMat);
		cv::waitKey();
		cv::destroyAllWindows();

	//this operation causes segmentation fault
    O[j]=(1+p)*I[j]+I[j+1]; 

    cv::imshow("inside convTri1Y's else 1.1",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
    j++; 
    cv::imshow("inside convTri1Y's else 2",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
    if(k==0) 
    	k=(h<=4) ? h-1 : 4;
    cv::imshow("inside convTri1Y's else 2.1",testMat);
		cv::waitKey();
		cv::destroyAllWindows();

	//segmenteation fault in here too
    for( ; j<k; j++ ) 
    	O[j]=I[j-1]+p*I[j]+I[j+1];
    cv::imshow("inside convTri1Y's else 3",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
    for( ; j<h-4; j+=4 ) 
    	STR(O[j],C4(1,0));
    cv::imshow("inside convTri1Y's else 4",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
    for( ; j<h-1; j++ ) 
    	O[j]=I[j-1]+p*I[j]+I[j+1];
    cv::imshow("inside convTri1Y's else 5",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
    O[j]=I[j-1]+(1+p)*I[j];
    cv::imshow("inside convTri1Y's else 5.1",testMat);
		cv::waitKey();
		cv::destroyAllWindows();
  }
  #undef C4
}

// convTri1( A, B, ns[0], ns[1], d, p, s );
// convConst('convTri1',I,12/r/(r+2)-2,s);
// convolve I by [1 p 1] filter (uses SSE)
void ColorChannel::convTri1( float *I, float *O, int h, int w, int d, float p, int s ) {
  const float nrm = 1.0f/((p+2)*(p+2)); int i, j, h0=h-(h%4);
  float *Il, *Im, *Ir, *T=(float*) malloc(h*sizeof(float));
  for( int d0=0; d0<d; d0++ ) for( i=s/2; i<w; i+=s ) 
  {
    Il=Im=Ir=I+i*h+d0*h*w; if(i>0) Il-=h; if(i<w-1) Ir+=h;
    for( j=0; j<h0; j+=4 )
      STR(T[j],MUL(nrm,ADD(ADD(LDu(Il[j]),MUL(p,LDu(Im[j]))),LDu(Ir[j]))));
    for( j=h0; j<h; j++ ) T[j]=nrm*(Il[j]+p*Im[j]+Ir[j]);
    convTri1Y(T,O,h,p,s); O+=h/s;
  }
  free(T);
}



/******************************************************/
//colorspace conversions
cv::Mat ColorChannel::rgbConvert(cv::Mat I)
{
    cv::Mat result;
    if (this->colorSpaceType == "luv")
				cvtColor(I, result, CV_BGR2Luv);
    else
        if (this->colorSpaceType == "hsv")
           cvtColor(I, result, CV_BGR2HSV);
        else
            if (this->colorSpaceType == "gray")
                cvtColor(I, result, CV_BGR2GRAY);
    return result;
}
