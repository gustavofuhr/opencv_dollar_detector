#include "ColorChannel.h"

void ColorChannel::readColorChannel(cv::FileNode colorNode)
{
	enabled = colorNode["enabled"];
	smooth = colorNode["smooth"];
	colorSpaceType = (cv::String)colorNode["colorSpace"];
	nChannels = colorNode["nChns"];
	padWith = (cv::String)colorNode["padWith"];
}

//convolutions taken from the convConst.cpp file
//I'm going to add just the ones that are actually called

//Probably need to change all of this or refactor some things

//this is the wrapper function that call the appropriate one
//this one could go away, if we refactor the other ones to
//operate on cvMat structures rather than just float*
cv::Mat ColorChannel::convolution(cv::Mat source, int radius, int s, int flag)
{
	//J = convConst('convTri',I,r,s);
	//nDims = mxGetNumberOfDimensions(prhs[1]);
	//ns = (int*) mxGetDimensions(prhs[1]);
	//d = (nDims == 3) ? ns[2] : 1;
	//s = (int) mxGetScalar(prhs[3]);
	//ms[0]=ns[0]/s; ms[1]=ns[1]/s; ms[2]=d;
	//B = (float*) mxMalloc(ms[0]*ms[1]*d*sizeof(float));
	float* O = (float*)malloc(source.rows/s*source.cols/s*3*sizeof(float));

	cv::Mat floatMat;
	source.convertTo(floatMat, CV_32FC1, 1.0/255.0);

	float* I;
	I = (float*)floatMat.data;

	switch(flag)
	{
		case CONV_TRI: 	
					triangleFilterConvolution(I, O, source.rows, source.cols, source.dims, radius, s);
					break;
		case CONV_TRI1: 
					int p = 12/radius/(radius+2)-2;
					convTri1(I, O, source.rows, source.cols, source.dims, p, s);
					break;
	}

	cv::Mat result(floatMat.rows, floatMat.cols, floatMat.type());
	result.data = (uchar*)O;

	//prints for debug
	std::cout << "result: rows = "<< result.rows << ", cols = " << result.cols << ", dims = " << result.dims << std::endl;
	std::cout << "value 0,0 = " << result.at<double>(0,0) << std::endl;
	std::cout << "value 0,100 = " << result.at<double>(0,100) << std::endl;
	std::cout << "value 0,719 = " << result.at<double>(0,719) << std::endl;
	std::cout << "value 1,0 = " << result.at<double>(1,0) << std::endl;
	std::cout << "value 100,0 = " << result.at<double>(100,0) << std::endl;
	std::cout << "value 10,10 = " << result.at<double>(10,10) << std::endl;
	std::cout << "value 100,100 = " << result.at<double>(100,100) << std::endl;
	std::cout << "value 300,300 = " << result.at<double>(300,300) << std::endl;
	std::cout << "value 570,710 = " << result.at<double>(570,710) << std::endl;
	std::cout << "inside convolution, before printing result" << std::endl;
	cv::imshow("inside convolution, printing floatMat (input)", floatMat);
	cv::imshow("inside convolution, image after conversion from float *O", result);
	cv::waitKey();				
	//end of debug section

	return result;
}

// convolve I by a 2rx1 triangle filter (uses SSE)
//aka convTri
void ColorChannel::triangleFilterConvolution( float *I, float *O, int h, int w, int d, int r, int s )
{
	r++; 
	float nrm = 1.0f/(r*r*r*r); 
	int i, j, k=(s-1)/2, h0, h1, w0;
	
	if(h%4==0) 
		h0=h1=h; 
	else 
	{ 
		h0=h-(h%4); 
		h1=h0+4; 
	} 
	w0=(w/s)*s;
	
	//the function uses alMalloc, dont know if it makes a difference
	//float *T=(float*) alMalloc(2*h1*sizeof(float),16), *U=T+h1;
	float *T=(float*) calloc(2*h1*sizeof(float), 16), *U=T+h1;

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
			float *Il, *Ir, *Im; Il=Ir=Im=I; 
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
			k++; 
			if(k==s) 
			{ 
				k=0; 
				convTriY(U,O,h,r-1,s); 
				O+=h/s; 
			}
		}
		I+=w*h;
	}
	//the function uses alFree, dont know if it makes a difference
	free(T);
}

// convolve one column of I by a 2rx1 triangle filter
void ColorChannel::convTriY( float *I, float *O, int h, int r, int s ) {
  r++; float t, u; int j, r0=r-1, r1=r+1, r2=2*h-r, h0=r+1, h1=h-r+1, h2=h;
  u=t=I[0]; for( j=1; j<r; j++ ) u+=t+=I[j]; u=2*u-t; t=0;
  if( s==1 ) {
    O[0]=u; j=1;
    for(; j<h0; j++) O[j] = u += t += I[r-j]  + I[r0+j] - 2*I[j-1];
    for(; j<h1; j++) O[j] = u += t += I[j-r1] + I[r0+j] - 2*I[j-1];
    for(; j<h2; j++) O[j] = u += t += I[j-r1] + I[r2-j] - 2*I[j-1];
  } else {
    int k=(s-1)/2; h2=(h/s)*s; if(h0>h2) h0=h2; if(h1>h2) h1=h2;
    if(++k==s) { k=0; *O++=u; } j=1;
    for(;j<h0;j++) { u+=t+=I[r-j] +I[r0+j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
    for(;j<h1;j++) { u+=t+=I[j-r1]+I[r0+j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
    for(;j<h2;j++) { u+=t+=I[j-r1]+I[r2-j]-2*I[j-1]; if(++k==s){ k=0; *O++=u; }}
  }
}

// convolve one column of I by [1 p 1] filter (uses SSE)
void ColorChannel::convTri1Y( float *I, float *O, int h, float p, int s ) {
	#define C4(m,o) ADD(ADD(LDu(I[m*j-1+o]),MUL(p,LDu(I[m*j+o]))),LDu(I[m*j+1+o]))
  int j=0, k=((~((size_t) O) + 1) & 15)/4, h2=(h-1)/2;
  if( s==2 ) 
  {
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
    O[j]=(1+p)*I[j]+I[j+1]; 
    j++; 
    if(k==0) 
    	k=(h<=4) ? h-1 : 4;
    for( ; j<k; j++ ) 
    	O[j]=I[j-1]+p*I[j]+I[j+1];
    for( ; j<h-4; j+=4 ) 
    	STR(O[j],C4(1,0));
    for( ; j<h-1; j++ ) 
    	O[j]=I[j-1]+p*I[j]+I[j+1];
    O[j]=I[j-1]+(1+p)*I[j];
  }
  #undef C4
}

// convTri1( A, B, ns[0], ns[1], d, p, s );
// convConst('convTri1',I,12/r/(r+2)-2,s);
// convolve I by [1 p 1] filter (uses SSE)
void ColorChannel::convTri1( float *I, float *O, int h, int w, int d, float p, int s )
{
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
