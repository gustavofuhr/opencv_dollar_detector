#include "ColorChannel.h"

//convolutions taken from the convConst.cpp file
//I'm going to add just the ones that are actually called

//Probably need to change all of this or refactor some things

//this is the wrapper function that call the appropriate one
//this one could go away, if we refactor the other ones to
//operate on cvMat structures rather than just float*
Mat ColorChannel::convolution(Mat source, int radius, int s, int flag)
{
	Mat floatMat(source.rows, source.cols, CV_32F);
	//src.convertTo(dst, CV_32F);
	source.convertTo(floatMat, CV_32F);
	float* I;
	int indexForI = 0;
	float* O;
	Mat result;
	int h = source.rows;
	int w = source.cols;
	int d = source.dims;

	/*for (int i=0; i < floatMat.rows; i++)
		for (int j=0; j < floatMat.cols; j++)
				I[indexForI++] = floatMat.at<float>(i,j);*/
	
	switch(flag)
	{
		case CONV_TRI: 	
					triangleFilterConvolution(I, O, h, w, d, radius, s);
					break;
		case CONV_TRI1: 
					int p = 12/radius/(radius+2)-2;
					convTri1(I, O, h, w, d, p, s);
					break;
	}
	result.data = (uchar*)O;
	return result;
}

// convolve I by a 2rx1 triangle filter (uses SSE)
//aka convTri
void ColorChannel::triangleFilterConvolution( float *I, float *O, int h, int w, int d, int r, int s ) {
  r++; float nrm = 1.0f/(r*r*r*r); int i, j, k=(s-1)/2, h0, h1, w0;
  if(h%4==0) h0=h1=h; else { h0=h-(h%4); h1=h0+4; } w0=(w/s)*s;
  float *T=(float*) malloc(2*h1*sizeof(float)), *U=T+h1;
  while(d-- > 0) {
    // initialize T and U
    for(j=0; j<h0; j+=4) STR(U[j], STR(T[j], LDu(I[j])));
    for(i=1; i<r; i++) for(j=0; j<h0; j+=4) INC(U[j],INC(T[j],LDu(I[j+i*h])));
    for(j=0; j<h0; j+=4) STR(U[j],MUL(nrm,(SUB(MUL(2,LD(U[j])),LD(T[j])))));
    for(j=0; j<h0; j+=4) STR(T[j],0);
    for(j=h0; j<h; j++ ) U[j]=T[j]=I[j];
    for(i=1; i<r; i++) for(j=h0; j<h; j++ ) U[j]+=T[j]+=I[j+i*h];
    for(j=h0; j<h; j++ ) { U[j] = nrm * (2*U[j]-T[j]); T[j]=0; }
    // prepare and convolve each column in turn
    for( i=0; i<w0; i++ ) {
      float *Il, *Ir, *Im; Il=Ir=Im=I; Im+=(i-1)*h;
      if( i<=r ) { Il+=(r-i)*h; Ir+=(r-1+i)*h; }
      else if( i<=w-r ) { Il-=(r+1-i)*h; Ir+=(r-1+i)*h; }
      else { Il-=(r+1-i)*h; Ir+=(2*w-r-i)*h; }
      if(i) for( j=0; j<h0; j+=4 ) {
        __m128 del = SUB(ADD(LDu(Il[j]),LDu(Ir[j])),MUL(2,LDu(Im[j])));
        INC(U[j], MUL(nrm,(INC(T[j],del))));
      }
      if(i) for( j=h0; j<h; j++ ) U[j]+=nrm*(T[j]+=Il[j]+Ir[j]-2*Im[j]);
      k++; if(k==s) { k=0; convTri1Y(U,O,h,r-1,s); O+=h/s; }
    }
    I+=w*h;
  }
  free(T);
}

// convolve one column of I by [1 p 1] filter (uses SSE)
void ColorChannel::convTri1Y( float *I, float *O, int h, float p, int s ) {
  #define C4(m,o) ADD(ADD(LDu(I[m*j-1+o]),MUL(p,LDu(I[m*j+o]))),LDu(I[m*j+1+o]))
  int j=0, k=((~((size_t) O) + 1) & 15)/4, h2=(h-1)/2;
  if( s==2 ) {
    for( ; j<k; j++ ) O[j]=I[2*j]+p*I[2*j+1]+I[2*j+2];
    for( ; j<h2-4; j+=4 ) STR(O[j],_mm_shuffle_ps(C4(2,1),C4(2,5),136));
    for( ; j<h2; j++ ) O[j]=I[2*j]+p*I[2*j+1]+I[2*j+2];
    if( h%2==0 ) O[j]=I[2*j]+(1+p)*I[2*j+1];
  } else {
    O[j]=(1+p)*I[j]+I[j+1]; j++; if(k==0) k=(h<=4) ? h-1 : 4;
    for( ; j<k; j++ ) O[j]=I[j-1]+p*I[j]+I[j+1];
    for( ; j<h-4; j+=4 ) STR(O[j],C4(1,0));
    for( ; j<h-1; j++ ) O[j]=I[j-1]+p*I[j]+I[j+1];
    O[j]=I[j-1]+(1+p)*I[j];
  }
  #undef C4
}

// convTri1( A, B, ns[0], ns[1], d, p, s );
// convConst('convTri1',I,12/r/(r+2)-2,s);
// convolve I by [1 p 1] filter (uses SSE)
void ColorChannel::convTri1( float *I, float *O, int h, int w, int d, float p, int s ) {
  const float nrm = 1.0f/((p+2)*(p+2)); int i, j, h0=h-(h%4);
  float *Il, *Im, *Ir, *T=(float*) malloc(h*sizeof(float));
  for( int d0=0; d0<d; d0++ ) for( i=s/2; i<w; i+=s ) {
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
Mat ColorChannel::rgbConvert(Mat I)
{
    Mat result;
    if (this->colorSpaceType == "luv")
        result = rgb2luv(I);
    else
        if (this->colorSpaceType == "hsv")
           result = rgb2hsv(I);
        else
            if (this->colorSpaceType == "gray")
                result = rgb2gray(I);
    return result;
}

Mat ColorChannel::rgb2luv(Mat I)
{
    Mat result;
    const double y0 = ((6.0/29)*(6.0/29)*(6.0/29));
    const double a = ((29.0/3)*(29.0/3)*(29.0/3));
    const double un = 0.197833, vn = 0.468331;
    double mr[3] = {0.430574, 0.222015, 0.020183};
    double mg[3] = {0.341550, 0.706655, 0.129553};
    double mb[3] = {0.178325, 0.071330, 0.939180};
    double maxi = 1.0/270, minu = -88*maxi, minv = -134*maxi;
    static double lTable[1064];
    double x, y, z, l;
    for (int i=0; i<1025; i++)
    {
        y = i/1024.0;
        if (y>y0)
            l = 116*pow(y,1.0/3.0)-16;
        else
            l = y * a;
        lTable[i] = l * maxi;
    }
    for (int i=1024; i<1064; i++)
        lTable[i] = lTable[i-1];
    for (int i = 0; i < I.rows*I.cols*I.channels(); i = i+3)
    {
        uchar b = I.data[i*I.step];
        uchar g = I.data[i*I.step+1];
        uchar r = I.data[i*I.step+2];
        x = mr[0]*r + mg[0]*g + mb[0]*b;
        y = mr[1]*r + mg[1]*g + mb[1]*b;
        z = mr[2]*r + mg[2]*g + mb[2]*b;
        l = lTable[(int)(y*1024)];
        result.data[i*I.step] = l; //l
        z = 1/(x + 15*y + 3*z + 1e-35);
        result.data[i*I.step+1] = l * (13*4*x*z - 13*un) - minu; //u
        result.data[i*I.step+2] = l * (13*9*y*z - 13*vn) - minv; //v
    }
    return result;
}

Mat ColorChannel::rgb2hsv(Mat I)
{
    //image will be represented as h,s,v
    Mat result;
    uchar h, s, v, maxValue=0, minValue=0;
    for (int i = 0; i < I.rows*I.cols*I.channels(); i = i+3)
    {
        uchar b = I.data[i*I.step];
        uchar g = I.data[i*I.step+1];
        uchar r = I.data[i*I.step+2];
        if (r == g && g == b)
        {
            result.data[i*I.step]   = 0; //hue
            result.data[i*I.step+1] = 0; //saturation
            result.data[i*I.step+2] = r; //value
        }
        else
        {
          if (r>=g && g>=b)
          {
              maxValue = r;
              if (g<b)
                  minValue = g;
              else
                  minValue = b;
              h = (g-b)/(maxValue-minValue)+6;
              if(h>=6)
                  h-=6;
          }
          else
            if(g>=r && g>=b)
            {
                maxValue = g;
                if(r<b)
                    minValue = r;
                else
                    minValue = b;
                h = (b-r)/(maxValue-minValue)+2;
            }
            else
            {
                maxValue = b;
                if (r<g)
                    minValue = r;
                else
                    minValue = g;
                h = (r-g)/(maxValue-minValue)+4;
            }
          result.data[i*I.step] = h*(1/6.0);
          result.data[i*I.step+1] = 1-minValue/maxValue;
          result.data[i*I.step+2] = maxValue;
        }
    }
    return result;
}

Mat ColorChannel::rgb2gray (Mat I)
{
    Mat grayImage;
    double redMultiplier   = 0.2989360213;
    double greenMultiplier = 0.5870430745;
    double blueMultiplier  = 0.1140209043;
    for (int i = 0; i < I.rows*I.cols*I.channels(); i = i+3)
    {
        grayImage.data[i] = I.data[i*I.step]*blueMultiplier +
        I.data[i*I.step+1]*greenMultiplier + I.data[i*I.step]*redMultiplier;
    }
    return grayImage;
}
