#include "ColorChannel.h"

//convolutions taken from convConst.cpp
Mat convolution(Mat source, int radius, int s, int flag)
{
	float* I = source.data;
	float* O;
	Mat result;

	switch(flag)
	{
		case CONV_TRI: convTri(I, O, source.rows, source.cols, source.dims, radius, s);
		case CONV_TRI1:
	}
	
	result.data = O;
	return result;
}

// convolve I by a 2rx1 triangle filter (uses SSE)
void convTri( float *I, float *O, int h, int w, int d, int r, int s ) {
  r++; float nrm = 1.0f/(r*r*r*r); int i, j, k=(s-1)/2, h0, h1, w0;
  if(h%4==0) h0=h1=h; else { h0=h-(h%4); h1=h0+4; } w0=(w/s)*s;
  float *T=(float*) alMalloc(2*h1*sizeof(float),16), *U=T+h1;
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
      k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    }
    I+=w*h;
  }
  alFree(T);
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

Mat rgb2luv(Mat I)
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

Mat rgb2hsv(Mat I)
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

Mat rgb2gray (Mat I)
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
