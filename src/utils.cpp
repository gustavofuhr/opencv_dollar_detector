#include "utils.h"

// wrapper functions if compiling from C/C++
inline void wrError(const char *errormsg) { throw errormsg; }
inline void* wrCalloc( size_t num, size_t size ) { return calloc(num,size); }
inline void* wrMalloc( size_t size ) { return malloc(size); }
inline void wrFree( void * ptr ) { free(ptr); }

/************************************************************************************************************/

// platform independent aligned memory allocation (see also alFree)
void* alMalloc( size_t size, int alignment ) {
  const size_t pSize = sizeof(void*), a = alignment-1;
  void *raw = wrMalloc(size + a + pSize);
  void *aligned = (void*) (((size_t) raw + pSize + a) & ~a);
  *(void**) ((size_t) aligned-pSize) = raw;
  return aligned;
}

// platform independent alignned memory de-allocation (see also alMalloc)
void alFree(void* aligned) {
  void* raw = *(void**)((char*)aligned-sizeof(void*));
  wrFree(raw);
}

/************************************************************************************************************/

// changes the way the rows, columns and channels are arranged to match the way piotr toolbox's procedures operate on
void cvMat2floatArray(cv::Mat source, float* result, int channels)
{
	float* tempFloat = (float*)source.data;
	int resultIndex=0;
	
	for (int channel=0; channel < channels; channel++)
		for (int column=0; column < source.cols; column++)
			for (int row=0; row < source.rows; row++)
				result[resultIndex++] = tempFloat[column*channels + row*source.cols*channels + channel];
}
// */

cv::Mat floatArray2cvMat(float* source, int rows, int cols, int channels)
{
  int type; 
  if (channels == 1)
    type = CV_32FC1;
  else
    type = CV_32FC3;
  
  float* tempFloat = (float*)malloc(rows*cols*channels*sizeof(float));
  int tempIndex=0;

  for (int channel=0; channel < channels; channel++)
    for (int column=0; column < cols; column++)
      for (int row=0; row < rows; row++)
        tempFloat[column*channels + row*cols*channels + channel] = source[tempIndex++];
  
  cv::Mat result(rows, cols, type);
  memcpy(result.data, tempFloat, rows*cols*channels*sizeof(float));

  free(tempFloat);

  return result;
}
// */

void features2floatArray (Info features, float* result, int rows, int cols, int colorChannels, int magChannels, int histChannels)
{
  int i, resultIndex=0;

  std::vector<cv::Mat> rgb;
  cv::split(features.image, rgb);
  for (int j=0; j < rgb.size(); j++)
  {
    cv::Mat tempMat;
    cv::transpose(rgb[j], tempMat);
    float* tempFloat = (float*)tempMat.data;
    for (i=0; i < rows*cols; i++)
      result[resultIndex++] = tempFloat[i];
  }

  cv::Mat tempMag;
  cv::transpose(features.gradientMagnitude, tempMag);
  float *floatMag = (float*)tempMag.data;
  for (i=0; i < rows*cols*magChannels; i++)
    result[resultIndex++] = floatMag[i];

  for (int j=0; j < histChannels; j++)
  {
    cv::Mat tempHist;
    cv::transpose(features.gradientHistogram[j], tempHist);
    float *floatHist = (float*)tempHist.data;
    for (i=0; i < rows*cols; i++)
      result[resultIndex++] = floatHist[i];
  }
}

/************************************************************************************************************/
// Convolutions

// convolve one column of I by a 2rx1 triangle filter
void convTriY( float *I, float *O, int h, int r, int s ) {
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
    k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    for( i=1; i<w0; i++ ) {
      float *Il=I+(i-1-r)*h; if(i<=r) Il=I+(r-i)*h; float *Im=I+(i-1)*h;
      float *Ir=I+(i-1+r)*h; if(i>w-r) Ir=I+(2*w-r-i)*h;
      for( j=0; j<h0; j+=4 ) {
        INC(T[j],ADD(LDu(Il[j]),LDu(Ir[j]),MUL(-2,LDu(Im[j]))));
        INC(U[j],MUL(nrm,LD(T[j])));
      }
      for( j=h0; j<h; j++ ) U[j]+=nrm*(T[j]+=Il[j]+Ir[j]-2*Im[j]);
      k++; if(k==s) { k=0; convTriY(U,O,h,r-1,s); O+=h/s; }
    }
    I+=w*h;
  }
  alFree(T);
}

// convolve one column of I by a [1 p 1] filter (uses SSE)
void convTri1Y( float *I, float *O, int h, float p, int s ) {
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

// convolve I by a [1 p 1] filter (uses SSE)
void convTri1( float *I, float *O, int h, int w, int d, float p, int s ) {
  const float nrm = 1.0f/((p+2)*(p+2)); int i, j, h0=h-(h%4);
  float *Il, *Im, *Ir, *T=(float*) alMalloc(h*sizeof(float),16);
  for( int d0=0; d0<d; d0++ ) for( i=s/2; i<w; i+=s ) {
    Il=Im=Ir=I+i*h+d0*h*w; if(i>0) Il-=h; if(i<w-1) Ir+=h;
    for( j=0; j<h0; j+=4 )
      STR(T[j],MUL(nrm,ADD(ADD(LDu(Il[j]),MUL(p,LDu(Im[j]))),LDu(Ir[j]))));
    for( j=h0; j<h; j++ ) T[j]=nrm*(Il[j]+p*Im[j]+Ir[j]);
    convTri1Y(T,O,h,p,s); O+=h/s;
  }
  alFree(T);
}

void convolution(float* source, float* result, int rows, int cols, int channels, int radius, int s)
{
  if (radius > 0 && radius <= 1)
  {
    int p = 12/radius/(radius+2)-2;
    convTri1(source, result, rows, cols, channels, p, s);
  }
  else
  {
    if (radius < 0)
      std::cout << " # Convolution error: invalid radius!\n";
    else
      convTri(source, result, rows, cols, channels, radius, s);
	}
}

/************************************************************************************************************/
// image resampling functions

// compute interpolation values for single column for resampling
void resampleCoef( int ha, int hb, int &n, int *&yas, int *&ybs, float *&wts, int bd[2], int pad=0 )
{
  const float s = float(hb)/float(ha), sInv = 1/s; float wt, wt0=float(1e-3)*s;
  bool ds=ha>hb; int nMax; bd[0]=bd[1]=0;
  if(ds) { n=0; nMax=ha+(pad>2 ? pad : 2)*hb; } else { n=nMax=hb; }
  // initialize memory
  wts = (float*)alMalloc(nMax*sizeof(float),16);
  yas = (int*)alMalloc(nMax*sizeof(int),16);
  ybs = (int*)alMalloc(nMax*sizeof(int),16);
  if( ds ) for( int yb=0; yb<hb; yb++ ) {
    // create coefficients for downsampling
    float ya0f=yb*sInv, ya1f=ya0f+sInv, W=0;
    int ya0=int(ceil(ya0f)), ya1=int(ya1f), n1=0;
    for( int ya=ya0-1; ya<ya1+1; ya++ ) {
      wt=s; if(ya==ya0-1) wt=(ya0-ya0f)*s; else if(ya==ya1) wt=(ya1f-ya1)*s;
      if(wt>wt0 && ya>=0) { ybs[n]=yb; yas[n]=ya; wts[n]=wt; n++; n1++; W+=wt; }
    }
    if(W>1) for( int i=0; i<n1; i++ ) wts[n-n1+i]/=W;
    if(n1>bd[0]) bd[0]=n1;
    while( n1<pad ) { ybs[n]=yb; yas[n]=yas[n-1]; wts[n]=0; n++; n1++; }
  } else for( int yb=0; yb<hb; yb++ ) {
    // create coefficients for upsampling
    float yaf = (float(.5)+yb)*sInv-float(.5); int ya=(int) floor(yaf);
    wt=1; if(ya>=0 && ya<ha-1) wt=1-(yaf-ya);
    if(ya<0) { ya=0; bd[0]++; } if(ya>=ha-1) { ya=ha-1; bd[1]++; }
    ybs[yb]=yb; yas[yb]=ya; wts[yb]=wt;
  }
}

// resample A using bilinear interpolation and and store result in B
void resample(float *A, float *B, int ha, int hb, int wa, int wb, int d, float r ) {
  int hn, wn, x, x1, y, z, xa, xb, ya; float *A0, *A1, *A2, *A3, *B0, wt, wt1;
  float *C = (float*) alMalloc((ha+4)*sizeof(float),16); 
  for(y=ha; y<ha+4; y++) C[y]=0;
  bool sse = !(size_t(A)&15) && !(size_t(B)&15);
  // get coefficients for resampling along w and h
  int *xas, *xbs, *yas, *ybs; float *xwts, *ywts; int xbd[2], ybd[2];
  resampleCoef( wa, wb, wn, xas, xbs, xwts, xbd, 0 );
  resampleCoef( ha, hb, hn, yas, ybs, ywts, ybd, 4 );
  if( wa==2*wb ) r/=2; if( wa==3*wb ) r/=3; if( wa==4*wb ) r/=4;
  r/=float(1+1e-6); for( y=0; y<hn; y++ ) ywts[y] *= r;
  // resample each channel in turn
  for( z=0; z<d; z++ ) for( x=0; x<wb; x++ ) {
    if(x==0) x1=0; xa=xas[x1]; xb=xbs[x1]; wt=xwts[x1]; wt1=1-wt; y=0;
    A0=A+z*ha*wa+xa*ha; A1=A0+ha, A2=A1+ha, A3=A2+ha; B0=B+z*hb*wb+xb*hb;
    // variables for SSE (simple casts to float)
    float *Af0, *Af1, *Af2, *Af3, *Bf0, *Cf, *ywtsf, wtf, wt1f;
    Af0=(float*) A0; Af1=(float*) A1; Af2=(float*) A2; Af3=(float*) A3;
    Bf0=(float*) B0; Cf=(float*) C;
    ywtsf=(float*) ywts; wtf=(float) wt; wt1f=(float) wt1;
    // resample along x direction (A -> C)
    #define FORs(X) if(sse) for(; y<ha-4; y+=4) STR(Cf[y],X);
    #define FORr(X) for(; y<ha; y++) C[y] = X;
    if( wa==2*wb ) {
      FORs( ADD(LDu(Af0[y]),LDu(Af1[y])) );
      FORr( A0[y]+A1[y] ); x1+=2;
    } else if( wa==3*wb ) {
      FORs( ADD(LDu(Af0[y]),LDu(Af1[y]),LDu(Af2[y])) );
      FORr( A0[y]+A1[y]+A2[y] ); x1+=3;
    } else if( wa==4*wb ) {
      FORs( ADD(LDu(Af0[y]),LDu(Af1[y]),LDu(Af2[y]),LDu(Af3[y])) );
      FORr( A0[y]+A1[y]+A2[y]+A3[y] ); x1+=4;
    } else if( wa>wb ) {
      int m=1; while( x1+m<wn && xb==xbs[x1+m] ) m++; float wtsf[4];
      for( int x0=0; x0<(m<4?m:4); x0++ ) wtsf[x0]=float(xwts[x1+x0]);
      #define U(x) MUL( LDu(*(Af ## x + y)), SET(wtsf[x]) )
      #define V(x) *(A ## x + y) * xwts[x1+x]
      if(m==1) { FORs(U(0));                     FORr(V(0)); }
      if(m==2) { FORs(ADD(U(0),U(1)));           FORr(V(0)+V(1)); }
      if(m==3) { FORs(ADD(U(0),U(1),U(2)));      FORr(V(0)+V(1)+V(2)); }
      if(m>=4) { FORs(ADD(U(0),U(1),U(2),U(3))); FORr(V(0)+V(1)+V(2)+V(3)); }
      #undef U
      #undef V
      for( int x0=4; x0<m; x0++ ) {
        A1=A0+x0*ha; wt1=xwts[x1+x0]; Af1=(float*) A1; wt1f=float(wt1); y=0;
        FORs(ADD(LD(Cf[y]),MUL(LDu(Af1[y]),SET(wt1f)))); FORr(C[y]+A1[y]*wt1);
      }
      x1+=m;
    } else {
      bool xBd = x<xbd[0] || x>=wb-xbd[1]; x1++;
      if(xBd) memcpy(C,A0,ha*sizeof(float));
      if(!xBd) FORs(ADD(MUL(LDu(Af0[y]),SET(wtf)),MUL(LDu(Af1[y]),SET(wt1f))));
      if(!xBd) FORr( A0[y]*wt + A1[y]*wt1 );
    }
    #undef FORs
    #undef FORr
    // resample along y direction (B -> C)
    if( ha==hb*2 ) {
		float r2 = r/2; int k=((~((size_t) B0) + 1) & 15)/4; y=0;
		for( ; y<k; y++ )  B0[y]=(C[2*y]+C[2*y+1])*r2;
		if(sse) for(; y<hb-4; y+=4) STR(Bf0[y],MUL((float)r2,_mm_shuffle_ps(ADD(
			LDu(Cf[2*y]),LDu(Cf[2*y+1])),ADD(LDu(Cf[2*y+4]),LDu(Cf[2*y+5])),136)));
		for( ; y<hb; y++ ) 
			B0[y]=(C[2*y]+C[2*y+1])*r2;
    } else if( ha==hb*3 ) {
      for(y=0; y<hb; y++) B0[y]=(C[3*y]+C[3*y+1]+C[3*y+2])*(r/3);
    } else if( ha==hb*4 ) {
      for(y=0; y<hb; y++) B0[y]=(C[4*y]+C[4*y+1]+C[4*y+2]+C[4*y+3])*(r/4);
    } else if( ha>hb ) {
      y=0;
      //if( sse && ybd[0]<=4 ) for(; y<hb; y++) // Requires SSE4
      //  STR1(Bf0[y],_mm_dp_ps(LDu(Cf[yas[y*4]]),LDu(ywtsf[y*4]),0xF1));
      #define U(o) C[ya+o]*ywts[y*4+o]
      if(ybd[0]==2) for(; y<hb; y++) { ya=yas[y*4]; B0[y]=U(0)+U(1); }
      if(ybd[0]==3) for(; y<hb; y++) { ya=yas[y*4]; B0[y]=U(0)+U(1)+U(2); }
      if(ybd[0]==4) for(; y<hb; y++) { ya=yas[y*4]; B0[y]=U(0)+U(1)+U(2)+U(3); }
      if(ybd[0]>4)  for(; y<hn; y++) { B0[ybs[y]] += C[yas[y]] * ywts[y]; }
      #undef U
    } else {
      for(y=0; y<ybd[0]; y++) B0[y] = C[yas[y]]*ywts[y];
      for(; y<hb-ybd[1]; y++) B0[y] = C[yas[y]]*ywts[y]+C[yas[y]+1]*(r-ywts[y]);
      for(; y<hb; y++)        B0[y] = C[yas[y]]*ywts[y];
    }
  }
  alFree(xas); alFree(xbs); alFree(xwts); alFree(C);
  alFree(yas); alFree(ybs); alFree(ywts);
}

/************************************************************************************************************/
// color space conversions using OpenCV functions
/*
cv::Mat rgbConvert(cv::Mat I, int colorSpace)
{
    cv::Mat result;
    cv::Mat uChar;

    if (colorSpace == LUV)
    {
        I.convertTo(uChar, CV_8UC3, 255.0);
        cvtColor(uChar, result, CV_BGR2Luv);
        result.convertTo(result, CV_32FC3, 1.0/255.0);
    }
    else
        if (colorSpace == HSV)
           cvtColor(I, result, CV_BGR2HSV);
        else
            if (colorSpace == GRAY)
                cvtColor(I, result, CV_BGR2GRAY);
            else
                return I;

    return result;
}

/************************************************************************************************************/
// colorspace conversions using Dóllar's mex file (rgbConvertMex.cpp)

// Constants for rgb2luv conversion and lookup table for y-> l conversion
template<class oT> oT* rgb2luv_setup( oT z, oT *mr, oT *mg, oT *mb,
  oT &minu, oT &minv, oT &un, oT &vn )
{
  // set constants for conversion
  const oT y0=(oT) ((6.0/29)*(6.0/29)*(6.0/29));
  const oT a= (oT) ((29.0/3)*(29.0/3)*(29.0/3));
  un=(oT) 0.197833; vn=(oT) 0.468331;
  mr[0]=(oT) 0.430574*z; mr[1]=(oT) 0.222015*z; mr[2]=(oT) 0.020183*z;
  mg[0]=(oT) 0.341550*z; mg[1]=(oT) 0.706655*z; mg[2]=(oT) 0.129553*z;
  mb[0]=(oT) 0.178325*z; mb[1]=(oT) 0.071330*z; mb[2]=(oT) 0.939180*z;
  oT maxi=(oT) 1.0/270; minu=-88*maxi; minv=-134*maxi;
  // build (padded) lookup table for y->l conversion assuming y in [0,1]
  static oT lTable[1064]; static bool lInit=false;
  if( lInit ) return lTable; oT y, l;
  for(int i=0; i<1025; i++) {
    y = (oT) (i/1024.0);
    l = y>y0 ? 116*(oT)pow((double)y,1.0/3.0)-16 : y*a;
    lTable[i] = l*maxi;
  }
  for(int i=1025; i<1064; i++) lTable[i]=lTable[i-1];
  lInit = true; return lTable;
}

// Convert from rgb to luv
template<class iT, class oT> void rgb2luv( iT *I, oT *J, int n, oT nrm ) {
  oT minu, minv, un, vn, mr[3], mg[3], mb[3];
  oT *lTable = rgb2luv_setup(nrm,mr,mg,mb,minu,minv,un,vn);
  oT *L=J, *U=L+n, *V=U+n; iT *R=I, *G=R+n, *B=G+n;
  for( int i=0; i<n; i++ ) {
    oT r, g, b, x, y, z, l;
    r=(oT)*R++; g=(oT)*G++; b=(oT)*B++;
    x = mr[0]*r + mg[0]*g + mb[0]*b;
    y = mr[1]*r + mg[1]*g + mb[1]*b;
    z = mr[2]*r + mg[2]*g + mb[2]*b;
    l = lTable[(int)(y*1024)];
    *(L++) = l; z = 1/(x + 15*y + 3*z + (oT)1e-35);
    *(U++) = l * (13*4*x*z - 13*un) - minu;
    *(V++) = l * (13*9*y*z - 13*vn) - minv;
  }
}

// Convert from rgb to luv using sse
template<class iT> void rgb2luv_sse( iT *I, float *J, int n, float nrm ) {
  const int k=256; float R[k], G[k], B[k];
  if( (size_t(R)&15||size_t(G)&15||size_t(B)&15||size_t(I)&15||size_t(J)&15)
    || n%4>0 ) { rgb2luv(I,J,n,nrm); return; }
  int i=0, i1, n1; float minu, minv, un, vn, mr[3], mg[3], mb[3];
  float *lTable = rgb2luv_setup(nrm,mr,mg,mb,minu,minv,un,vn);
  while( i<n ) {
    n1 = i+k; if(n1>n) n1=n; float *J1=J+i; float *R1, *G1, *B1;
    // convert to floats (and load input into cache)
    if( typeid(iT) != typeid(float) ) {
      R1=R; G1=G; B1=B; iT *Ri=I+i, *Gi=Ri+n, *Bi=Gi+n;
      for( i1=0; i1<(n1-i); i1++ ) {
        R1[i1] = (float) *Ri++; G1[i1] = (float) *Gi++; B1[i1] = (float) *Bi++;
      }
    } else { R1=((float*)I)+i; G1=R1+n; B1=G1+n; }
    // compute RGB -> XYZ
    for( int j=0; j<3; j++ ) {
      __m128 _mr, _mg, _mb, *_J=(__m128*) (J1+j*n);
      __m128 *_R=(__m128*) R1, *_G=(__m128*) G1, *_B=(__m128*) B1;
      _mr=SET(mr[j]); _mg=SET(mg[j]); _mb=SET(mb[j]);
      for( i1=i; i1<n1; i1+=4 ) *(_J++) = ADD( ADD(MUL(*(_R++),_mr),
        MUL(*(_G++),_mg)),MUL(*(_B++),_mb));
    }
    { // compute XZY -> LUV (without doing L lookup/normalization)
      __m128 _c15, _c3, _cEps, _c52, _c117, _c1024, _cun, _cvn;
      _c15=SET(15.0f); _c3=SET(3.0f); _cEps=SET(1e-35f);
      _c52=SET(52.0f); _c117=SET(117.0f), _c1024=SET(1024.0f);
      _cun=SET(13*un); _cvn=SET(13*vn);
      __m128 *_X, *_Y, *_Z, _x, _y, _z;
      _X=(__m128*) J1; _Y=(__m128*) (J1+n); _Z=(__m128*) (J1+2*n);
      for( i1=i; i1<n1; i1+=4 ) {
        _x = *_X; _y=*_Y; _z=*_Z;
        _z = RCP(ADD(_x,ADD(_cEps,ADD(MUL(_c15,_y),MUL(_c3,_z)))));
        *(_X++) = MUL(_c1024,_y);
        *(_Y++) = SUB(MUL(MUL(_c52,_x),_z),_cun);
        *(_Z++) = SUB(MUL(MUL(_c117,_y),_z),_cvn);
      }
    }
    { // perform lookup for L and finalize computation of U and V
      for( i1=i; i1<n1; i1++ ) J[i1] = lTable[(int)J[i1]];
      __m128 *_L, *_U, *_V, _l, _cminu, _cminv;
      _L=(__m128*) J1; _U=(__m128*) (J1+n); _V=(__m128*) (J1+2*n);
      _cminu=SET(minu); _cminv=SET(minv);
      for( i1=i; i1<n1; i1+=4 ) {
        _l = *(_L++);
        *(_U++) = SUB(MUL(_l,*_U),_cminu);
        *(_V++) = SUB(MUL(_l,*_V),_cminv);
      }
    }
    i = n1;
  }
}

// Convert from rgb to hsv
template<class iT, class oT> void rgb2hsv( iT *I, oT *J, int n, oT nrm ) {
  oT *H=J, *S=H+n, *V=S+n;
  iT *R=I, *G=R+n, *B=G+n;
  for(int i=0; i<n; i++) {
    const oT r=(oT)*(R++), g=(oT)*(G++), b=(oT)*(B++);
    oT h, s, v, minv, maxv;
    if( r==g && g==b ) {
      *(H++) = 0; *(S++) = 0; *(V++) = r*nrm; continue;
    } else if( r>=g && r>=b ) {
      maxv = r; minv = g<b ? g : b;
      h = (g-b)/(maxv-minv)+6; if(h>=6) h-=6;
    } else if( g>=r && g>=b ) {
      maxv = g; minv = r<b ? r : b;
      h = (b-r)/(maxv-minv)+2;
    } else {
      maxv = b; minv = r<g ? r : g;
      h = (r-g)/(maxv-minv)+4;
    }
    h*=(oT) (1/6.0); s=1-minv/maxv; v=maxv*nrm;
    *(H++) = h; *(S++) = s; *(V++) = v;
  }
}

// Convert from rgb to gray
template<class iT, class oT> void rgb2gray( iT *I, oT *J, int n, oT nrm ) {
  oT *GR=J; iT *R=I, *G=R+n, *B=G+n; int i;
  oT mr=(oT).2989360213*nrm, mg=(oT).5870430745*nrm, mb=(oT).1140209043*nrm;
  for(i=0; i<n; i++) *(GR++)=(oT)*(R++)*mr + (oT)*(G++)*mg + (oT)*(B++)*mb;
}

// Convert from rgb (double) to gray (float)
template<> void rgb2gray( double *I, float *J, int n, float nrm ) {
  float *GR=J; double *R=I, *G=R+n, *B=G+n; int i;
  double mr=.2989360213*nrm, mg=.5870430745*nrm, mb=.1140209043*nrm;
  for(i=0; i<n; i++) *(GR++) = (float) (*(R++)*mr + *(G++)*mg + *(B++)*mb);
}

// Copy and normalize only
template<class iT, class oT> void normalize( iT *I, oT *J, int n, oT nrm ) {
  for(int i=0; i<n; i++) *(J++)=(oT)*(I++)*nrm;
}

// Convert rgb to various colorspaces
void rgbConvert(float *I, float* J, int n, int d, int flag, float nrm) 
{
  //float *J = (float*) wrMalloc(n*(flag==0 ? (d==1?1:d/3) : d)*sizeof(float));
  int i, n1=d*(n<1000?n/10:100); float thr = float(1.001);
  if(flag>1 && nrm==1) for(i=0; i<n1; i++) if(I[i]>thr)
    wrError("For floats all values in I must be smaller than 1.");
  bool useSse = n%4==0;
  if( flag==2 && useSse )
    for(i=0; i<d/3; i++) rgb2luv_sse(I+i*n*3,(float*)(J+i*n*3),n,(float)nrm);
    //for(i=0; i<d/3; i++) rgb2luv(I+i*n*3,(float*)(J+i*n*3),n,(float)nrm);
  else if( (flag==0 && d==1) || flag==1 ) normalize(I,J,n*d,nrm);
  else if( flag==0 ) for(i=0; i<d/3; i++) rgb2gray(I+i*n*3,J+i*n*1,n,nrm);
  else if( flag==2 ) for(i=0; i<d/3; i++) rgb2luv(I+i*n*3,J+i*n*3,n,nrm);
  else if( flag==3 ) for(i=0; i<d/3; i++) rgb2hsv(I+i*n*3,J+i*n*3,n,nrm);
  else wrError("Unknown flag.");
  //return J;
}

/************************************************************************************************************/

// pad A by [pt,pb,pl,pr] and store result in B
void imPad(float *A, float *B, int h, int w, int d, int pt, int pb, int pl, int pr, int flag, float val)
{
  int h1=h+pt, hb=h1+pb, w1=w+pl, wb=w1+pr, x, y, z, mPad;
  int ct=0, cb=0, cl=0, cr=0;
  if(pt<0) { ct=-pt; pt=0; } if(pb<0) { h1+=pb; cb=-pb; pb=0; }
  if(pl<0) { cl=-pl; pl=0; } if(pr<0) { w1+=pr; cr=-pr; pr=0; }
  int *xs, *ys; x=pr>pl?pr:pl; y=pt>pb?pt:pb; mPad=x>y?x:y;
  bool useLookup = ((flag==2 || flag==3) && (mPad>h || mPad>w))
    || (flag==3 && (ct || cb || cl || cr ));
  // helper macro for padding
  #define PAD(XL,XM,XR,YT,YM,YB) \
  for(x=0;  x<pl; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XL+cl)*h+YT+ct]; \
  for(x=0;  x<pl; x++) for(y=pt; y<h1; y++) B[x*hb+y]=A[(XL+cl)*h+YM+ct]; \
  for(x=0;  x<pl; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XL+cl)*h+YB-cb]; \
  for(x=pl; x<w1; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XM+cl)*h+YT+ct]; \
  for(x=pl; x<w1; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XM+cl)*h+YB-cb]; \
  for(x=w1; x<wb; x++) for(y=0;  y<pt; y++) B[x*hb+y]=A[(XR-cr)*h+YT+ct]; \
  for(x=w1; x<wb; x++) for(y=pt; y<h1; y++) B[x*hb+y]=A[(XR-cr)*h+YM+ct]; \
  for(x=w1; x<wb; x++) for(y=h1; y<hb; y++) B[x*hb+y]=A[(XR-cr)*h+YB-cb];
  // build lookup table for xs and ys if necessary
  if( useLookup ) {
    xs = (int*) malloc(wb*sizeof(int)); int h2=(pt+1)*2*h;
    ys = (int*) malloc(hb*sizeof(int)); int w2=(pl+1)*2*w;
    if( flag==2 ) {
      for(x=0; x<wb; x++) { z=(x-pl+w2)%(w*2); xs[x]=z<w ? z : w*2-z-1; }
      for(y=0; y<hb; y++) { z=(y-pt+h2)%(h*2); ys[y]=z<h ? z : h*2-z-1; }
    } else if( flag==3 ) {
      for(x=0; x<wb; x++) xs[x]=(x-pl+w2)%w;
      for(y=0; y<hb; y++) ys[y]=(y-pt+h2)%h;
    }
  }
  // pad by appropriate value
  for( z=0; z<d; z++ ) {
    // copy over A to relevant region in B
    for( x=0; x<w-cr-cl; x++ )
      memcpy(B+(x+pl)*hb+pt,A+(x+cl)*h+ct,sizeof(float)*(h-ct-cb));
    // set boundaries of B to appropriate values
    if( flag==0 && val!=0 ) { // "constant"
      for(x=0;  x<pl; x++) for(y=0;  y<hb; y++) B[x*hb+y]=val;
      for(x=pl; x<w1; x++) for(y=0;  y<pt; y++) B[x*hb+y]=val;
      for(x=pl; x<w1; x++) for(y=h1; y<hb; y++) B[x*hb+y]=val;
      for(x=w1; x<wb; x++) for(y=0;  y<hb; y++) B[x*hb+y]=val;
    } else if( useLookup ) { // "lookup"
      PAD( xs[x], xs[x], xs[x], ys[y], ys[y], ys[y] );
    } else if( flag==1 ) {  // "replicate"
      PAD( 0, x-pl, w-1, 0, y-pt, h-1 );
    } else if( flag==2 ) { // "symmetric"
      PAD( pl-x-1, x-pl, w+w1-1-x, pt-y-1, y-pt, h+h1-1-y );
    } else if( flag==3 ) { // "circular"
      PAD( x-pl+w, x-pl, x-pl-w, y-pt+h, y-pt, y-pt-h );
    }
    A += h*w;  B += hb*wb;
  }
  if( useLookup ) { alFree(xs); alFree(ys); }
  #undef PAD
}

cv::Mat padImage(cv::Mat source, int channels, int *pad, int padSize, int type)
{
	cv::Mat result;
	int padTop, padBottom, padLeft, padRight;

	switch (padSize)
	{
		case 1:	
			padTop = padBottom = padLeft = padRight = pad[0];
		break;
		case 2:
			padTop = padBottom = pad[0];
			padLeft = padRight = pad[1];
		break;
		case 4:
			padTop = pad[0];
			padBottom = pad[1];
			padLeft = pad[2];
			padRight = pad[3];
		break;
		default:
			std::cout << " # padImage error: padSize must be 1, 2 or 4." << std::endl;
	}

	int newRows = source.rows+padTop+padBottom;
	int newCols = source.cols+padLeft+padRight;

	float *O, *I;
	O = (float*) malloc((newRows)*(newCols)*channels*sizeof(float));

	//I = cvImage2floatArray(source, channels);
  
  int imgIndex=0;
  I = (float*)malloc(source.rows*source.cols*channels*sizeof(float));
  std::vector<cv::Mat> splitChannels;
  cv::split(source, splitChannels);
  for (int j=0; j < splitChannels.size(); j++)
  {
    cv::Mat tempMat;
    cv::transpose(splitChannels[j], tempMat);
    float* tempFloat = (float*)tempMat.data;
    for (int k=0; k < source.rows*source.cols; k++)
        I[imgIndex++] = tempFloat[k];
  }
  // */

	double val = (double)type;

	imPad(I, O, source.rows, source.cols, channels, padTop, padBottom, padLeft, padRight, type, val);

  /*
  // old conversion
  float* tempdata = (float*)malloc(newRows*newCols*channels*sizeof(float));
  floatArray2cvData(O, tempdata, newRows, newCols, channels);
  int matType = CV_32FC1;
  if (channels == 3)
    matType = CV_32FC3;
  cv::Mat tempMat(newRows, newCols, matType);
  tempMat.data = (uchar*)tempdata;
  tempMat.copyTo(result);
  tempMat.release();
  free(tempdata);
  // old conversion */

  // improved conversion
  result = floatArray2cvMat(O, newRows, newCols, channels);

  free(I);
  free(O);

	return result;
}

/************************************************************************************************************/

std::vector<std::string> getDataSetFileNames(std::string directory)
{
  std::vector<std::string> result, temp;
  dirent* de;
  DIR* dp;

  dp = opendir(directory.empty() ? "." : directory.c_str());

  if (dp)
  {
    de = readdir(dp);
    while (de != NULL)
    {
      temp.push_back(std::string(de->d_name));
      de = readdir(dp);
    }
    closedir(dp);

    for (int i=0; i < temp.size(); i++)
    {
      if (temp[i].at(0) != '.')
        result.push_back(temp[i]);
    }

    std::sort(result.begin(), result.end());
  }

  return result;
}

/************************************************************************************************************/

cv::Point2f imagePoint2groundPlanePoint(float imageU, float imageV, float imageZ, cv::Mat homography)
{
  cv::Mat inverseH;
  invert(homography, inverseH);

  float x = inverseH.at<float>(0,0)*imageU +
        inverseH.at<float>(0,1)*imageV +  
        inverseH.at<float>(0,2)*imageZ;

  float y = inverseH.at<float>(1,0)*imageU +
        inverseH.at<float>(1,1)*imageV +
        inverseH.at<float>(1,2)*imageZ;

  float z = inverseH.at<float>(2,0)*imageU +
        inverseH.at<float>(2,1)*imageV +
        inverseH.at<float>(2,2)*imageZ;

  x = x/z;
  y = y/z;

  cv::Point2f result;
  result.x = x;
  result.y = y;

  return result;
} 

cv::Point worldPoint2imagePoint(float worldX, float worldY, float worldZ, cv::Mat homography)
{
  float u = homography.at<float>(0,0)*worldX +
        homography.at<float>(0,1)*worldY +  
        homography.at<float>(0,2)*worldZ;

  float v = homography.at<float>(1,0)*worldX +
        homography.at<float>(1,1)*worldY +
        homography.at<float>(1,2)*worldZ;

  float z = homography.at<float>(2,0)*worldX +
        homography.at<float>(2,1)*worldY +
        homography.at<float>(2,2)*worldZ;

  u = u/z;
  v = v/z;

  cv::Point result;
  result.x = u;
  result.y = v;

  return result;
}

double findWorldHeight(int u, int bottom_v, int top_v, cv::Mat_<float> &projection, cv::Mat_<float> &homography)
{
  double height;

  float p2 = projection.at<float>(1,2);
  float p3 = projection.at<float>(2,2);

  cv::Mat inverseH;
  invert(homography, inverseH);

  float h1 = inverseH.at<float>(2,0);
  float h2 = inverseH.at<float>(2,1);
  float h3 = inverseH.at<float>(2,2);

  height = -(bottom_v - top_v)/(p2*(h3 + h1*u + h2*bottom_v) - p3*top_v*(h3 + h1*u + h2*bottom_v));

  return height;
}

std::vector<cv::Point2f> findGroundPlaneAndImageIntersectionPoints(int imageWidth, int imageHeight, int boundingBoxImageWidth, int boundingBoxImageHeight, 
  float maxPedestrianHeight, cv::Mat_<float> &projection, cv::Mat_<float> &homography)
{
  std::vector<cv::Point2f> points(4);
  bool foundP1=false, foundP2=false;

  // finding top left point
  int curV = boundingBoxImageHeight;
  int curU = boundingBoxImageWidth/2;
  double curHeight = maxPedestrianHeight+1.0;
  while (curHeight > maxPedestrianHeight && curV <= imageHeight)
  {
    curHeight = findWorldHeight(curU, curV, curV-boundingBoxImageHeight, projection, homography);

    if (curHeight <= maxPedestrianHeight)
    {
      foundP1 = true;
      cv::Point2f p1 = imagePoint2groundPlanePoint(0, curV, 1.0, homography);
      points[0] = p1;

      // debug
      //std::cout << "p1_i=(" << 0 << "," << curV << "), p1_w=(" << p1.x << "," << p1.y << ")\n";
    }

    curV = curV + 1;   
  }

  // finding top right point
  curU = imageWidth - boundingBoxImageWidth/2;
  curV = boundingBoxImageHeight;
  curHeight = maxPedestrianHeight+1.0;
  while (curHeight > maxPedestrianHeight && curV <= imageHeight)
  {
    curHeight = findWorldHeight(imageWidth-boundingBoxImageWidth, curV, curV-boundingBoxImageHeight, projection, homography);

    if (curHeight <= maxPedestrianHeight)
    {
      foundP2 = true;
      cv::Point2f p2 = imagePoint2groundPlanePoint(imageWidth - boundingBoxImageWidth, curV, 1.0, homography);
      points[1] = p2;

      //std::cout << "p2_i=(" << imageWidth - boundingBoxImageWidth << "," << curV << "), p2_w=(" << p2.x << "," << p2.y << ")\n";
    }

    curV = curV + 1;   
  }

  // finding bottom left point
  cv::Point2f p3 = imagePoint2groundPlanePoint(0, imageHeight, 1.0, homography);
  points[2] = p3;
  //std::cout << "p3_i=(" << 0 << "," << imageHeight << "), p3_w=(" << p3.x << "," << p3.y << ")\n";

  // finding bottom right point
  cv::Point2f p4 = imagePoint2groundPlanePoint(imageWidth-boundingBoxImageWidth, imageHeight, 1.0, homography);
  points[3] = p4;
  //std::cout << "p4_i=(" << imageWidth-boundingBoxImageWidth << "," << imageHeight << "), p4_w=(" << p4.x << "," << p4.y << ")\n";

  return points;
}
 
std::vector<cv::Point2f> trimGroundPlanesBottomPoints(int imageWidth, int imageHeight, int boundingBoxImageWidth, int boundingBoxImageHeight, 
  double targetPedestrianWorldHeight, std::vector<cv::Point2f> groundPlaneLimits, cv::Mat_<float> &projection, cv::Mat_<float> &homography)
{
  std::vector<cv::Point2f> newGroundPlaneLimits(4);
  newGroundPlaneLimits[0] = groundPlaneLimits[0];
  newGroundPlaneLimits[1] = groundPlaneLimits[1];
  cv::Point topLeftLimitImagePoint = worldPoint2imagePoint(groundPlaneLimits[0].x, groundPlaneLimits[0].y, 0.0, homography);
  cv::Point topRightLimitImagePoint = worldPoint2imagePoint(groundPlaneLimits[1].x, groundPlaneLimits[1].y, 0.0, homography);
  double curBBworldHeight = targetPedestrianWorldHeight - 1.0;
  int curBBimageHeight = boundingBoxImageHeight;
  int octaves = 0;
  bool done = false;
  int curV = imageHeight;

  while (curV > topLeftLimitImagePoint.y && !done)
  {
    curBBworldHeight = findWorldHeight(0, curV, curV-curBBimageHeight, projection, homography);
    std::cout << "curBBworldHeight=" << curBBworldHeight << std::endl;
    if (curBBworldHeight >= targetPedestrianWorldHeight)
    {
      done = true;
      cv::Point2f newBottomLeft = imagePoint2groundPlanePoint(0, curV, 1.0, homography);
      newGroundPlaneLimits[2] = newBottomLeft;

      // debug
      std::cout << "after " << octaves << " octaves, bottomLeft=(" << newBottomLeft.x << "," << newBottomLeft.y << ")\n";  
      // debug */
    }
    else
    {
      octaves++;
      curBBimageHeight = curBBimageHeight * 2;
      if (curBBimageHeight > imageHeight)
      {
        octaves = 0;
        curBBimageHeight = boundingBoxImageHeight;
        curV = curV - 1;
      }
    }
  }
  if (!done)
    newGroundPlaneLimits[2] = groundPlaneLimits[2];
  else
    done = false;

  curBBworldHeight = targetPedestrianWorldHeight - 1.0;
  curBBimageHeight = boundingBoxImageHeight;
  octaves = 0;
  curV = imageHeight;
  while (curV > topRightLimitImagePoint.y && !done)
  {
    curBBworldHeight = findWorldHeight(imageWidth-boundingBoxImageWidth, curV, curV-curBBimageHeight, projection, homography);
    std::cout << "curBBworldHeight=" << curBBworldHeight << std::endl;
    if (curBBworldHeight >= targetPedestrianWorldHeight)
    {
      done = true;
      cv::Point2f newBottomRight = imagePoint2groundPlanePoint(imageWidth-boundingBoxImageWidth, curV, 1.0, homography);
      newGroundPlaneLimits[3] = newBottomRight;
      // debug
      std::cout << "after " << octaves << " octaves, bottomRight=(" << newBottomRight.x << "," << newBottomRight.y << ")\n";  
      // debug */
    }
    else
    {
      octaves++;
      curBBimageHeight = curBBimageHeight * 2;
      if (curBBimageHeight > imageHeight)
      {
        octaves = 0;
        curBBimageHeight = boundingBoxImageHeight;
        curV = curV - 1;
      }
    }
  }
  if (!done)
    newGroundPlaneLimits[3] = groundPlaneLimits[3];


  return newGroundPlaneLimits;
}

double findLastNecessaryScaleInAPoint(int u, int v, int imageHeight, int boundingBoxImageHeight, double maxPedestrianWorldHeight, 
  cv::Mat_<float> &projection, cv::Mat_<float> &homography)
{
  double lastScale;
  bool done = false;
  double curBBworldHeight = 0.0;
  int curBBimageHeight = boundingBoxImageHeight;

  while (curBBimageHeight < imageHeight && !done)
  {
    curBBworldHeight = findWorldHeight(u, v, v-curBBimageHeight, projection, homography);

    if (curBBworldHeight >= maxPedestrianWorldHeight)
    {
      done = true;
    }
    else
      curBBimageHeight++;
  }

  lastScale = (double)boundingBoxImageHeight/curBBimageHeight;

  //std::cout << "With curBBimageHeight=" << curBBimageHeight << " and curBBworldHeight=" << curBBworldHeight << ", lastScale=" << lastScale << std::endl;

  return lastScale;
}

int findNecessaryNumberOfOctaves(int imageWidth, int imageHeight, int boundingBoxImageWidth, int boundingBoxImageHeight, 
  double minPedestrianWorldHeight, double maxPedestrianWorldHeight, cv::Mat_<float> &projection, cv::Mat_<float> &homography)
{
  int necessaryOctaves;
  int octavesLeft=0, octavesRight=0;
  double curBBworldHeight = minPedestrianWorldHeight - 1.0;
  int curBBimageHeight = boundingBoxImageHeight;
  bool done = false;
  int curV = imageHeight;
  while (curBBimageHeight <= imageHeight && !done)
  {
    curBBworldHeight = findWorldHeight(0, curV, curV-curBBimageHeight, projection, homography);
    std::cout << "curBBworldHeight=" << curBBworldHeight << std::endl;
    if (curBBworldHeight >= minPedestrianWorldHeight)
    {
      done = true;
      
      if (curBBworldHeight > maxPedestrianWorldHeight)
        octavesLeft--;

      // debug
      std::cout << octavesLeft << " octaves on the bottom left corner\n";   
      // debug */
    }
    else
    {
      octavesLeft++;
      curBBimageHeight = curBBimageHeight * 2;
    }
  }
  done = false;

  curBBworldHeight = minPedestrianWorldHeight - 1.0;
  curBBimageHeight = boundingBoxImageHeight;
  curV = imageHeight;
  while (curBBimageHeight <= imageHeight && !done)
  {
    curBBworldHeight = findWorldHeight(imageWidth-boundingBoxImageWidth, curV, curV-curBBimageHeight, projection, homography);
    std::cout << "curBBworldHeight=" << curBBworldHeight << std::endl;
    if (curBBworldHeight >= minPedestrianWorldHeight)
    {
      done = true;

      if (curBBworldHeight > maxPedestrianWorldHeight)
        octavesRight--;

      // debug
      std::cout << octavesRight << " octaves on the bottom right corner\n";  
      // debug */
    }
    else
    {
      octavesRight++;
      curBBimageHeight = curBBimageHeight * 2;
    }
  }

  necessaryOctaves = octavesLeft;
  if (octavesRight > necessaryOctaves)
    necessaryOctaves = octavesRight;

  return necessaryOctaves;
}

cv::Mat world2image(cv::Mat &w_point, cv::Mat_<float> &P) {
  
  cv::Mat wp;
  wp   = P*w_point;

  for (int i = 0; i < w_point.cols; i++) {
    // normalize the point
    wp.at<float>(0,i) = wp.at<float>(0,i)/wp.at<float>(2,i);
    wp.at<float>(1,i) = wp.at<float>(1,i)/wp.at<float>(2,i);
  }

  cv::Mat res = wp(cv::Range(0,2), cv::Range::all());
  return res;
}

BoundingBox wcoord2bbox(cv::Point2f w_point, cv::Mat_<float> &P, float w_height, float aspect_ratio) {
  BoundingBox bb;

  // create the foot and head points
  cv::Mat_<float> w_feet_point(4, 1, 0.0);
  w_feet_point(0) = w_point.x;
  w_feet_point(1) = w_point.y;
  w_feet_point(2) = 0.0;
  w_feet_point(3) = 1.0;

  cv::Mat_<float> w_head_point(4, 1, 0.0);
  w_feet_point.copyTo(w_head_point);
  w_head_point(2) = w_height;

  // project them to create the bbox
  cv::Mat_<float> i_feet_point = world2image(w_feet_point, P);
  cv::Mat_<float> i_head_point = world2image(w_head_point, P);

  // compute the height and width in the image
  float i_height = fabs(i_feet_point(1) - i_head_point(1));
  float i_width  = i_height * aspect_ratio;

  // the central x will be the middle point between the line
  // from the feet to the head
  float ix_middle = (i_feet_point(0) + i_head_point(0))/2.0;

  bb.topLeftPoint.x = ix_middle - i_width/2.0;
  bb.topLeftPoint.y = i_head_point(1);
  bb.width          = i_width;
  bb.height         = i_height;
  bb.world_height   = w_height;

  return bb;
}

void print_fmatrix(const std::string &title, const cv::Mat &m) {

  std::cout << title << std::endl;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      std::cout << m.at<float>(i,j) << " ";
    }
    std::cout << std::endl;
  }

}

void print_dmatrix(const std::string &title, const cv::Mat &m) {

  std::cout << title << std::endl;
  for (int i = 0; i < m.rows; ++i) {
    for (int j = 0; j < m.cols; ++j) {
      std::cout << m.at<double>(i,j) << " ";
    }
    std::cout << std::endl;
  }

}

double log_base_n(double y, double base) 
{
  return log(y)/log(base);
}

/************************************************************************************************************/
