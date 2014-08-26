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

float* cvMat2floatArray(cv::Mat source, int channels)
{
	float* result = (float*)malloc(source.rows*source.cols*channels*sizeof(float));;
	float* tempFloat;
	cv::Mat tempMat;
	int resultIndex=0;

	//first, we need the type conversion
	if (channels == 1)
		source.convertTo(tempMat, CV_32FC1, 1.0/255.0);
	else
		source.convertTo(tempMat, CV_32FC3, 1.0/255.0);

	tempFloat = (float*)tempMat.data;

	//the next step is changing the way the rows, columns and channels are arranged
	for (int channel=0; channel < channels; channel++)
		for (int column=0; column < source.cols; column++)
			for (int row=0; row < source.rows; row++)
				result[resultIndex++] = tempFloat[column*channels + row*source.cols*channels + channel];

	return result;
}

cv::Mat floatArray2cvMat(float* source, int rows, int cols, int channels)
{
	cv::Mat result;
	if (channels == 1)
		result.convertTo(result, CV_32FC1);
	else
		result.convertTo(result, CV_32FC3);
	float* tempFloat = (float*)malloc(rows*cols*channels*sizeof(float));
	int tempIndex=0;

	for (int channel=0; channel < channels; channel++)
		for (int column=0; column < cols; column++)
			for (int row=0; row < rows; row++)
				tempFloat[column*channels + row*cols*channels +channel] = source[tempIndex++];

	result.data = (uchar*)tempFloat;

	return result;
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
//aka convTri
void triangleFilterConvolution( float *I, float *O, int h, int w, int d, int r, int s )
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
	
	float *T=(float*) alMalloc(2*h1*sizeof(float),16), *U=T+h1;

	while(d-- > 0) 
	{
		// initialize T and U
		for(j=0; j<h0; j+=4) 
		{
			STR(U[j], STR(T[j], LDu(I[j])));
		}
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
	alFree(T);
}

// convolve one column of I by [1 p 1] filter (uses SSE)
void convTri1Y( float *I, float *O, int h, float p, int s ) {
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
void convTri1( float *I, float *O, int h, int w, int d, float p, int s )
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

cv::Mat convolution(cv::Mat source, int channels, int radius, int s, int flag)
{

	// debug
	std::cout << "inside convolution" << std::endl;

	float* O = (float*)malloc(source.rows/s*source.cols/s*channels*sizeof(float));

	// debug
	std::cout << "inside convolution, after first malloc" << std::endl;

	float* I;
	I = cvMat2floatArray(source, channels);

	// debug
	std::cout << "inside convolution, after first conversion" << std::endl;

	switch(flag)
	{
		case CONV_TRI: 	

					// debug
					std::cout << "inside convolution, before CONV_TRI" << std::endl;

					triangleFilterConvolution(I, O, source.rows, source.cols, channels, radius, s);

					// debug
					std::cout << "inside convolution, after CONV_TRI" << std::endl;

					break;
		case CONV_TRI1: 
					int p = 12/radius/(radius+2)-2;
					convTri1(I, O, source.rows, source.cols, channels, p, s);
					break;
	}

	// debug
	std::cout << "inside convolution, before last conversion" << std::endl;

	cv::Mat result;
	result = floatArray2cvMat(O, source.rows, source.cols, channels); 

	// debug
	std::cout << "inside convolution, after last conversion" << std::endl;

	return result;
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
void imResample(float *A, float *B, int ha, int hb, int wa, int wb, int d, float r ) {
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

// I1=imResampleMex(I,sz1(1),sz1(2),1);
cv::Mat resample(cv::Mat source, int h, int w, float nrm, int channels)
{

	// debug
	//std::cout << "inside resample" << std::endl;

	cv::Mat result;

	// debug
	//std::cout << "inside resample, before conversion number one" << std::endl;

	float *I = cvMat2floatArray(source, channels);

	// debug
	//std::cout << "inside resample, after conversion number one" << std::endl;

	float *O = (float*)malloc(h*w*channels*sizeof(float));

	/* //debug
	cv::Mat testMat = floatArray2cvMat(I, source.rows, source.cols);
	cv::imshow("resample input", source);
	cv::imshow("resample input2", testMat);
	cv::waitKey();
	*/

	// debug
	//std::cout << "before imResample, source.rows = " << source.rows << ", h = " << h << ", source.cols = " << source.cols << ", w = " << w << std::endl;

	// resample((float*)A, (float*)B, ns[0], ms[0], ns[1], ms[1], nCh, nrm);
	// ns = (int*) mxGetDimensions(prhs[0]);
	// nCh=(nDims==2) ? 1 : ns[2];
	// nrm=(double)mxGetScalar(prhs[3]);
	// ms[0]=(int)mxGetScalar(prhs[1]); ms[1]=(int)mxGetScalar(prhs[2]); ms[2]=nCh;
	imResample(I, O, source.rows, h, source.cols, w, channels, nrm);

	// debug
	//std::cout << "inside resample, after imResample" << std::endl;

	result = floatArray2cvMat(O, h, w, channels); // three channels
	return result;
}