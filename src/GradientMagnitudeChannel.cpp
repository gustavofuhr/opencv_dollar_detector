#include "GradientMagnitudeChannel.h"

void GradientMagnitudeChannel::readGradientMagnitude(cv::FileNode gradMagNode)
{
	enabled = gradMagNode["enabled"];
	colorChannelIndex = gradMagNode["colorChn"];
	normalizationRadius = gradMagNode["normRad"];
	normalizationConstant = gradMagNode["normConst"];
	full = gradMagNode["full"];
}

// platform independent aligned memory allocation (see also alFree)
void* alMalloc( size_t size, int alignment ) {
  const size_t pSize = sizeof(void*), a = alignment-1;
  void *raw = malloc(size + a + pSize);
  void *aligned = (void*) (((size_t) raw + pSize + a) & ~a);
  *(void**) ((size_t) aligned-pSize) = raw;
  return aligned;
}

// platform independent alignned memory de-allocation (see also alMalloc)
void alFree(void* aligned) {
  void* raw = *(void**)((char*)aligned-sizeof(void*));
  free(raw);
}

// [M,O] = gradMag( I, channel, full ) - see gradientMag.m
std::vector<cv::Mat> GradientMagnitudeChannel::mGradMag(cv::Mat I, int channel)
{
  int c, d; 
	float *M, *O=0;
  M = (float*)alMalloc(I.rows*I.cols, 16);
	std::vector<cv::Mat> resultMatrix;

	//checkArgs procedure is called but it is not actually needed
	//probably just need to test some of the parameters, if that
	//checkArgs(nl,pl,nr,pr,1,2,3,3,&h,&w,&d,mxSINGLE_CLASS,(void**)&I);

	if (I.rows>=2 && I.cols>=2)
	{
		//for now, the gradMag procedure will be unchanged so we put the raw data of the image in a float pointer to be used there 
    cv::Mat floatMat;
    I.convertTo(floatMat, CV_32FC1, 1.0/255.0);
    float* If;
    //If = (float*)floatMat.data;
    If = (float*)I.data;

    cv::Mat test(I.rows, I.cols, I.type());
    test.data = (uchar*)If;

		if (channel>0 && channel<=I.dims)
		{
			If += I.rows*I.cols*(channel-1); 
			d=1;  
		}
		//i dont think this next if statement will be necessary
		/*if (nl>=2)
			Mat pl; //pl[0] = mxCreateMatrix3(h,w,1,mxSINGLE_CLASS,0,(void**)&M);*/

    /*cv::destroyAllWindows();
    cv::imshow("input image inside mGradMag", I);
    cv::imshow("float image inside mGradMag", floatMat);
    cv::imshow("test image inside mGradMag", test);
    cv::waitKey();
    cv::destroyAllWindows();*/

    std::cout << "inside mGradMag, before gradMag" << std::endl;

    cv::imshow("image passed to gradMag from mGradMag", test);
    cv::waitKey();

		//call to the actual function: 
		//void gradMag(float*, float*, float*, int, int, int, bool);
    gradMag(If, M, O, I.rows, I.cols, 3, full>0);

     std::cout << "inside mGradMag, after gradMag" << std::endl;

		//next, we assign the values of M and O to the matrix thats going to be returned
    cv::Mat matM;
    matM.data = (uchar*)M;
    cv::Mat matO;
    matO.data = (uchar*)O;
    resultMatrix.push_back(matM);
    resultMatrix.push_back(matO);
    int f = resultMatrix[0].cols;
	}
	else
	{
		//error: matrix I must be at least 2x2
	}
	return resultMatrix;
}

// gradMagNorm( M, S, norm ) - operates on M - see gradientMag.m
//gradientMex('gradientMagNorm',M,S,normConst);
// normalize gradient magnitude at each location (uses sse) 
cv::Mat GradientMagnitudeChannel::gradMagNorm(float* M, float* S, int h, int w) {
	cv::Mat resultingM;
	__m128 *_M, *_S, _norm; 
	int i=0, n=h*w, n4=n/4;
  	_S = (__m128*) S; 
	_M = (__m128*) M; 
	_norm = SET(normalizationConstant);
  	bool sse = !(size_t(M)&15) && !(size_t(S)&15);
  	if(sse) 
	{ 
		for(; i<n4; i++) 
			*_M++=MUL(*_M,RCP(ADD(*_S++,_norm))); 
		i*=4; 
	}
  	for(; i<n; i++) 
		M[i] /= (S[i] + normalizationConstant);
	resultingM.data = (uchar*)M;
	return resultingM;
}

// compute gradient magnitude and orientation at each location
// gradMag(I, M, O, h, w, d, full>0 );
void GradientMagnitudeChannel::gradMag(float *I, float *M, float *O, int h, int w, int d, bool full) {
  int x, y, y1, c, h4, s; float *Gx, *Gy, *M2; __m128 *_Gx, *_Gy, *_M2, _m;
  float *acost = acosTable(), acMult=10000.0f;
  // allocate memory for storing one column of output (padded so h4%4==0)
  h4=(h%4==0) ? h : h-(h%4)+4; s=d*h4*sizeof(float);
  M2=(float*) alMalloc(s,16); _M2=(__m128*) M2;
  Gx=(float*) alMalloc(s,16); _Gx=(__m128*) Gx;
  Gy=(float*) alMalloc(s,16); _Gy=(__m128*) Gy;
  // compute gradient magnitude and orientation for each column
  for( x=0; x<w; x++ ) {
    // compute gradients (Gx, Gy) with maximum squared magnitude (M2)
    for(c=0; c<d; c++) 
    {
      cv::destroyAllWindows();
      std::cout << "inside gradMag, before assignment from I to test.data" << std::endl;
      cv::Mat test(h, w, 16);
      test.data = (uchar*)I;
      std::cout << "inside gradMag, before printing image to be passed to grad1" << std::endl;
      cv::imshow("image passed to grad1 from gradMag", test);
      cv::waitKey();
      //the last call before the crash is x = 0, c = 2, but it crashes at c =1 if i try to print the input image
      std::cout << "inside gradMag, before calling grad1, x = " << x << ", c = " << c << std::endl;
      grad1( I+x*h+c*w*h, Gx+c*h4, Gy+c*h4, h, w, x );
      std::cout << "inside gradMag, after calling grad1" << std::endl;
      for( y=0; y<h4/4; y++ ) {
        y1=h4/4*c+y;
        _M2[y1]=ADD(MUL(_Gx[y1],_Gx[y1]),MUL(_Gy[y1],_Gy[y1]));
        if( c==0 ) continue; 
    _m = CMPGT( _M2[y1], _M2[y] );
        _M2[y] = OR( AND(_m,_M2[y1]), ANDNOT(_m,_M2[y]) );
        _Gx[y] = OR( AND(_m,_Gx[y1]), ANDNOT(_m,_Gx[y]) );
        _Gy[y] = OR( AND(_m,_Gy[y1]), ANDNOT(_m,_Gy[y]) );
      }
    }
    // compute gradient mangitude (M) and normalize Gx
    for( y=0; y<h4/4; y++ ) {
      _m = MMIN( RCPSQRT(_M2[y]), SET(1e10f) );
      _M2[y] = RCP(_m);
      if(O) _Gx[y] = MUL( MUL(_Gx[y],_m), SET(acMult) );
      if(O) _Gx[y] = XOR( _Gx[y], AND(_Gy[y], SET(-0.f)) );
    };
    memcpy( M+x*h, M2, h*sizeof(float) );
    // compute and store gradient orientation (O) via table lookup
    if( O!=0 ) for( y=0; y<h; y++ ) O[x*h+y] = acost[(int)Gx[y]];
    if( O!=0 && full ) {
      y1=((~size_t(O+x*h)+1)&15)/4; y=0;
      for( ; y<y1; y++ ) O[y+x*h]+=(Gy[y]<0)*PI;
      for( ; y<h-4; y+=4 ) STRu( O[y+x*h],
        ADD( LDu(O[y+x*h]), AND(CMPLT(LDu(Gy[y]),SET(0.f)),SET(PI)) ) );
      for( ; y<h; y++ ) O[y+x*h]+=(Gy[y]<0)*PI;
    }
  }
}

// compute x and y gradients for just one column (uses sse)
void GradientMagnitudeChannel::grad1( float *I, float *Gx, float *Gy, int h, int w, int x ) {
  /*std::cout << "inside grad1, before assignment from I to test.data" << std::endl;
  cv::Mat test(h, w, 16);
  test.data = (uchar*)I;
  std::cout << "inside grad1, before printing image" << std::endl;
  cv::imshow("testing grad1 input", test);
  cv::waitKey();
  std::cout << "inside grad1, after printing image" << std::endl;*/
  int y, y1; float *Ip, *In, r; __m128 *_Ip, *_In, *_G, _r;
  // compute column of Gx
  Ip=I-h; In=I+h; r=.5f;
  if(x==0) { r=1; Ip+=h; } else if(x==w-1) { r=1; In-=h; }
  if( h<4 || h%4>0 || (size_t(I)&15) || (size_t(Gx)&15) ) {
    for( y=0; y<h; y++ ) *Gx++=(*In++-*Ip++)*r;
  } else {
    _G=(__m128*) Gx; _Ip=(__m128*) Ip; _In=(__m128*) In; _r = SET(r);
    std::cout << "inside grad1's else, before the for, h = " << h << ", w = " << w << ", x = " << x << std::endl;
    //if i dont try to imshow the input image in the beggining, the program crashes here
    for(y=0; y<h; y+=4) *_G++=MUL(SUB(*_In++,*_Ip++),_r);
    std::cout << "inside grad1's else, after the for" << std::endl;
  }
  // compute column of Gy
  #define GRADY(r) *Gy++=(*In++-*Ip++)*r;
  Ip=I; In=Ip+1;
  // GRADY(1); Ip--; for(y=1; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
  y1=((~((size_t) Gy) + 1) & 15)/4; if(y1==0) y1=4; if(y1>h-1) y1=h-1;
  GRADY(1); Ip--; for(y=1; y<y1; y++) GRADY(.5f);
  _r = SET(.5f); _G=(__m128*) Gy;
  for(; y+4<h-1; y+=4, Ip+=4, In+=4, Gy+=4)
    *_G++=MUL(SUB(LDu(*In),LDu(*Ip)),_r);
  for(; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
  #undef GRADY
}

// build lookup table a[] s.t. a[x*n]~=acos(x) for x in [-1,1]
float* GradientMagnitudeChannel::acosTable() {
  const int n=10000, b=10; int i;
  static float a[n*2+b*2]; static bool init=false;
  float *a1=a+n+b; if( init ) return a1;
  for( i=-n-b; i<-n; i++ )   a1[i]=PI;
  for( i=-n; i<n; i++ )      a1[i]=float(acos(i/float(n)));
  for( i=n; i<n+b; i++ )     a1[i]=0;
  for( i=-n-b; i<n/10; i++ ) if( a1[i] > PI-1e-6f ) a1[i]=PI-1e-6f;
  init=true; return a1;
}
