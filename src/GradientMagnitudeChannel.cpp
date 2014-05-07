#include "GradientMagnitudeChannel.h"

void GradientMagnitudeChannel::readGradientMagnitude(cv::FileNode gradMagNode)
{
	enabled = gradMagNode["enabled"];
	colorChannelIndex = gradMagNode["colorChn"];
	normalizationRadius = gradMagNode["normRad"];
	normalizationConstant = gradMagNode["normConst"];
	full = gradMagNode["full"];
}

// [M,O] = gradMag( I, channel, full ) - see gradientMag.m
cv::Mat* GradientMagnitudeChannel::mGradMag(cv::Mat I, int channel)
{
	int c, d; 
	float *M, *O=0;
	cv::Mat resultMatrix[2];

	//checkArgs procedure is called but it is not actually needed
	//probably just need to test some of the parameters, if that
	//checkArgs(nl,pl,nr,pr,1,2,3,3,&h,&w,&d,mxSINGLE_CLASS,(void**)&I);

	if (I.rows>=2 && I.cols>=2)
	{
		//for now, the gradMag procedure will be unchanged so we put the raw data of the image in a float pointer to be used there 
		float *If = (float*)I.data;
		if (channel>0 && channel<=I.dims)
		{
			If += I.rows*I.cols*(channel-1); 
			d=1;  
		}
		//i dont think this next if statement will be necessary
		/*if (nl>=2)
			Mat pl; //pl[0] = mxCreateMatrix3(h,w,1,mxSINGLE_CLASS,0,(void**)&M);*/

		//call to the actual function: 
		//void gradMag(float*, float*, float*, int, int, int, bool);
		gradMag(If, M, O, I.rows, I.cols, I.dims, full>0 );

		//next, we assign the values of M and O to the matrix thats going to be returned
		resultMatrix[0].data = (uchar*)M;
		resultMatrix[1].data = (uchar*)O;
	}
	else
	{
		//error: matrix I must be at least 2x2
		resultMatrix[0].data = NULL;
		resultMatrix[1].data = NULL;
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
void GradientMagnitudeChannel::gradMag( float *I, float *M, float *O, int h, int w, int d, bool full ) {
  int x, y, y1, c, h4, s; float *Gx, *Gy, *M2; __m128 *_Gx, *_Gy, *_M2, _m;
  float *acost = acosTable(), acMult=10000.0f;
  // allocate memory for storing one column of output (padded so h4%4==0)
  h4=(h%4==0) ? h : h-(h%4)+4; s=d*h4*sizeof(float);
  M2=(float*) malloc(s); _M2=(__m128*) M2;
  Gx=(float*) malloc(s); _Gx=(__m128*) Gx;
  Gy=(float*) malloc(s); _Gy=(__m128*) Gy;
  // compute gradient magnitude and orientation for each column
  for( x=0; x<w; x++ ) {
    // compute gradients (Gx, Gy) with maximum squared magnitude (M2)
    for(c=0; c<d; c++) {
      grad1( I+x*h+c*w*h, Gx+c*h4, Gy+c*h4, h, w, x );
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
  int y, y1; float *Ip, *In, r; __m128 *_Ip, *_In, *_G, _r;
  // compute column of Gx
  Ip=I-h; In=I+h; r=.5f;
  if(x==0) { r=1; Ip+=h; } else if(x==w-1) { r=1; In-=h; }
  if( h<4 || h%4>0 || (size_t(I)&15) || (size_t(Gx)&15) ) {
    for( y=0; y<h; y++ ) *Gx++=(*In++-*Ip++)*r;
  } else {
    _G=(__m128*) Gx; _Ip=(__m128*) Ip; _In=(__m128*) In; _r = SET(r);
    for(y=0; y<h; y+=4) *_G++=MUL(SUB(*_In++,*_Ip++),_r);
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
