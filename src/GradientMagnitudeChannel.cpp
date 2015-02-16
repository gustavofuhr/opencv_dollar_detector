#include "GradientMagnitudeChannel.h"

void GradientMagnitudeChannel::readGradientMagnitude(cv::FileNode gradMagNode)
{
	enabled = gradMagNode["enabled"];
	colorChannelIndex = gradMagNode["colorChn"];
	normalizationRadius = gradMagNode["normRad"];
	normalizationConstant = gradMagNode["normConst"];
	full = gradMagNode["full"];
  nChannels = gradMagNode["nChns"];
}


// compute x and y gradients for just one column (uses sse)
void grad1( float *I, float *Gx, float *Gy, int h, int w, int x ) {
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
float* acosTable() {
  const int n=10000, b=10; int i;
  static float a[n*2+b*2]; static bool init=false;
  float *a1=a+n+b; if( init ) return a1;
  for( i=-n-b; i<-n; i++ )   a1[i]=PI;
  for( i=-n; i<n; i++ )      a1[i]=float(acos(i/float(n)));
  for( i=n; i<n+b; i++ )     a1[i]=0;
  for( i=-n-b; i<n/10; i++ ) if( a1[i] > PI-1e-6f ) a1[i]=PI-1e-6f;
  init=true; return a1;
}

// compute gradient magnitude and orientation at each location (uses sse)
void gradMag( float *I, float *M, float *O, int h, int w, int d, bool full ) {
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
    for(c=0; c<d; c++) {
      grad1( I+x*h+c*w*h, Gx+c*h4, Gy+c*h4, h, w, x );
      for( y=0; y<h4/4; y++ ) {
        y1=h4/4*c+y;
        _M2[y1]=ADD(MUL(_Gx[y1],_Gx[y1]),MUL(_Gy[y1],_Gy[y1]));
        if( c==0 ) continue; _m = CMPGT( _M2[y1], _M2[y] );
        _M2[y] = OR( AND(_m,_M2[y1]), ANDNOT(_m,_M2[y]) );
        _Gx[y] = OR( AND(_m,_Gx[y1]), ANDNOT(_m,_Gx[y]) );
        _Gy[y] = OR( AND(_m,_Gy[y1]), ANDNOT(_m,_Gy[y]) );
      }
    }
    // compute gradient mangitude (M) and normalize Gx
    for( y=0; y<h4/4; y++ ) {
      _m = SSEMIN( RCPSQRT(_M2[y]), SET(1e10f) );
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
  alFree(Gx); alFree(Gy); alFree(M2);
}

// gradMagNorm( M, S, norm ) - operates on M - see gradientMag.m
// gradientMex('gradientMagNorm',M,S,normConst);
// normalize gradient magnitude at each location (uses sse)
void GradientMagnitudeChannel::gradMagNorm(float *M, float *S, int h, int w) {
  float norm = normalizationConstant;
  __m128 *_M, *_S, _norm; int i=0, n=h*w, n4=n/4;
  _S = (__m128*) S; _M = (__m128*) M; _norm = SET(norm);
  bool sse = !(size_t(M)&15) && !(size_t(S)&15);
  if(sse) { for(; i<n4; i++) *_M++=MUL(*_M,RCP(ADD(*_S++,_norm))); i*=4; }
  for(; i<n; i++) M[i] /= (S[i] + norm);
}

// M=gradientMex('gradientMag',I,channel,full); ou [M,O] = gradientMex('gradientMag',I,channel,full);
// if(!strcmp(action,"gradientMag")) mGradMag(nl,pl,nr,pr);
std::vector<float*> GradientMagnitudeChannel::mGradMag(float* source, int rows, int cols, int channel)
{
  /*
  int h, w, d, c, full; 
  float *I, *M, *O=0;
  checkArgs(nl,pl,nr,pr,1,2,3,3,&h,&w,&d,mxSINGLE_CLASS,(void**)&I);
  */

  /*
  void checkArgs( int nl, mxArray *pl[], int nr, const mxArray *pr[], int nl0, int nl1, int nr0, int nr1, int *h, int *w, int *d, mxClassID id, void **I )
  {
    const int *dims; int nDims;
    if( nl<nl0 || nl>nl1 ) mexErrMsghTxt("Incorrect number of outputs.");
    if( nr<nr0 || nr>nr1 ) mexErrMsgTxt("Incorrect number of inputs.");
    nDims = mxGetNumberOfDimensions(pr[0]); 
    dims = mxGetDimensions(pr[0]);
    *h=dims[0]; *w=dims[1]; 
    *d=(nDims==2) ? 1 : dims[2]; 
    *I = mxGetPr(pr[0]); // I recebe o valor de pr[0] (que é o I dentro do Matlab)
    if( nDims!=2 && nDims!=3 ) mexErrMsgTxt("I must be a 2D or 3D array.");
    if( mxGetClassID(pr[0])!=id ) mexErrMsgTxt("I has incorrect type.");
  }
  */

  // h = I.rows, w = I.cols
  // nDims =  mxGetNumberOfDimensions(pr[0]);
  // dims = mxGetDimensions(pr[0]);
  // if (nDims == 2) d = 1; else d = dims[2];
  // c = pr[1];


  int h = rows, w = cols, c = channel, d = 3; 
  float *M, *O=0;
  float *I = source;
  std::vector<float*> result(2);
  
  /*
  if(h<2 || w<2) mexErrMsgTxt("I must be at least 2x2.");
  c = (int) mxGetScalar(pr[1]); 
  full = (int) mxGetScalar(pr[2]);
  */
	if (h>=2 && w>=2)
	{
    // if( c>0 && c<=d ) { I += h*w*(c-1); d=1; }
		if (c>0 && c<=d)
		{
      // does this modify the input?
			I += h*w*(c-1); 
      d=1;
		}

    // pl[0] = mxCreateMatrix3(h,w,1,mxSINGLE_CLASS,0,(void**)&M);
    M = (float*)malloc(h*w*1*sizeof(float));

    // if(nl>=2) pl[1] = mxCreateMatrix3(h,w,1,mxSINGLE_CLASS,0,(void**)&O);
    O = (float*)malloc(h*w*1*sizeof(float));

    // gradMag(I, M, O, h, w, d, full>0 );
		// call to the actual function: gradMag(I, M, O, h, w, d, full>0 );
		// void gradMag(float*, float*, float*, int, int, int, bool);
    gradMag(I, M, O, h, w, d, full>0);

		// next, we assign the values of M and O to the matrix thats going to be returned
    result[0] = M;
    result[1] = O;
	}
	else
	{
		//error: matrix I must be at least 2x2
    std::cout << " # mGradMag error: provided matrix should have two dimensions!" << std::endl;
	}

	return result;
}