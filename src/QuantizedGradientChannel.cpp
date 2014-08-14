#include "QuantizedGradientChannel.h"

void QuantizedGradientChannel::readGradientHistogram(cv::FileNode histNode)
{
	enabled = histNode["enabled"];
	orientationChannels = histNode["nOrients"];
	useSoftBinning = histNode["softBin"];
	useHogNormalization = histNode["useHog"];
	clipHog = histNode["clipHog"];
  binSize = 8;
}

//Compute oriented gradient histograms
//H=gradHist(M,O,[...]) - see gradientHist.m
//H=gradientHist(M,O,binSize,p.nOrients,p.softBin,p.useHog,p.clipHog,full);
cv::Mat QuantizedGradientChannel::mGradHist(cv::Mat gradMag, cv::Mat gradOri, int full)
{
	cv::Mat result;
	float *M, *O, *H;
	
	//for now, we will make this trasnformations to avoid having to rewrite the other procedures
	M = cvMat2floatArray(gradMag);
	O = cvMat2floatArray(gradOri);

	//checkArgs procedure is called but it is not actually needed
	//probably just need to test some of the parameters, if that
	//checkArgs(nl,pl,nr,pr,1,3,2,8,&h,&w,&d,mxSINGLE_CLASS,(void**)&M);

	int h = gradMag.rows;
	int w = gradMag.cols;

	//next there's a bunch of tests to see which parameters were given and which are to be default, i wont put that here for now at least

  std::cout << "inside mGradHist, before if" << std::endl;

	if (orientationChannels == 0)
		result.data = NULL;
	else
	{
		if (useHogNormalization == 0)
    {
      std::cout << "inside mGradHist, before calling gradHist" << std::endl;
			gradHist(M, O, H, h, w, binSize, orientationChannels, useSoftBinning, full);
    }
		else
		{
			if (useHogNormalization == 1)
				hog(M, O, H, h, w, binSize, orientationChannels, useSoftBinning, full, clipHog );
			else
				fhog(M, O, H, h, w, binSize, orientationChannels, useSoftBinning, clipHog );
		}
		//the resulting histogram matrix is our return value
		result = floatArray2cvMat(H, h, w, CV_32FC3);
	}
	return result;
}

// compute nOrients gradient histograms per bin x bin block of pixels
void QuantizedGradientChannel::gradHist( float *M, float *O, float *H, int h, int w,
  int bin, int nOrients, int softBin, int full)
{
  //bin is zero right now, which causes a crash
  std::cout << "inside gradHist, bin = " << bin << std::endl;
  const int hb=h/bin, wb=w/bin, h0=hb*bin, w0=wb*bin, nb=wb*hb;
  const float s=(float)bin, sInv=1/s, sInv2=1/s/s;
  float *H0, *H1, *M0, *M1; int x, y; int *O0, *O1;
  O0=(int*)malloc(h*sizeof(int)); M0=(float*) malloc(h*sizeof(float));
  O1=(int*)malloc(h*sizeof(int)); M1=(float*) malloc(h*sizeof(float));
  // main loop
  for( x=0; x<w0; x++ ) {
    // compute target orientation bins for entire column - very fast
    gradQuantize(O+x*h,M+x*h,O0,O1,M0,M1,nb,h0,sInv2,nOrients,full,softBin>=0);

    if( softBin<0 && softBin%2==0 ) {
      // no interpolation w.r.t. either orienation or spatial bin
      H1=H+(x/bin)*hb;
      #define GH H1[O0[y]]+=M0[y]; y++;
      if( bin==1 )      for(y=0; y<h0;) { GH; H1++; }
      else if( bin==2 ) for(y=0; y<h0;) { GH; GH; H1++; }
      else if( bin==3 ) for(y=0; y<h0;) { GH; GH; GH; H1++; }
      else if( bin==4 ) for(y=0; y<h0;) { GH; GH; GH; GH; H1++; }
      else for( y=0; y<h0;) { for( int y1=0; y1<bin; y1++ ) { GH; } H1++; }
      #undef GH

    } else if( softBin%2==0 || bin==1 ) {
      // interpolate w.r.t. orientation only, not spatial bin
      H1=H+(x/bin)*hb;
      #define GH H1[O0[y]]+=M0[y]; H1[O1[y]]+=M1[y]; y++;
      if( bin==1 )      for(y=0; y<h0;) { GH; H1++; }
      else if( bin==2 ) for(y=0; y<h0;) { GH; GH; H1++; }
      else if( bin==3 ) for(y=0; y<h0;) { GH; GH; GH; H1++; }
      else if( bin==4 ) for(y=0; y<h0;) { GH; GH; GH; GH; H1++; }
      else for( y=0; y<h0;) { for( int y1=0; y1<bin; y1++ ) { GH; } H1++; }
      #undef GH

    } else {
      // interpolate using trilinear interpolation
      float ms[4], xyd, xb, yb, xd, yd, init; __m128 _m, _m0, _m1;
      bool hasLf, hasRt; int xb0, yb0;
      if( x==0 ) { init=(0+.5f)*sInv-0.5f; xb=init; }
      hasLf = xb>=0; xb0 = hasLf?(int)xb:-1; hasRt = xb0 < wb-1;
      xd=xb-xb0; xb+=sInv; yb=init; y=0;
      // macros for code conciseness
      #define GHinit yd=yb-yb0; yb+=sInv; H0=H+xb0*hb+yb0; xyd=xd*yd; \
        ms[0]=1-xd-yd+xyd; ms[1]=yd-xyd; ms[2]=xd-xyd; ms[3]=xyd;
      #define GH(H,ma,mb) H1=H; STRu(*H1,ADD(LDu(*H1),MUL(ma,mb)));
      // leading rows, no top bin
      for( ; y<bin/2; y++ ) {
        yb0=-1; GHinit;
        if(hasLf) { H0[O0[y]+1]+=ms[1]*M0[y]; H0[O1[y]+1]+=ms[1]*M1[y]; }
        if(hasRt) { H0[O0[y]+hb+1]+=ms[3]*M0[y]; H0[O1[y]+hb+1]+=ms[3]*M1[y]; }
      }
      // main rows, has top and bottom bins, use SSE for minor speedup
      if( softBin<0 ) for( ; ; y++ ) {
        yb0 = (int) yb; if(yb0>=hb-1) break; GHinit; _m0=SET(M0[y]);
        if(hasLf) { _m=SET(0,0,ms[1],ms[0]); GH(H0+O0[y],_m,_m0); }
        if(hasRt) { _m=SET(0,0,ms[3],ms[2]); GH(H0+O0[y]+hb,_m,_m0); }
      } else for( ; ; y++ ) {
        yb0 = (int) yb; if(yb0>=hb-1) break; GHinit;
        _m0=SET(M0[y]); _m1=SET(M1[y]);
        if(hasLf) { _m=SET(0,0,ms[1],ms[0]);
          GH(H0+O0[y],_m,_m0); GH(H0+O1[y],_m,_m1); }
        if(hasRt) { _m=SET(0,0,ms[3],ms[2]);
          GH(H0+O0[y]+hb,_m,_m0); GH(H0+O1[y]+hb,_m,_m1); }
      }
      // final rows, no bottom bin
      for( ; y<h0; y++ ) {
        yb0 = (int) yb; GHinit;
        if(hasLf) { H0[O0[y]]+=ms[0]*M0[y]; H0[O1[y]]+=ms[0]*M1[y]; }
        if(hasRt) { H0[O0[y]+hb]+=ms[2]*M0[y]; H0[O1[y]+hb]+=ms[2]*M1[y]; }
      }
      #undef GHinit
      #undef GH
    }
  }
  free(O0); free(O1); free(M0); free(M1);
  // normalize boundary bins which only get 7/8 of weight of interior bins
  if( softBin%2!=0 ) for( int o=0; o<nOrients; o++ ) {
    x=0; for( y=0; y<hb; y++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    y=0; for( x=0; x<wb; x++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    x=wb-1; for( y=0; y<hb; y++ ) H[o*nb+x*hb+y]*=8.f/7.f;
    y=hb-1; for( x=0; x<wb; x++ ) H[o*nb+x*hb+y]*=8.f/7.f;
  }
}

// helper for gradHist, quantize O and M into O0, O1 and M0, M1 (uses sse)
void QuantizedGradientChannel::gradQuantize( float *O, float *M, int *O0, int *O1, float *M0, float *M1, int nb, int n, float norm, int nOrients, int full, bool interpolate )
{
  // assumes all *OUTPUT* matrices are 4-byte aligned
  int i, o0, o1; float o, od, m;
  __m128i _o0, _o1, *_O0, *_O1; __m128 _o, _od, _m, *_M0, *_M1;
  // define useful constants
  const float oMult=(float)nOrients/(full?2*PI:PI); const int oMax=nOrients*nb;
  const __m128 _norm=SET(norm), _oMult=SET(oMult), _nbf=SET((float)nb);
  const __m128i _oMax=SET(oMax), _nb=SET(nb);
  // perform the majority of the work with sse
  _O0=(__m128i*) O0; _O1=(__m128i*) O1; _M0=(__m128*) M0; _M1=(__m128*) M1;
  if( interpolate ) for( i=0; i<=n-4; i+=4 ) {
    _o=MUL(LDu(O[i]),_oMult); _o0=CVT(_o); _od=SUB(_o,CVT(_o0));
    _o0=CVT(MUL(CVT(_o0),_nbf)); _o0=AND(CMPGT(_oMax,_o0),_o0); *_O0++=_o0;
    _o1=ADD(_o0,_nb); _o1=AND(CMPGT(_oMax,_o1),_o1); *_O1++=_o1;
    _m=MUL(LDu(M[i]),_norm); *_M1=MUL(_od,_m); *_M0++=SUB(_m,*_M1); _M1++;
  } else for( i=0; i<=n-4; i+=4 ) {
    _o=MUL(LDu(O[i]),_oMult); _o0=CVT(ADD(_o,SET(.5f)));
    _o0=CVT(MUL(CVT(_o0),_nbf)); _o0=AND(CMPGT(_oMax,_o0),_o0); *_O0++=_o0;
    *_M0++=MUL(LDu(M[i]),_norm); *_M1++=SET(0.f); *_O1++=SET(0);
  }
  // compute trailing locations without sse
  if( interpolate ) for( i; i<n; i++ ) {
    o=O[i]*oMult; o0=(int) o; od=o-o0;
    o0*=nb; if(o0>=oMax) o0=0; O0[i]=o0;
    o1=o0+nb; if(o1==oMax) o1=0; O1[i]=o1;
    m=M[i]*norm; M1[i]=od*m; M0[i]=m-M1[i];
  } else for( i; i<n; i++ ) {
    o=O[i]*oMult; o0=(int) (o+.5f);
    o0*=nb; if(o0>=oMax) o0=0; O0[i]=o0;
    M0[i]=M[i]*norm; M1[i]=0; O1[i]=0;
  }
}

// compute HOG features
void QuantizedGradientChannel::hog( float *M, float *O, float *H, int h, int w, int binSize, int nOrients, int softBin, int full, float clip )
{
  float *N, *R; const int hb=h/binSize, wb=w/binSize, nb=hb*wb;
  // compute unnormalized gradient histograms
  R = (float*) calloc(wb*hb*nOrients,sizeof(float));
  gradHist( M, O, R, h, w, binSize, nOrients, softBin, full );
  // compute block normalization values
  N = hogNormMatrix( R, nOrients, hb, wb, binSize );
  // perform four normalizations per spatial block
  hogChannels( H, R, N, hb, wb, nOrients, clip, 0 );
  free(N); free(R);
}

// compute FHOG features
void QuantizedGradientChannel::fhog( float *M, float *O, float *H, int h, int w, int binSize, int nOrients, int softBin, float clip )
{
  const int hb=h/binSize, wb=w/binSize, nb=hb*wb, nbo=nb*nOrients;
  float *N, *R1, *R2; int o, x;
  // compute unnormalized constrast sensitive histograms
  R1 = (float*) calloc(wb*hb*nOrients*2,sizeof(float));
  gradHist( M, O, R1, h, w, binSize, nOrients*2, softBin, true );
  // compute unnormalized contrast insensitive histograms
  R2 = (float*) calloc(wb*hb*nOrients,sizeof(float));
  for( o=0; o<nOrients; o++ ) for( x=0; x<nb; x++ )
    R2[o*nb+x] = R1[o*nb+x]+R1[(o+nOrients)*nb+x];
  // compute block normalization values
  N = hogNormMatrix( R2, nOrients, hb, wb, binSize );
  // normalized histograms and texture channels
  hogChannels( H+nbo*0, R1, N, hb, wb, nOrients*2, clip, 1 );
  hogChannels( H+nbo*2, R2, N, hb, wb, nOrients*1, clip, 1 );
  hogChannels( H+nbo*3, R1, N, hb, wb, nOrients*2, clip, 2 );
  free(N); free(R1); free(R2);
}

// HOG helper: compute 2x2 block normalization values (padded by 1 pixel)
float* QuantizedGradientChannel::hogNormMatrix( float *H, int nOrients, int hb, int wb, int bin ) {
  float *N, *N1, *n; int o, x, y, dx, dy, hb1=hb+1, wb1=wb+1;
  float eps = 1e-4f/4/bin/bin/bin/bin; // precise backward equality
  N = (float*) calloc(hb1*wb1,sizeof(float)); N1=N+hb1+1;
  for( o=0; o<nOrients; o++ ) for( x=0; x<wb; x++ ) for( y=0; y<hb; y++ )
    N1[x*hb1+y] += H[o*wb*hb+x*hb+y]*H[o*wb*hb+x*hb+y];
  for( x=0; x<wb-1; x++ ) for( y=0; y<hb-1; y++ ) {
    n=N1+x*hb1+y; *n=1/float(sqrt(n[0]+n[1]+n[hb1]+n[hb1+1]+eps)); }
  x=0;     dx= 1; dy= 1; y=0;                  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=0;     dx= 1; dy= 0; for(y=0; y<hb1; y++)  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=0;     dx= 1; dy=-1; y=hb1-1;              N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=wb1-1; dx=-1; dy= 1; y=0;                  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=wb1-1; dx=-1; dy= 0; for( y=0; y<hb1; y++) N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  x=wb1-1; dx=-1; dy=-1; y=hb1-1;              N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  y=0;     dx= 0; dy= 1; for(x=0; x<wb1; x++)  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  y=hb1-1; dx= 0; dy=-1; for(x=0; x<wb1; x++)  N[x*hb1+y]=N[(x+dx)*hb1+y+dy];
  return N;
}

// HOG helper: compute HOG or FHOG channels
void QuantizedGradientChannel::hogChannels( float *H, const float *R, const float *N,
  int hb, int wb, int nOrients, float clip, int type )
{
  #define GETT(blk) t=R1[y]*N1[y-(blk)]; if(t>clip) t=clip; c++;
  const float r=.2357f; int o, x, y, c; float t;
  const int nb=wb*hb, nbo=nOrients*nb, hb1=hb+1;
  for( o=0; o<nOrients; o++ ) for( x=0; x<wb; x++ ) {
    const float *R1=R+o*nb+x*hb, *N1=N+x*hb1+hb1+1;
    float *H1 = (type<=1) ? (H+o*nb+x*hb) : (H+x*hb);
    if( type==0) for( y=0; y<hb; y++ ) {
      // store each orientation and normalization (nOrients*4 channels)
      c=-1; GETT(0); H1[c*nbo+y]=t; GETT(1); H1[c*nbo+y]=t;
      GETT(hb1); H1[c*nbo+y]=t; GETT(hb1+1); H1[c*nbo+y]=t;
    } else if( type==1 ) for( y=0; y<hb; y++ ) {
      // sum across all normalizations (nOrients channels)
      c=-1; GETT(0); H1[y]+=t*.5f; GETT(1); H1[y]+=t*.5f;
      GETT(hb1); H1[y]+=t*.5f; GETT(hb1+1); H1[y]+=t*.5f;
    } else if( type==2 ) for( y=0; y<hb; y++ ) {
      // sum across all orientations (4 channels)
      c=-1; GETT(0); H1[c*nb+y]+=t*r; GETT(1); H1[c*nb+y]+=t*r;
      GETT(hb1); H1[c*nb+y]+=t*r; GETT(hb1+1); H1[c*nb+y]+=t*r;
    }
  }
  #undef GETT
}


