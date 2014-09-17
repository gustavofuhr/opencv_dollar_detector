#include "mex.h"
#include <cmath>

// using namespace std;
typedef unsigned int uint32;

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  // get inputs
  float *chns = (float*) mxGetData(prhs[0]);
  mxArray *trees = (mxArray*) prhs[1];
  const int shrink = (int) mxGetScalar(prhs[2]);
  const int modelHt = (int) mxGetScalar(prhs[3]);
  const int modelWd = (int) mxGetScalar(prhs[4]);
  const int stride = (int) mxGetScalar(prhs[5]);
  const float cascThr = (float) mxGetScalar(prhs[6]);

  // extract relevant fields from trees
  float *thrs = (float*) mxGetData(mxGetField(trees,0,"thrs"));
  float *hs = (float*) mxGetData(mxGetField(trees,0,"hs"));
  uint32 *fids = (uint32*) mxGetData(mxGetField(trees,0,"fids"));
  uint32 *child = (uint32*) mxGetData(mxGetField(trees,0,"child"));
  const int treeDepth = mxGetField(trees,0,"treeDepth")==NULL ? 0 :
    (int) mxGetScalar(mxGetField(trees,0,"treeDepth"));

  // get dimensions and constants
  const mwSize *chnsSize = mxGetDimensions(prhs[0]);
  const int height = (int) chnsSize[0];
  const int width = (int) chnsSize[1];
  const int nChns = mxGetNumberOfDimensions(prhs[0])<=2 ? 1 : (int) chnsSize[2];
  const mwSize *fidsSize = mxGetDimensions(mxGetField(trees,0,"fids"));
  const int nTreeNodes = (int) fidsSize[0];
  const int nTrees = (int) fidsSize[1];
  const int height1 = (int) ceil(float(height*shrink-modelHt+1)/stride);
  const int width1 = (int) ceil(float(width*shrink-modelWd+1)/stride);

  // construct cids array
  int nFtrs = modelHt/shrink*modelWd/shrink*nChns;
  uint32 *cids = new uint32[nFtrs]; int m=0;
  for( int z=0; z<nChns; z++ )
    for( int c=0; c<modelWd/shrink; c++ )
      for( int r=0; r<modelHt/shrink; r++ )
        cids[m++] = z*width*height + c*height + r;


  mexPrintf("shrink=%d, modelHt=%d, modelWd=%d, stride=%d, cascThr=%f, treeDepth=%d\n", shrink, modelHt, modelWd, stride, cascThr, treeDepth);
  mexPrintf("height=%d, width=%d, nChns=%d, nTreeNodes=%d, nTrees=%d, height1=%d, width1=%d, nFtrs=%d\n\n", height, width, nChns, nTreeNodes, nTrees, height1, width1, nFtrs);

  mexPrintf("\nFirst ten elements of chns:\n");
  for (int i=0; i < 10; i++)
    mexPrintf(" %f", chns[i]);
  mexPrintf("\n");

  mexPrintf("\nFirst ten elements of thrs:\n");
  for (int i=0; i < 10; i++)
    mexPrintf(" %f", thrs[i]);
  mexPrintf("\n");

  mexPrintf("\nFirst ten elements of hs:\n");
  for (int i=0; i < 10; i++)
    mexPrintf(" %f", hs[i]);
  mexPrintf("\n");

  mexPrintf("\nFirst ten elements of fids:\n");
  for (int i=0; i < 10; i++)
    mexPrintf(" %d", fids[i]);
  mexPrintf("\n");

  mexPrintf("\nFirst ten elements of child:\n");
  for (int i=0; i < 10; i++)
    mexPrintf(" %d", child[i]);
  mexPrintf("\n");

  mexPrintf("\nFirst ten elements of cids:\n");
  for (int i=0; i < 10; i++)
    mexPrintf(" %d", cids[i]);
  mexPrintf("\n");

  /*
  shrink=4, modelHt=128, modelWd=64, stride=4, cascThr=-1.000000, treeDepth=2
  height=152, width=186, nChns=10, nTreeNodes=7, nTrees=2048, height1=121, width1=171, nFtrs=5120


  First ten elements of thrs:
   0.190968 0.273173 0.064298 0.000000 0.000000 0.000000 0.000000 0.063845 0.083150 0.303183

  First ten elements of hs:
   0.530309 -0.520309 0.530309 -0.520309 0.530309 -0.520309 0.530309 -0.439302 0.449302 -0.439302

  First ten elements of fids:
   2328 2392 3271 0 0 0 0 2301 2222 2697

  First ten elements of child:
   2 4 6 0 0 0 0 2 4 6
  */

  plhs[0] = mxCreateNumericMatrix(m,5,mxDOUBLE_CLASS,mxREAL);
}


