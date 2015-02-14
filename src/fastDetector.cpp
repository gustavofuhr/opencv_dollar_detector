BB_Array Detector::applyFastCalibratedDetectorToFrame(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, 
	float *hs, uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns, int imageWidth, int imageHeight, 
	cv::Mat_<float> &homography, cv::Mat_<float> &projection)
{
	BB_Array result;

	std::vector<float*> scales_chns(opts.pPyramid.computedScales, NULL);
	std::vector<uint32*> scales_cids(opts.pPyramid.computedScales, NULL);

	// pre-compute the way we access the features for each scale
	for (int i=0; i < opts.pPyramid.computedScales; i++) {
		int height = pyramid[i].image.rows;
		int width = pyramid[i].image.cols;

		int channels = opts.pPyramid.pChns.pColor.nChannels + opts.pPyramid.pChns.pGradMag.nChannels + opts.pPyramid.pChns.pGradHist.nChannels;
		float* chns = (float*)malloc(height*width*channels*sizeof(float));
		features2floatArray(pyramid[i], chns, height, width,  opts.pPyramid.pChns.pColor.nChannels, opts.pPyramid.pChns.pGradMag.nChannels, opts.pPyramid.pChns.pGradHist.nChannels);
		scales_chns[i] = chns;


		int nFtrs = modelHt/shrink*modelWd/shrink*nChns;
	  	uint32 *cids = new uint32[nFtrs]; int m=0;
	  	
	  	for( int z=0; z<nChns; z++ ) {
	    	for( int cc=0; cc<modelWd/shrink; cc++ ) {
	      		for( int rr=0; rr<modelHt/shrink; rr++ ) {
	        		cids[m++] = z*width*height + cc*height + rr;
	        	}
	        }
	    }

	    scales_cids[i] = cids;
	}

	std::vector<int> rs, cs; std::vector<float> hs1;
	for(int col=0; col<pyramid[0].image.cols-modelWd; col++) 
	{
		for(int row=modelHt; row<pyramid[0].image.rows; row++) 
		{
			int c = col;
			int r = 0;
			double  curBBworldHeight=findWorldHeight(col+(modelWd/2), row, r, projection, homography);

			while (curBBworldHeight > targetPedestrianHeight && r >= row-modelHt)
			{
				r++;
				curBBworldHeight = findWorldHeight(col+(modelWd/2), row, r, projection, homography);
			}

			int curBBimageHeight = row-r;
			float targetScale = opts.modelDs[0]/curBBimageHeight;

			int s=0;
			if (r != 0)
			{
				double curDiff = 1.0;
				double newDiff = opts.pPyramid.scales[s] - targetScale;

				while (newDiff < curDiff)
				{
					s++;
					curDiff = newDiff;
					newDiff = opts.pPyramid.scales[s] - targetScale;
				}
			}

			// at the end of this, i know the values of c, r and s(i)
			int height = pyramid[s].image.rows;	

			float h=0, *chns1=scales_chns[s]+(r*stride/shrink) + (c*stride/shrink)*height;	
			if( treeDepth==1 ) {
		      // specialized case for treeDepth==1
		      for( int t = 0; t < nTrees; t++ ) {
		        uint32 offset=t*nTreeNodes, k=offset, k0=0;
		        getChild(chns1,cids,fids,thrs,offset,k0,k);
		        h += hs[k]; if( h<=cascThr ) break;
		      }
		    } else if( treeDepth==2 ) {
		      // specialized case for treeDepth==2
		      for( int t = 0; t < nTrees; t++ ) {
		        uint32 offset=t*nTreeNodes, k=offset, k0=0;
		        getChild(chns1,cids,fids,thrs,offset,k0,k);
		        getChild(chns1,cids,fids,thrs,offset,k0,k);
		        h += hs[k]; if( h<=cascThr ) break;
		      }
		    } else if( treeDepth>2) {
		      // specialized case for treeDepth>2
		      for( int t = 0; t < nTrees; t++ ) {
		        uint32 offset=t*nTreeNodes, k=offset, k0=0;
		        for( int i=0; i<treeDepth; i++ )
		          getChild(chns1,cids,fids,thrs,offset,k0,k);
		        h += hs[k]; if( h<=cascThr ) break;
		      }
		    } else {
		      // general case (variable tree depth)
		      for( int t = 0; t < nTrees; t++ ) {
		        uint32 offset=t*nTreeNodes, k=offset, k0=k;
		        while( child[k] ) {
		          float ftr = chns1[cids[fids[k]]];
		          k = (ftr<thrs[k]) ? 1 : 0;
		          k0 = k = child[k0]-k+offset;
		        }
		        h += hs[k]; if( h<=cascThr ) break;
		      }
	    	}
	    	if(h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); } // detecção	

		}
	}

	for(int j=0; j<m; j++ )
	{
		BoundingBox bb = pyramidRowColumn2BoundingBox(rs[j], cs[j],  modelHt, modelWd, i, stride) ;

		// bb.topLeftPoint.x = cs[j]*stride;
		// bb.topLeftPoint.x = (bb.topLeftPoint.x+shift[1])/opts.pPyramid.scales_w[i];
		// bb.topLeftPoint.y = rs[j]*stride;
		// bb.topLeftPoint.y = (bb.topLeftPoint.y+shift[0])/opts.pPyramid.scales_h[i];
		// bb.height = opts.modelDs[0]/opts.pPyramid.scales[i];
		// bb.width = opts.modelDs[1]/opts.pPyramid.scales[i];
		bb.score = hs1[j];
		bb.scale = i;
		result.push_back(bb);
	}

	cs.clear();
	rs.clear();
	hs1.clear();

	return result;
}