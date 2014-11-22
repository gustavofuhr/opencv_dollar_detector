#include "Detector.h"

// i dont know if its gonna be needed but this is start
void Detector::exportDetectorModel(cv::String fileName)
{
	cv::FileStorage xml;
	
	xml.open(fileName, cv::FileStorage::WRITE);

	xml << "opts" << "{";
		xml << "pPyramid" << "{";
			xml << "pChns" << "{";
				xml << "shrink" << opts.pPyramid.pChns.shrink;
				xml << "pColor" << "{";
					xml << "enabled" << opts.pPyramid.pChns.pColor.enabled;
				xml << "}";
			xml << "}";
		xml << "}";
	xml << "stride" << opts.stride;
	xml << "}";
	
	//xml << "clf" << this->clf;

	xml.release();
}

// reads the detector model from the xml model
// for now, it must be like this since the current model was not written by this program 
// this will change after we are set on a class structure
void Detector::importDetectorModel(cv::String fileName)
{
	cv::FileStorage xml;

	xml.open(fileName, cv::FileStorage::READ);

	if (!xml.isOpened())
	{
		std::cerr << " # Failed to open " << fileName << std::endl;
	}
	else
	{
		opts.readOptions(xml["detector"]["opts"]);
	
		xml["detector"]["clf"]["fids"] >> fids;
		xml["detector"]["clf"]["child"] >> child;
		xml["detector"]["clf"]["thrs"] >> thrs;
		xml["detector"]["clf"]["hs"] >> hs;
		xml["detector"]["clf"]["weights"] >> weights;
		xml["detector"]["clf"]["depth"] >> depth;
		xml["detector"]["clf"]["errs"] >> errs;
		xml["detector"]["clf"]["losses"] >> losses;		
		xml["detector"]["clf"]["treeDepth"] >> treeDepth;	

		timeSpentInDetection = 0;

		xml.release();
	}
}

void showDetections(cv::Mat I, BB_Array detections, cv::String windowName)
{
	cv::Mat img = I.clone();
	for (int j = 0; j<detections.size(); j++) 
		detections[j].plot(img, cv::Scalar(0,255,0));

	cv::imshow(windowName, img);
}

void printDetections(BB_Array detections, int frameIndex)
{
	std::cout << "Detections in frame " << frameIndex << ":\n";
	for (int i=0; i < detections.size(); i++)
		std::cout << detections[i].toString(frameIndex) << std::endl;
	std::cout << std::endl;
}

// this procedure was just copied verbatim
inline void getChild(float *chns1, uint32 *cids, uint32 *fids, float *thrs, uint32 offset, uint32 &k0, uint32 &k)
{
  float ftr = chns1[cids[fids[k]]];
  k = (ftr<thrs[k]) ? 1 : 2;
  k0=k+=k0*2; k+=offset;
}

/*
// this is not done yet!
BB_Array Detector::applyCalibratedDetectorToFrame(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, 
												  float *hs, uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns, cv::Mat homography)
{
	BB_Array result;

	int lastRow = (int)ceil(float(pyramid[0].image.rows*shrink-modelHt+1)/stride);
	int lastCol = (int)ceil(float(pyramid[0].image.cols*shrink-modelWd+1)/stride);

	int channels = opts.pPyramid.pChns.pColor.nChannels + opts.pPyramid.pChns.pGradMag.nChannels + opts.pPyramid.pChns.pGradHist.nChannels;

	std::vector<int> rs, cs; std::vector<float> hs1;
	for( int c=0; c<lastCol; c++ ) 
	{
		for( int r=0; r<lastRow; r++ ) 
		{
			cv::Point groundPoint = imagePoint2worldPoint(c, r+modelHt, 1, homography);
			cv::Point topPoint = imagePoint2worldPoint(c, r, 1, homography);

			//float boundingBoxWorldHeight = difference in height between topPoint and groundPoint;

			if (boundingBoxWorldHeight >= minPedestrianHeight && boundingBoxWorldHeight <= maxPedestrianHeight)
			{
				// discover which of the scales is the correct one
				int scaleIndex = findBestScale(boundingBoxWorldHeight, std::vector<Info> pyramid);
				float* chns = (float*)malloc(height*width*channels*sizeof(float));
				chns = features2floatArray(pyramid[scaleIndex], chns, height, width,  opts.pPyramid.pChns.pColor.nChannels, opts.pPyramid.pChns.pGradMag.nChannels, opts.pPyramid.pChns.pGradHist.nChannels);

				// the rest of the detection works the same way:

				// construct cids array
			  	int nFtrs = modelHt/shrink*modelWd/shrink*nChns;
			  	uint32 *cids = new uint32[nFtrs]; int m=0;
			  	for( int z=0; z<nChns; z++ )
			    	for( int c=0; c<modelWd/shrink; c++ )
			      		for( int r=0; r<modelHt/shrink; r++ )
			        		cids[m++] = z*width*height + c*height + r;

			    float h=0, *chns1=chns+(r*stride/shrink) + (c*stride/shrink)*height;
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
				if(h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }

				delete [] cids;
				free(chns);
				m=cs.size();

				// shift=(modelDsPad-modelDs)/2-pad;
				double shift[2];
				shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
				shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

				for(int j=0; j<m; j++ )
				{
					BoundingBox bb;
					bb.topLeftPoint.x = cs[j]*stride;
					bb.topLeftPoint.x = (bb.topLeftPoint.x+shift[1])/opts.pPyramid.scales_w[i];
					bb.topLeftPoint.y = rs[j]*stride;
					bb.topLeftPoint.y = (bb.topLeftPoint.y+shift[0])/opts.pPyramid.scales_h[i];
					bb.height = opts.modelDs[0]/opts.pPyramid.scales[i];
					bb.width = opts.modelDs[1]/opts.pPyramid.scales[i];
					bb.score = hs1[j];
					bb.scale = i;
					result.push_back(bb);
				}

				cs.clear();
				rs.clear();
				hs1.clear();
			}
		}	
	}
}
*/

BB_Array Detector::applyDetectorToFrame(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, float *hs, 
										uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns)
{
	BB_Array result;


	printf("Number of scales: %d\n", opts.pPyramid.computedScales);
	// this became a simple loop because we will apply just one detector here, 
	// to apply multiple detector models you need to create multiple Detector objects. 
	for (int i = 0; i < opts.pPyramid.computedScales; i++)
	{
		// in the original file: *chnsSize = mxGetDimensions(P.data{i});
		// const int height = (int) chnsSize[0];
  		// const int width = (int) chnsSize[1];
  		// const int nChns = mxGetNumberOfDimensions(prhs[0])<=2 ? 1 : (int) chnsSize[2];
		int height = pyramid[i].image.rows;
		int width = pyramid[i].image.cols;
		fprintf("Size of the downsampled image: %d x %d\n", width, height);

		int height1 = (int)ceil(float(height*shrink-modelHt+1)/stride);
		int width1 = (int)ceil(float(width*shrink-modelWd+1)/stride);
		fprintf("Size of the patche: %d x %d\n", width, height);

		int channels = opts.pPyramid.pChns.pColor.nChannels + opts.pPyramid.pChns.pGradMag.nChannels + opts.pPyramid.pChns.pGradHist.nChannels;
		float* chns = (float*)malloc(height*width*channels*sizeof(float));
		features2floatArray(pyramid[i], chns, height, width,  opts.pPyramid.pChns.pColor.nChannels, opts.pPyramid.pChns.pGradMag.nChannels, opts.pPyramid.pChns.pGradHist.nChannels);
		
		/*
		// debug: read chns from file
	  	cv::Mat scalei;
	  	std::string scaleName;

	  	scaleName += "scale";
	  	if (i < 9)
	  		scaleName += "0";
	  	std::ostringstream scaleNumber;
        scaleNumber << (i+1);
	  	scaleName += scaleNumber.str();

		xml["pyramid"][scaleName] >> scalei;
		//float* floatScale = cvImage2floatArray(scalei, 1);
		//printElements(floatScale, scalei.rows, scaleName + " read from xml file");
		free(chns);
		chns = (float*) malloc((ch1Size+ch2Size+ch3Size)*sizeof(float));
		chns = cvImage2floatArray(scalei, 1);
		// debug */

		// construct cids array
	  	int nFtrs = modelHt/shrink*modelWd/shrink*nChns;
	  	uint32 *cids = new uint32[nFtrs]; int m=0;
	  	for( int z=0; z<nChns; z++ )
	    	for( int c=0; c<modelWd/shrink; c++ )
	      		for( int r=0; r<modelHt/shrink; r++ )
	        		cids[m++] = z*width*height + c*height + r;

		/*
		// debug: prints values of several variables, all of these return correct results
		// shrink=4, modelHt=128, modelWd=64, stride=4, cascThr=-1.000000, treeDepth=2
		// height=152, width=186, nChns=10, nTreeNodes=7, nTrees=2048, height1=121, width1=171, nFtrs=5120
		std::cout << "shrink=" << shrink << ", modelHt=" << modelHt << ", modelWd=" << modelWd << ", stride=" << stride << ", cascThr=" << cascThr << ", treeDepth=" << treeDepth <<  ", modelDs=(" <<
			opts.modelDs[0] << "," << opts.modelDs[1] << ")" << std::endl;
		std::cout << "height=" << height << ", width=" << width << ", nChns=" << nChns <<  ", nTreeNodes=" << nTreeNodes << ", nTrees=" << nTrees << ", height1=" << height1 << 
			", width1=" << width1 << ", nFtrs=" << nFtrs << std::endl;
		// debug */

		// apply classifier to each patch
  		std::vector<int> rs, cs; std::vector<float> hs1;
  		for( int c=0; c<width1; c++ ) 
  		{
  			for( int r=0; r<height1; r++ ) 
  			{
			    float h=0, *chns1=chns+(r*stride/shrink) + (c*stride/shrink)*height;
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
		    if(h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }
		  }
		}
		delete [] cids;
		free(chns);
		m=cs.size();

		// shift=(modelDsPad-modelDs)/2-pad;
		double shift[2];
		shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
		shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

		for(int j=0; j<m; j++ )
		{
			BoundingBox bb;
			bb.topLeftPoint.x = cs[j]*stride;
			bb.topLeftPoint.x = (bb.topLeftPoint.x+shift[1])/opts.pPyramid.scales_w[i];
			bb.topLeftPoint.y = rs[j]*stride;
			bb.topLeftPoint.y = (bb.topLeftPoint.y+shift[0])/opts.pPyramid.scales_h[i];
			bb.height = opts.modelDs[0]/opts.pPyramid.scales[i];
			bb.width = opts.modelDs[1]/opts.pPyramid.scales[i];
			bb.score = hs1[j];
			bb.scale = i;
			result.push_back(bb);
		}

		cs.clear();
		rs.clear();
		hs1.clear();
	}

	return result;
}



BB_Array Detector::applyDetectorToFrameSmart(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, float *hs, 
										uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns)
{
	BB_Array result;

	// this became a simple loop because we will apply just one detector here, 
	// to apply multiple detector models you need to create multiple Detector objects. 
	for (int i = 0; i < opts.pPyramid.computedScales; i++)
	{
		// in the original file: *chnsSize = mxGetDimensions(P.data{i});
		// const int height = (int) chnsSize[0];
  		// const int width = (int) chnsSize[1];
  		// const int nChns = mxGetNumberOfDimensions(prhs[0])<=2 ? 1 : (int) chnsSize[2];
		int height = pyramid[i].image.rows;
		int width = pyramid[i].image.cols;

		int height1 = (int)ceil(float(height*shrink-modelHt+1)/stride);
		int width1 = (int)ceil(float(width*shrink-modelWd+1)/stride);

		int channels = opts.pPyramid.pChns.pColor.nChannels + opts.pPyramid.pChns.pGradMag.nChannels + opts.pPyramid.pChns.pGradHist.nChannels;
		float* chns = (float*)malloc(height*width*channels*sizeof(float));
		features2floatArray(pyramid[i], chns, height, width,  opts.pPyramid.pChns.pColor.nChannels, opts.pPyramid.pChns.pGradMag.nChannels, opts.pPyramid.pChns.pGradHist.nChannels);
		
		/*
		// debug: read chns from file
	  	cv::Mat scalei;
	  	std::string scaleName;

	  	scaleName += "scale";
	  	if (i < 9)
	  		scaleName += "0";
	  	std::ostringstream scaleNumber;
        scaleNumber << (i+1);
	  	scaleName += scaleNumber.str();

		xml["pyramid"][scaleName] >> scalei;
		//float* floatScale = cvImage2floatArray(scalei, 1);
		//printElements(floatScale, scalei.rows, scaleName + " read from xml file");
		free(chns);
		chns = (float*) malloc((ch1Size+ch2Size+ch3Size)*sizeof(float));
		chns = cvImage2floatArray(scalei, 1);
		// debug */

		// construct cids array
	  	int nFtrs = modelHt/shrink*modelWd/shrink*nChns;
	  	uint32 *cids = new uint32[nFtrs]; int m=0;
	  	for( int z=0; z<nChns; z++ )
	    	for( int c=0; c<modelWd/shrink; c++ )
	      		for( int r=0; r<modelHt/shrink; r++ )
	        		cids[m++] = z*width*height + c*height + r;

		/*
		// debug: prints values of several variables, all of these return correct results
		// shrink=4, modelHt=128, modelWd=64, stride=4, cascThr=-1.000000, treeDepth=2
		// height=152, width=186, nChns=10, nTreeNodes=7, nTrees=2048, height1=121, width1=171, nFtrs=5120
		std::cout << "shrink=" << shrink << ", modelHt=" << modelHt << ", modelWd=" << modelWd << ", stride=" << stride << ", cascThr=" << cascThr << ", treeDepth=" << treeDepth <<  ", modelDs=(" <<
			opts.modelDs[0] << "," << opts.modelDs[1] << ")" << std::endl;
		std::cout << "height=" << height << ", width=" << width << ", nChns=" << nChns <<  ", nTreeNodes=" << nTreeNodes << ", nTrees=" << nTrees << ", height1=" << height1 << 
			", width1=" << width1 << ", nFtrs=" << nFtrs << std::endl;
		// debug */

		// apply classifier to each patch
  		std::vector<int> rs, cs; std::vector<float> hs1;
  		for( int c=0; c<width1; c++ ) 
  		{
  			for( int r=0; r<height1; r++ ) 
  			{
			    float h=0, *chns1=chns+(r*stride/shrink) + (c*stride/shrink)*height;
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
		    	if(h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }
		  	}
		}
		delete [] cids;
		free(chns);
		m=cs.size();

		// shift=(modelDsPad-modelDs)/2-pad;
		double shift[2];
		shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
		shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

		for(int j=0; j<m; j++ )
		{
			BoundingBox bb;
			bb.topLeftPoint.x = cs[j]*stride;
			bb.topLeftPoint.x = (bb.topLeftPoint.x+shift[1])/opts.pPyramid.scales_w[i];
			bb.topLeftPoint.y = rs[j]*stride;
			bb.topLeftPoint.y = (bb.topLeftPoint.y+shift[0])/opts.pPyramid.scales_h[i];
			bb.height = opts.modelDs[0]/opts.pPyramid.scales[i];
			bb.width = opts.modelDs[1]/opts.pPyramid.scales[i];
			bb.score = hs1[j];
			bb.scale = i;
			result.push_back(bb);
		}

		cs.clear();
		rs.clear();
		hs1.clear();
	}

	return result;
}




//bb = acfDetect1(P.data{i},Ds{j}.clf,shrink,modelDsPad(1),modelDsPad(2),opts.stride,opts.cascThr);
void Detector::acfDetect(std::vector<std::string> imageNames, std::string dataSetDirectoryName, int firstFrame, int lastFrame)
{
	int shrink = opts.pPyramid.pChns.shrink;
	int modelHt = opts.modelDsPad[0];
	int modelWd = opts.modelDsPad[1];
	int stride = opts.stride;
	float cascThr = opts.cascadeThreshold;

	cv::Mat tempThrs;
	cv::transpose(this->thrs, tempThrs);
	float *thrs = (float*)tempThrs.data;

	cv::Mat tempHs;
	cv::transpose(this->hs, tempHs);
	float *hs = (float*)tempHs.data;
	
	cv::Mat tempFids;
	cv::transpose(this->fids, tempFids);
	uint32 *fids = (uint32*) tempFids.data;
	
	cv::Mat tempChild;
	cv::transpose(this->child, tempChild);
	uint32 *child = (uint32*) tempChild.data;

	// const mwSize *fidsSize = mxGetDimensions(mxGetField(trees,0,"fids"));
	// const int nTreeNodes = (int) fidsSize[0];
 	// const int nTrees = (int) fidsSize[1];
	int nTreeNodes = this->fids.rows;
	int nTrees = this->fids.cols;
	
	int treeDepth = this->treeDepth;
	int nChns = opts.pPyramid.pChns.pColor.nChannels + opts.pPyramid.pChns.pGradMag.nChannels + opts.pPyramid.pChns.pGradHist.nChannels; 

	for (int i = firstFrame; i < imageNames.size() && i < lastFrame; i++)
	{
		clock_t frameStart = clock();

		// this conversion is necessary, so we don't apply this transformation multiple times, which would break the image inside chnsPyramid
		cv::Mat image = cv::imread(dataSetDirectoryName + '/' + imageNames[i]);
		cv::Mat I;
		// which one of these conversions is best?
		//image.convertTo(I, CV_32FC3, 1.0/255.0);
		cv::normalize(image, I, 0.0, 1.0, cv::NORM_MINMAX, CV_32FC3);

		// compute feature pyramid
		std::vector<Info> framePyramid;
		framePyramid = opts.pPyramid.computeMultiScaleChannelFeaturePyramid(I);
	
		clock_t detectionStart = clock();

		BB_Array frameDetections = applyDetectorToFrame(framePyramid, shrink, modelHt, modelWd, stride, cascThr, thrs, hs, fids, child, nTreeNodes, nTrees, treeDepth, nChns);
		detections.push_back(frameDetections);
		frameDetections.clear(); //doesn't seem to make a difference

		clock_t detectionEnd = clock();
		timeSpentInDetection = timeSpentInDetection + (double(detectionEnd - detectionStart) / CLOCKS_PER_SEC);

		/*	
		// debug: shows detections before suppression
		cv::imshow("source image", I);
		showDetections(I, detections[i], "detections before suppression");
		// debug */

		detections[i] = nonMaximalSuppression(detections[i]);

		
		// debug: shows detections after suppression
		showDetections(I, detections[i], "detections after suppression");
		//printDetections(detections[i], i);
		cv::waitKey();
		// debug */

		
		// debug: saves image with embedded detections
		for (int j = 0; j<detections[i].size(); j++) 
			detections[i][j].plot(image, cv::Scalar(0,255,0));
		cv::imwrite("../datasets/results/"+imageNames[i], image);
		// debug */

		
		// experimental: do i need to clear these?
		for (int j=0; j < opts.pPyramid.computedScales; j++)
		{
			framePyramid[j].image.release();
			framePyramid[j].gradientMagnitude.release();
			framePyramid[j].gradientHistogram.clear();
		}
		image.release();
		I.release();
		// experimental */

		clock_t frameEnd = clock();
		double elapsed_secs = double(frameEnd - frameStart) / CLOCKS_PER_SEC;

		std::cout << "Frame " << i+1 << " of " << imageNames.size() << " was processed in " << elapsed_secs << " seconds.\n"; 
	}
}

// for each i suppress all j st j>i and area-overlap>overlap
BB_Array nmsMax(BB_Array source, bool greedy, double overlapArea, cv::String overlapDenominator)
{
	BB_Array result;
	BB_Array sortedArray;
	// bool discarded[source.size()];
	bool *discarded = (bool*)malloc(source.size()*sizeof(bool));

	for (int i=0; i < source.size(); i++)
	{
		sortedArray.push_back(source[i]);
		discarded[i] = false;
	}
 
	std::sort(sortedArray.begin(), sortedArray.begin()+sortedArray.size());
	
	for (int i = 0; i < sortedArray.size(); i++)
	{
		if (!greedy || !discarded[i]) // continue only if its not greedy or result[i] was not yet discarded
		{
			for (int j = i+1; j < sortedArray.size(); j++)
			{
				if (discarded[j] == false) // continue this iteration only if result[j] was not yet discarded
				{
					double xei, xej, xmin, xsMax, iw;
					double yei, yej, ymin, ysMax, ih;
					xei = sortedArray[i].topLeftPoint.x + sortedArray[i].width;
					xej = sortedArray[j].topLeftPoint.x + sortedArray[j].width;
					xmin = xej;			
					if (xei < xej)
						xmin = xei;
					xsMax = sortedArray[i].topLeftPoint.x;
					if (sortedArray[j].topLeftPoint.x > sortedArray[i].topLeftPoint.x)
						xsMax = sortedArray[j].topLeftPoint.x;
					iw = xmin - xsMax;
					yei = sortedArray[i].topLeftPoint.y + sortedArray[i].height;
					yej = sortedArray[j].topLeftPoint.y + sortedArray[j].height;
					ymin = yej;			
					if (yei < yej)
						ymin = yei;
					ysMax = sortedArray[i].topLeftPoint.y;
					if (sortedArray[j].topLeftPoint.y > sortedArray[i].topLeftPoint.y)
						ysMax = sortedArray[j].topLeftPoint.y;
					ih = ymin - ysMax;
					if (iw > 0 && ih > 0)
					{
						double o = iw * ih;
						double u;
						if (overlapDenominator == "union")
							u = sortedArray[i].height*sortedArray[i].width + sortedArray[j].height*sortedArray[j].width-o;
						else if (overlapDenominator == "min")
						{
							u = sortedArray[i].height*sortedArray[i].width;
							if (sortedArray[i].height*sortedArray[i].width > sortedArray[j].height*sortedArray[j].width)
								u = sortedArray[j].height*sortedArray[j].width;
						}
						o = o/u;
						if (o > overlapArea) // sortedArray[j] is no longer needed (is discarded)
							discarded[j] = true;
					}
				}
			}	
		}
	}
	
	// result keeps only the bounding boxes that were not discarded
	for (int i=0; i < sortedArray.size(); i++)
		if (!discarded[i])
			result.push_back(sortedArray[i]);

	free(discarded);

	return result;
}

BB_Array Detector::nonMaximalSuppression(BB_Array bbs)
{
	BB_Array result;

	//keep just the bounding boxes with scores higher than the threshold
	for (int i=0; i < bbs.size(); i++)
		if (bbs[i].score > opts.suppressionThreshold)
			result.push_back(bbs[i]);

	// bbNms would apply resize to the bounding boxes now
	// but our models dont use that, so it will be suppressed
		
	// since we just apply one detector model at a time,
	// our separate attribute would always be false
	// so the next part is simpler, nms1 follows
	
	// if there are too many bounding boxes,
	// he splits into two arrays and recurses, merging afterwards
	// this will be done if necessary
	
	// run actual nms on given bbs
	// other types might be added later
	switch (opts.suppressionType)
	{
		case MAX:
			result = nmsMax(result, false, opts.overlapArea, opts.overlapDenominator);
		break;
		case MAXG:
			result = nmsMax(result, true, opts.overlapArea, opts.overlapDenominator);
		break;
		case MS:
			// not yet implemented
		break;
		case COVER:
			// not yet implemented
		break;	
	}

	return result;
}