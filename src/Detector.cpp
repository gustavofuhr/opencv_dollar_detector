#include "Detector.h"

//i dont know if its gonna be needed but this is start
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

//reads the detector model from the xml model
//for now, it must be like this since the current model
//was not written by this program this will change after we are 
//set on a class structure
void Detector::importDetectorModel(cv::String fileName)
{
	cv::FileStorage xml;

	xml.open(fileName, cv::FileStorage::READ);

	if (!xml.isOpened())
	{
		std::cerr << "Failed to open " << fileName << std::endl;
	}
	else
	{
		opts.readOptions(xml["detector"]["opts"]);

		clf.readClf(xml["detector"]["clf"]);		

		xml.release();
	}
}

//std::cout << "k=" << k << ", fids[k]=" << fids[k] << ", cids[fids[k]]=" << cids[fids[k]] << ", chns1[cids[fids[k]]]=" << chns1[cids[fids[k]]] << ", thrs[k]=" << thrs[k] << std::endl;

//this procedure was just copied verbatim
inline void getChild(float *chns1, uint32 *cids, uint32 *fids, float *thrs, uint32 offset, uint32 &k0, uint32 &k)
{
  float ftr = chns1[cids[fids[k]]];
  k = (ftr<thrs[k]) ? 1 : 2;
  k0=k+=k0*2; k+=offset;
}

//bb = acfDetect1(P.data{i},Ds{j}.clf,shrink,modelDsPad(1),modelDsPad(2),opts.stride,opts.cascThr);
void Detector::acfDetect(cv::Mat image)
{
	// this is necessary, so we don't apply this transformation multiple times which would break the image inside chnsPyramid
	cv::Mat I;
	image.convertTo(I, CV_32FC3, 1.0/255.0);

	totalDetections = 0;

	// compute feature pyramid
	opts.pPyramid.computeMultiScaleChannelFeaturePyramid(I);


	// debug: test values of scales and scales hw
	std::cout << "printing scales:" << std::endl;
	for (int i=0; i < opts.pPyramid.computedScales; i++)
		std::cout << i << ": " << opts.pPyramid.scales[i] << std::endl;

	std::cout << std::endl << "printing scaleshw:" << std::endl;
	for (int i=0; i < opts.pPyramid.computedScales; i++)
		std::cout << i << ": (" << opts.pPyramid.scales_h[i] << "," << opts.pPyramid.scales_w[i] << ")" << std::endl;
	// debug */


	// debug: test size of uint32
	std::cout << std::endl << "sizeof(uint32)=" << sizeof(uint32) << ", sizeof(uint32_t)=" << sizeof(uint32_t) << std::endl << std::endl; 
	// debug */

	/*
	// debug: read pyramid from file
	cv::FileStorage xml;
	xml.open("../opencv_dollar_detector/pyramid.xml", cv::FileStorage::READ);
	if (!xml.isOpened())
	{
		std::cerr << "Failed to open pyramid.xml" << std::endl;
		std::cin.get();
	}
	// debug */

	const int shrink = opts.pPyramid.pChns.shrink;
	const int modelHt = opts.modelDsPad[0];
	const int modelWd = opts.modelDsPad[1];
	const int stride = opts.stride;
	const float cascThr = opts.cascadeThreshold;

	float *thrs = cvImage2floatArray(clf.thrs, 1);
	float *hs = cvImage2floatArray(clf.hs, 1);
	
	cv::Mat tempFids;
	cv::transpose(clf.fids, tempFids);
	uint32 *fids = (uint32*) tempFids.data;
	
	cv::Mat tempChild;
	cv::transpose(clf.child, tempChild);
	uint32 *child = (uint32*) tempChild.data;

	// const mwSize *fidsSize = mxGetDimensions(mxGetField(trees,0,"fids"));
	// const int nTreeNodes = (int) fidsSize[0];
 	// const int nTrees = (int) fidsSize[1];
	const int nTreeNodes = clf.fids.rows;
	const int nTrees = clf.fids.cols;
	
	const int treeDepth = clf.treeDepth;
	const int nChns = opts.pPyramid.pChns.pColor.nChannels + opts.pPyramid.pChns.pGradMag.nChannels + opts.pPyramid.pChns.pGradHist.nChannels; 


	// this became a simple loop because we will apply just one detector here, 
	// to apply multiple detector models you need to create multiple Detector objects. 
	for (int i = 0; i < opts.pPyramid.computedScales; i++)
	{
		// in the original file: *chnsSize = mxGetDimensions(P.data{i});
		// const int height = (int) chnsSize[0];
  		// const int width = (int) chnsSize[1];
  		// const int nChns = mxGetNumberOfDimensions(prhs[0])<=2 ? 1 : (int) chnsSize[2];
		const int height = opts.pPyramid.computedChannels[i].image.rows;
		const int width = opts.pPyramid.computedChannels[i].image.cols;

		const int height1 = (int)ceil(float(height*shrink-modelHt+1)/stride);
		const int width1 = (int)ceil(float(width*shrink-modelWd+1)/stride);

		float* chns;
		chns = features2floatArray(opts.pPyramid.computedChannels[i], height, width, 3, 1, 6);
		
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
		
		/*
		// debug: print input matrices
		int rows = opts.pPyramid.computedChannels[i].image.rows;

		print_20_elements(chns, "channel 0");
		print_100_elements(chns, rows, "channel 0");
		print_100_elements(&chns[height*width], rows, "channel 1");
		print_100_elements(&chns[height*width*2], rows, "channel 2");
		print_100_elements(&chns[height*width*3], rows, "channel 3");
		print_100_elements(&chns[height*width*4], rows, "channel 4");
		print_100_elements(&chns[height*width*5], rows, "channel 5");
		print_100_elements(&chns[height*width*6], rows, "channel 6");
		print_100_elements(&chns[height*width*7], rows, "channel 7");
		print_100_elements(&chns[height*width*8], rows, "channel 8");
		print_100_elements(&chns[height*width*9], rows, "channel 9");

		print_20_elements(thrs, "thrs");
		print_20_elements(hs, "hs");

		print_20i_elements(fids, "fids");
		print_20i_elements(child, "child");
		print_20i_elements(cids, "cids");

		std::cin.get();
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
		    if(h>cascThr) { 
		    	// debug
		    	std::cout << "detection! scale=" << i << ", c=" << c << ", r=" << r << ", h=" << h << std::endl;

		    	cs.push_back(c); rs.push_back(r); hs1.push_back(h); 
		    }
		  }
		}
		delete [] cids;
		free(chns);
		m=cs.size();

		// std::cin.get();

		// debug
		// std::cout << "acfDetect, after loop" << std::endl;

		// shift=(modelDsPad-modelDs)/2-pad;
		double shift[2];
		shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
		shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

		// debug
		std::cout << "stride=" << stride << ", shift=(" << shift[0] << "," << shift[1] << "), scaleshw[i]=(" << opts.pPyramid.scales_h[i] << "," << opts.pPyramid.scales_w[i] << ")" <<
		    ", scales[i]=" << opts.pPyramid.scales[i] << std::endl;

		for(int j=0; j<m; j++ )
		{
			BoundingBox bb;
			bb.firstPoint.x = cs[j]*stride;
			bb.firstPoint.x = (bb.firstPoint.x+shift[1])/opts.pPyramid.scales_w[i];
			bb.firstPoint.y = rs[j]*stride;
			bb.firstPoint.y = (bb.firstPoint.y+shift[0])/opts.pPyramid.scales_h[i];
			bb.height = opts.modelDs[0]/opts.pPyramid.scales[i];
			bb.width = opts.modelDs[1]/opts.pPyramid.scales[i];
			bb.score = hs1[j];
			bb.scale = i;
			detections.push_back(bb);

			// debug
			totalDetections++;
			std::cout << "totalDetections: "<<totalDetections << ", fp=("<<bb.firstPoint.x <<","<< bb.firstPoint.y<< "), height="<<bb.height << ", width="<<bb.width << ", score="<<bb.score << ", scale="<<bb.scale << std::endl;
		}

		// debug
		//std::cin.get();
		std::cout << std::endl;
		// debug */

		cs.clear();
		rs.clear();
		hs1.clear();

		// causing segmentation fault
		// bbs = bbNms(bbs, m);

		// debug
		// std::cout << "acfDetect, end of loop, i=" << i << ", computedScales=" << opts.pPyramid.computedScales << std::endl;

	}


	// debug: print all detections
	cv::Mat img = cv::imread("../opencv_dollar_detector/frame0254.png");
	for (int i = 0; i<detections.size(); ++i) 
	{
		std::cout << "Detection fp: " << detections[i].firstPoint << std::endl;
		std::cout << "Detection (w, h): " << detections[i].width << "  " << detections[i].height << std::endl;
		std::cout << "Score: "  << detections[i].score << std::endl;
		std::cout << "Scale: "  << detections[i].scale << std::endl;

		detections[i].plot(img, cv::Scalar(0,255,0));
	}
	cv::imshow("all detections", img);
	cv::waitKey();
	// debug */

	// debug
	// xml.release();
}

BB_Array Detector::bbNms(BB_Array bbs, int size)
{
	BB_Array result;
	int j;

	//keep just the bounding boxes with scores higher than the threshold
	for (int i=0; i < size; i++)
	{
		if (bbs[i].score > opts.pNms.threshold)
		{
			result.push_back(bbs[i]);
			j++;
		}
	}

	// bbNms would apply resize to the bounding boxes now
	// but our models dont use that, so it will be suppressed
		
	// since we just apply one detector model at a time,
	// our separate attribute would always be false
	// so the next part is simpler, nms1 follows
	
	// if there are too many bounding boxes,
	// he splits into two arrays and recurses, merging afterwards
	// this will be done if needed
	
	// run actual nms on given bbs
	// other types might be added later
	if (opts.pNms.type == "maxg")
		result = nmsMax(result, size, true);

	return result;
}

// for each i suppress all j st j>i and area-overlap>overlap
BB_Array Detector::nmsMax(BB_Array source, int size, bool greedy)
{
	BB_Array result;
	BB_Array temp = source;
	
	// sort the result array by score
	// std::sort(std::begin(temp), std::end(temp));
	
	for (int i = 0; i < size; i++)
	{
		// continue only if its not greedy or result[i] was not yet discarded
		for (int j = i+1; j < size; j++)
		{
			// continue this loop only if result[j] was not yet discarded
			double xei, xej, xmin, xsMax, iw;
			double yei, yej, ymin, ysMax, ih;
			xei = source[i].firstPoint.x + source[i].width;
			xej = source[j].firstPoint.x + source[j].width;
			xmin = xej;			
			if (xei < xej)
				xmin = xei;
			xsMax = source[i].firstPoint.x;
			if (source[j].firstPoint.x > source[i].firstPoint.x)
				xsMax = source[j].firstPoint.x;
			iw = xmin - xsMax;
			yei = source[i].firstPoint.y + source[i].height;
			yej = source[j].firstPoint.y + source[j].height;
			ymin = yej;			
			if (yei < yej)
				ymin = yei;
			ysMax = source[i].firstPoint.y;
			if (source[j].firstPoint.y > source[i].firstPoint.y)
				ysMax = source[j].firstPoint.y;
			ih = ymin - ysMax;
			if (iw > 0 && ih > 0)
			{
				double o = iw * ih;
				double u;
				if (opts.pNms.ovrDnm == "union")
					u = source[i].height*source[i].width + source[j].height*source[j].width-o;
				else if (opts.pNms.ovrDnm == "min")
				{
					u = source[i].height*source[i].width;
					if (source[i].height*source[i].width > source[j].height*source[j].width)
						u = source[j].height*source[j].width;
				}
				o = o/u;
				if (o > opts.pNms.overlap);
					// source[j] is no longer needed (is discarded)
			}
		}	
	}
	
	// result keeps only the bounding boxes that were not discarded

	return result;
}
