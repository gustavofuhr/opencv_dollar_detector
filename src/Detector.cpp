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

//this procedure was just copied verbatim
void Detector::getChild( float *chns1, uint32_t *cids, uint32_t *fids, float *thrs, uint32_t offset, uint32_t &k0, uint32_t &k )
{
  	float ftr = chns1[cids[fids[k]]];
  	k = (ftr<thrs[k]) ? 1 : 2;
  	k0=k+=k0*2;
	k+=offset;
}

//bb = acfDetect1(P.data{i},Ds{j}.clf,shrink,modelDsPad(1),modelDsPad(2),opts.stride,opts.cascThr);
void Detector::acfDetect(cv::Mat image)
{
	// this is necessary, so we don't apply this transformation multiple times which would break the image inside chnsPyramid
	cv::Mat I;
	image.convertTo(I, CV_32FC3, 1.0/255.0);

	totalDetections = 0;

	//compute feature pyramid
	opts.pPyramid.computeMultiScaleChannelFeaturePyramid(I);

	//this became a simple loop because we will apply just one detector here, 
	//to apply multiple detector models you need to create multiple Detector objects. 
	for (int i = 0; i < opts.pPyramid.computedScales; i++)
	{

		// debug
		std::cout << "acfDetect, start of loop, i=" << i << ", computedScales=" << opts.pPyramid.computedScales << std::endl;

		// mxGetData(P.data{i});
		float* chns;
		float* ch1 = cvMat2floatArray(opts.pPyramid.computedChannels[i].image, 3);
		int ch1Size = opts.pPyramid.computedChannels[i].image.rows * opts.pPyramid.computedChannels[i].image.cols;
		float* ch2 = cvMat2floatArray(opts.pPyramid.computedChannels[i].gradientMagnitude, 1);
		int ch2Size = opts.pPyramid.computedChannels[i].gradientMagnitude.rows * opts.pPyramid.computedChannels[i].gradientMagnitude.cols;
		float* ch3 = cvMat2floatArray(opts.pPyramid.computedChannels[i].gradientHistogram, 1);
		int ch3Size = opts.pPyramid.computedChannels[i].gradientHistogram.rows * opts.pPyramid.computedChannels[i].gradientHistogram.cols;
		chns = (float*) malloc(ch1Size+ch2Size+ch3Size);

		// debug
		std::cout << "acfDetect, before memcpy, ch1Size=" << ch1Size << ", ch2Size=" << ch2Size << ", ch3Size=" << ch3Size << std::endl;

		memcpy(chns, ch1, ch1Size);

		// debug
		std::cout << "acfDetect, after memcpy 1" << std::endl;

		memcpy(&chns[ch1Size], ch2, ch2Size);

		// debug
		std::cout << "acfDetect, after memcpy 2" << std::endl;

		memcpy(&chns[ch1Size+ch2Size], ch3, ch3Size);

		// debug
		std::cout << "acfDetect, after memcpy" << std::endl;

		const int shrink = opts.pPyramid.pChns.shrink;
		const int modelHt = opts.modelDsPad[0];
		const int modelWd = opts.modelDsPad[1];
		const int stride = opts.stride;
		const float cascThr = opts.cascadeThreshold;

		float *thrs = cvMat2floatArray(clf.thrs, 1);
		float *hs = cvMat2floatArray(clf.hs, 1);
		uint32_t *fids = (uint32_t*) clf.fids.data;
		uint32_t *child = (uint32_t*) clf.child.data;
		const int treeDepth = clf.treeDepth;

		// in the original file: *chnsSize = mxGetDimensions(P.data{i});
		// still need to check these values properly
		const int height = opts.pPyramid.computedScales;
		// second dimension of chns
		const int width = opts.pPyramid.channelTypes;
		const int nChns = 3; 

		//nTreeNodes size of the first dimension of fids
		const int nTreeNodes = clf.fids.rows;
		//nTrees size of the second dimension of fids
		const int nTrees = clf.fids.cols;
		const int height1 = (int)ceil(float(height*shrink-modelHt+1/stride));
		const int width1 = (int)ceil(float(width*shrink-modelWd+1/stride));

		//The number of color channels
		int nChannels = opts.pPyramid.pChns.pColor.nChannels;

		//construct cids array
		int nFtrs = modelHt/shrink * modelWd/shrink * nChannels;
		uint32_t *cids = new uint32_t[nFtrs];
		int m = 0;
		for (int z = 0; z<nChannels; z++)
			for (int c = 0; c<modelWd / shrink; c++)
				for (int r = 0; r<modelHt / shrink; r++)
					cids[m++] = z*width*height + c*height + r;

		// debug
		std::cout << "acfDetect, after cids" << std::endl;

		// apply classifier to each patch
		std::vector<int> rs, cs; std::vector<float> hs1;
		int c, r;
		for (c = 0; c<width1; c++) 
			for (r = 0; r<height1; r++) 
			{
				float h = 0, *chns1 = chns + (r*stride/shrink) + (c*stride/shrink)*height;
				if (treeDepth == 1) 
				{
					// specialized case for treeDepth==1
					for (int t = 0; t < nTrees; t++) 
					{
						uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
						getChild(chns1, cids, fids, thrs, offset, k0, k);
						h += hs[k]; if (h <= cascThr) break;
					}
				}	
				else if (treeDepth == 2) 
				{
					// specialized case for treeDepth==2
					for (int t = 0; t < nTrees; t++) {
						uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
						getChild(chns1, cids, fids, thrs, offset, k0, k);
						getChild(chns1, cids, fids, thrs, offset, k0, k);
						h += hs[k]; if (h <= cascThr) break;
					}
				}
				else if (treeDepth>2) 
				{
					// specialized case for treeDepth>2
					for (int t = 0; t < nTrees; t++) 
					{
						uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
						for (int i = 0; i<treeDepth; i++)
							getChild(chns1, cids, fids, thrs, offset, k0, k);
						h += hs[k]; if (h <= cascThr) break;
					}
				}
				else 
				{
					// general case (variable tree depth)
					for (int t = 0; t < nTrees; t++) 
					{
						uint32_t offset = t*nTreeNodes, k = offset, k0 = k;
						while (child[k]) 
						{
							float ftr = chns1[cids[fids[k]]];
							k = (ftr<thrs[k]) ? 1 : 0;
							k0 = k = child[k0] - k + offset;
						}
						h += hs[k]; if (h <= cascThr) break;
					}
				}
				if (h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }
			}
			delete [] cids;
			m=cs.size();

			// debug
			std::cout << "acfDetect, after loop" << std::endl;

			double shift[2];
			shift[0] = (modelHt-opts.modelDs[0])/2-opts.pPyramid.pad[0];
			shift[1] = (modelWd-opts.modelDs[1])/2-opts.pPyramid.pad[1];

			BB_Array bbs;
			for(int j=0; j<m; j++ )
			{
				BoundingBox bb;
				bb.firstPoint.x = cs[j]*stride;
				bb.firstPoint.x = (bb.firstPoint.x+shift[1])/opts.pPyramid.scaleshw[j].x;
				bb.firstPoint.y = rs[j]*stride;
				bb.firstPoint.y = (bb.firstPoint.y+shift[0])/opts.pPyramid.scaleshw[j].y;
				bb.height = modelHt/opts.pPyramid.scales[j];
				bb.width = modelWd/opts.pPyramid.scales[j];
				bb.score = hs1[j];
				bbs.push_back(bb);
				totalDetections++;
				std::cout << "a boundingBox was added, totalDetections: " << totalDetections << std::endl;
			}

			bbs = bbNms(bbs, m);

			detections.push_back(bbs);

			// debug
			std::cout << "acfDetect, end of loop, i=" << i << ", computedScales=" << opts.pPyramid.computedScales << std::endl;

	}
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
