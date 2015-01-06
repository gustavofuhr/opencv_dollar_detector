#include "Detector.h"

#include <fstream>


OddConfig::OddConfig(std::string config_file) :
	resizeImage(1.0),
	firstFrame(0),
	lastFrame(99999),
	saveFrames(false),
	saveLog(false),
	useCalibration(false)
{

	std::ifstream in_file;
	in_file.open(config_file);

	std::string token;
	while (in_file >> token) {
		if (token == "resizeImage") in_file >> resizeImage;
		else if (token == "firstFrame") in_file >> firstFrame;
		else if (token == "lastFrame") in_file >> lastFrame;
		else if (token == "detectorFileName") in_file >> detectorFileName;
		else if (token == "dataSetDirectory") in_file >> dataSetDirectory;
		else if (token == "saveFrames") {
			std::string sbool;
			in_file >> sbool;
			saveFrames = (sbool == "true");
		}
		else if (token == "outputFolder") in_file >> outputFolder;
		else if (token == "saveLog") {
			std::string sbool;
			in_file >> sbool;
			saveLog = (sbool == "true");
		}
		else if (token == "logFilename") in_file >> logFilename;
		else if (token == "useCalibration")  {
			std::string sbool;
			in_file >> sbool;
			useCalibration = (sbool == "true");
		}
		else if (token == "calibrationP") {
			float *dP = new float[12];
			for (int i=0; i<12; ++i) {
				in_file >> dP[i];
				std::cout << dP[i] << std::endl;
			}

			calibrationP = new cv::Mat_<float>(3, 4, dP);
		}
		else if (token == "supressionThreshold") in_file >> supressionThreshold;
		else if (token == "classifierThreshold") in_file >> classifierThreshold;
		else {
			std::cout << "Token not recognized!" << std::endl;
		}
	}
}




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


inline int getCandidateImageHeight(int u, int v, cv::Mat_<float> &P3, cv::Mat_<float> &inv_H_3rd_line, float worldHeight) {
	double imageHeight = v - (v + P3(1)*worldHeight*(inv_H_3rd_line(2) + inv_H_3rd_line(0)*u + inv_H_3rd_line(1)*v))/(P3(2)*worldHeight*(inv_H_3rd_line(2) + inv_H_3rd_line(0)*u + inv_H_3rd_line(1)*v) + 1);

	return (int)imageHeight;
}

BB_Array* Detector::generateCandidatesFaster(int imageHeight, int imageWidth, int shrink, cv::Mat_<float> &P, double *maxHeight,
							cv::Mat &im_debug, float meanHeight/* = 1.7m*/, float stdHeight/* = 0.1m*/, float factorStdHeight/* = 2.0*/) 
{

	// there is a set of parameters here that are hard coded, but should
	// be read from a file or something...
	cv::Mat_<float> P3 = P.col(2);

	std::cout << "P3: " << P3 << std::endl;
	float aspectRatio = 0.41;
	float minImageHeight = 80;

	float stepHeight = 100;
	int totalCandidates = 0;

	BB_Array *candidates = new BB_Array();
	double max_h = 0;

	// assemble the third line of the inverse of the homography
	cv::Mat_<float> H(3, 3, 0.0);
	H(0,0) = P(0,0); H(1,0) = P(1,0); H(2,0) = P(2,0);
	H(0,1) = P(0,1); H(1,1) = P(1,1); H(2,1) = P(2,1);
	H(0,2) = P(0,3); H(1,2) = P(1,3); H(2,2) = P(2,3);
	
	std::cout << "H: " << H << std::endl;
	cv::Mat_<float> H_inv = H.inv();

	// std::cout << "H_inv " << H_inv << std::endl;
	// cv::Mat_<float> inv_H_3rd_line = H_inv.row(2);

	// std::cout << " inv_H_3rd_line  " <<  inv_H_3rd_line  << std::endl;

	// create foot points using the pixels of the image
	for (int u = 0; u < imageWidth; u+=shrink) {
		for (int v = minImageHeight; v < imageHeight; v+=shrink ) {

			float Xw = (H_inv(0,0)*u + H_inv(0,1)*v + H_inv(0,2))/(H_inv(2,0)*u + H_inv(2,1)*v + H_inv(2,2));
			float Yw = (H_inv(1,0)*u + H_inv(1,1)*v + H_inv(1,2))/(H_inv(2,0)*u + H_inv(2,1)*v + H_inv(2,2));

			// now create candidates at different heights
			for (float h = -stdHeight * factorStdHeight; h <= stdHeight * factorStdHeight; h+= stepHeight) {
				float wHeight = meanHeight + h;

				int head_v = (int)((Xw*P(1,0) + Yw*P(1,1) + wHeight*P(1,2) + P(1,3))/(Xw*P(2,0) + Yw*P(2,1) + wHeight*P(2,2) + P(2,3)));
				int i_height = v - head_v;

				if (i_height >= minImageHeight) {
					int head_u = (int)((Xw*P(0,0) + Yw*P(0,1) + wHeight*P(0,2) + P(0,3))/(Xw*P(2,0) + Yw*P(2,1) + wHeight*P(2,2) + P(2,3)));

					BoundingBox bb;

					int i_width = i_height*aspectRatio;

				    bb.topLeftPoint.x = (int)(((u + head_u)/2.0) - (i_width/2.0));
				    bb.topLeftPoint.y = head_v;
				    bb.width          = i_width;
				    bb.height         = i_height;
				    bb.world_height   = wHeight;

				    if (bb.topLeftPoint.x >= 0 && bb.topLeftPoint.x+bb.width < imageWidth && 
						    bb.topLeftPoint.y >= 0 && bb.topLeftPoint.y+bb.height < imageHeight &&
						    bb.height >= minImageHeight) {
						if (bb.height > max_h) 
							max_h = bb.height;
						candidates->push_back(bb);
					}
				}

				totalCandidates++;

			}
		}
	}
	
	
	std::cout << "Total candidates: " << totalCandidates << std::endl;
 
	*maxHeight = max_h;

	return candidates;
}





BB_Array* Detector::generateCandidates(int imageHeight, int imageWidth, cv::Mat_<float> &P, double *maxHeight,
							float meanHeight/* = 1.7m*/, float stdHeight/* = 0.1m*/, float factorStdHeight/* = 2.0*/) 
{

	// there is a set of parameters here that are hard coded, but should
	// be read from a file or something...
	cv::Mat_<float> area(2,2);
	area(0,0) = -20000;
	area(0,1) =  20000;
	area(1,0) = -20000;
	area(1,1) =  20000;

	float step = 200;
	float aspectRatio = 0.41;

	float stepHeight = 100;

	float minImageHeight = 80;

	int totalCandidates = 0;

	BB_Array *candidates = new BB_Array();
	double max_h = 0;
	for (float x = area(0,0); x < area(0,1); x += step) {
		for (float y = area(1,0); y < area(1,1); y += step) {
			// for all points in the ground plane (according to the area), I try all
			// the heights that I need to

			for (float h = -stdHeight * factorStdHeight; h <= stdHeight * factorStdHeight; h+= stepHeight) {

				cv::Mat_<float> w_feet_point(4, 1, 0.0);
			    w_feet_point(0) = x;
			    w_feet_point(1) = y;
			    w_feet_point(2) = 0.0;
			    w_feet_point(3) = 1.0;
				cv::Mat_<float> i_feet_point = world2image(w_feet_point, P);

				if (i_feet_point(0) >= 0 && i_feet_point(0) < imageWidth &&
					i_feet_point(1) >= 0 && i_feet_point(1) < imageHeight) {


				float wHeight = meanHeight + h;

				cv::Point2f wPoint(x, y);

				


				BoundingBox bb = wcoord2bbox(wPoint, P, wHeight, aspectRatio);

				// only put if it is inside the visible image
				if (bb.topLeftPoint.x >= 0 && bb.topLeftPoint.x+bb.width < imageWidth && 
					    bb.topLeftPoint.y >= 0 && bb.topLeftPoint.y+bb.height < imageHeight &&
					    bb.height >= minImageHeight) {
					if (bb.height > max_h) 
						max_h = bb.height;
					candidates->push_back(bb);
				}

				totalCandidates++;
				}
			}
			
		}
	}

	std::cout << "Total candidates: " << totalCandidates << std::endl;
 
	*maxHeight = max_h;

	return candidates;
}

/*
    This function return the scale in which the pyramid will be more fitted
*/
int Detector::findClosestScaleFromBbox(std::vector<Info> &pyramid, BoundingBox &bb,
												int modelHeight, int imageHeight)
{
	// actually here is the size of the the image that changes, the model stays the same
	// to see the best fit for the bounding box, one must find the relation between the original
	// and then find the closest to the size of the bounding box
	float min_dist = imageHeight;
	int i_min = -1;

	for (int i = 0; i < opts.pPyramid.computedScales; i++) {
		float diff = fabs(opts.modelDs[0]/opts.pPyramid.scales[i] - bb.height);

		if (diff < min_dist) {
			i_min = i;
			min_dist = diff;
		}
		else break;

	}

	return i_min;

}

void Detector::bbTopLeft2PyramidRowColumn(int *r, int *c, BoundingBox &bb, int modelHt, int modelWd, int ith_scale, int stride) {


	double s1, s2, sw, sh, tlx, tly;

	s1 = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
	s2 = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

	sw = opts.pPyramid.scales_w[ith_scale];
	sh = opts.pPyramid.scales_h[ith_scale];

	tlx = (double)bb.topLeftPoint.x;
	tly = (double)bb.topLeftPoint.y;

	double fc = (sw*tlx - s2)/(double)stride;
	double fr = (sh*tly - s1)/(double)stride;

	*r = (int)fr;
	*c = (int)fc;
}


BoundingBox Detector::pyramidRowColumn2BoundingBox(int r, int c,  int modelHt, int modelWd, int ith_scale, int stride) {

	double shift[2];
	shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
	shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

	BoundingBox bb;
	bb.topLeftPoint.x = c*stride;
	bb.topLeftPoint.x = (bb.topLeftPoint.x+shift[1])/opts.pPyramid.scales_w[ith_scale];
	bb.topLeftPoint.y = r*stride;
	bb.topLeftPoint.y = (bb.topLeftPoint.y+shift[0])/opts.pPyramid.scales_h[ith_scale];
	bb.height = opts.modelDs[0]/opts.pPyramid.scales[ith_scale];
	bb.width = opts.modelDs[1]/opts.pPyramid.scales[ith_scale];
	bb.scale = ith_scale;

	return bb;
}


BB_Array Detector::applyDetectorToFrameSmart(std::vector<Info> pyramid, BB_Array* bbox_candidates, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, float *hs, 
										uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns, int imageWidth, int imageHeight, cv::Mat_<float> &P, cv::Mat &debug_image)
{
	BB_Array result;


	// std::cout << "stride: " << stride << std::endl;
	// std::cout << "shrink: " << shrink << std::endl;


	
	std::cout << "Number of candidates: " << bbox_candidates->size() << std::endl;

	// create one candidate only for debug
	// BB_Array bbox_candidates;
	// BoundingBox bb;
	// bb.topLeftPoint.x = 486;
	// bb.topLeftPoint.y = 148;
	// bb.width = 57;
	// bb.height = 104;
	// bbox_candidates.push_back(bb);


	// plot the candidates, only for DEBUG
	// cv::waitKey(1000);
	// for (int i = 0; i < bbox_candidates.size(); ++i) {
		
	// 	std::cout << bbox_candidates[i].topLeftPoint << std::endl;
	//  	bbox_candidates[i].plot(debug_image, cv::Scalar(0, 255, 0));
	//  	cv::imshow("candidates", debug_image);
	//  	cv::waitKey(0);
	// }
	//cv::waitKey(40);

	std::cout << "opts.pPyramid.computedScales: " << opts.pPyramid.computedScales << std::endl;

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




	float max_h = -1000;

	for (int i = 0; i < bbox_candidates->size(); ++i) {

		// see which scale is best suited to the candidate
		int ith_scale = findClosestScaleFromBbox(pyramid, (*bbox_candidates)[i], modelHt, imageHeight);
		
		int height = pyramid[ith_scale].image.rows;
		int width = pyramid[ith_scale].image.cols;
		
		// std::cout << "Original bbox " << bbox_candidates[i].topLeftPoint.x << " " << bbox_candidates[i].topLeftPoint.y << " size: " 
		// 	<< bbox_candidates[i].width << " x " << bbox_candidates[i].height << std::endl;
		// int row, col;
		// bbTopLeft2PyramidRowColumn(&row, &col, bbox_candidates[i], modelHt, modelWd, ith_scale, stride);
		// std::cout << "r " << row << " c " << col << std::endl;
		// BoundingBox testbb = pyramidRowColumn2BoundingBox(row, col, modelHt, modelWd, ith_scale, stride);
		// std::cout << "New bbox " << testbb.topLeftPoint.x << " " << testbb.topLeftPoint.y << " size: " 
		// 	<< testbb.width << " x " << testbb.height << std::endl;
		// testbb.plot(debug_image, cv::Scalar(0, 0, 255));
		// bbox_candidates[i].plot(debug_image, cv::Scalar(0, 255, 0));
		// testbb.plot(debug_image, cv::Scalar(0, 0, 255));
		// cv::imshow("candidates", debug_image);
		// cv::waitKey(0);
		// exit(42);	

		// r and c are defined by the candidate itself
		int r, c;
		bbTopLeft2PyramidRowColumn(&r, &c, (*bbox_candidates)[i], modelHt, modelWd, ith_scale, stride);

		float h=0, *chns1=scales_chns[ith_scale]+(r*stride/shrink) + (c*stride/shrink)*height;
	    
	    if( treeDepth==1 ) {
	      // specialized case for treeDepth==1
	      for( int t = 0; t < nTrees; t++ ) {
	        uint32 offset=t*nTreeNodes, k=offset, k0=0;
	        getChild(chns1,scales_cids[ith_scale],fids,thrs,offset,k0,k);
	        h += hs[k]; if( h<=cascThr ) break;
	      }
	    } else if( treeDepth==2 ) {
	      // specialized case for treeDepth==2
	      for( int t = 0; t < nTrees; t++ ) {
	        uint32 offset=t*nTreeNodes, k=offset, k0=0;

	        getChild(chns1,scales_cids[ith_scale],fids,thrs,offset,k0,k);
	        getChild(chns1,scales_cids[ith_scale],fids,thrs,offset,k0,k);
	        
	        h += hs[k]; if( h<=cascThr ) break;
	      }
	    } else if( treeDepth>2) {
	      // specialized case for treeDepth>2
	      for( int t = 0; t < nTrees; t++ ) {
	        uint32 offset=t*nTreeNodes, k=offset, k0=0;
	        for( int i=0; i<treeDepth; i++ )
	          getChild(chns1,scales_cids[ith_scale],fids,thrs,offset,k0,k);
	        h += hs[k]; if( h<=cascThr ) break;
	      }
	    } else {
	      // general case (variable tree depth)
	      for( int t = 0; t < nTrees; t++ ) {
	        uint32 offset=t*nTreeNodes, k=offset, k0=k;
	        while( child[k] ) {
	          float ftr = chns1[scales_cids[ith_scale][fids[k]]];
	          k = (ftr<thrs[k]) ? 1 : 0;
	          k0 = k = child[k0]-k+offset;
	        }
	        h += hs[k]; if( h<=cascThr ) break;
	      }
	    }

	    if (h > max_h) {
	    	max_h = h;
	    }


	    if(h>config.classifierThreshold){
			// std::cout << h << std::endl;
			// std::cout << "hey" << std::endl;
			//cv::imshow("results", debug_image);
			BoundingBox detection((*bbox_candidates)[i]);
			detection.score = h;
			detection.scale = ith_scale;

	    	result.push_back(detection);
	    	// bbox_candidates[i].plot(debug_image, cv::Scalar(0, 255, 0));
	    	//cv::waitKey(100);
	    }
	
	}

	// std::cout << "max h: " << max_h << std::endl;
	// std::cout << "cascThr: " << cascThr << std::endl;
	
	//cv::waitKey(0);

	//free the memory used to pre-allocate indexes
	for (int i=0; i < opts.pPyramid.computedScales; i++) {
		free(scales_chns[i]);
		free(scales_cids[i]);
	}

	return result;
}




BB_Array Detector::applyDetectorToFrame(std::vector<Info> pyramid, int shrink, int modelHt, int modelWd, int stride, float cascThr, float *thrs, float *hs, 
										uint32 *fids, uint32 *child, int nTreeNodes, int nTrees, int treeDepth, int nChns)
{
	BB_Array result;


	printf("Number of scales: %d\n", opts.pPyramid.computedScales);
	printf("Model size: %d x %d\n", modelWd, modelHt);
	// this became a simple loop because we will apply just one detector here, 
	// to apply multiple detector models you need to create multiple Detector objects. 
	float max_h = -10000;

	int n_candidates = 0;

	for (int i = 0; i < opts.pPyramid.computedScales; i++)
	{
		// in the original file: *chnsSize = mxGetDimensions(P.data{i});
		// const int height = (int) chnsSize[0];
  		// const int width = (int) chnsSize[1];
  		// const int nChns = mxGetNumberOfDimensions(prhs[0])<=2 ? 1 : (int) chnsSize[2];
		int height = pyramid[i].image.rows;
		int width = pyramid[i].image.cols;
		printf("Size of the downsampled image: %d x %d\n", width, height);
		printf("Aspect ratio of the image: %f\n", height/(float)width);

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
	  	for( int z=0; z<nChns; z++ ) {
	    	for( int c=0; c<modelWd/shrink; c++ ) {
	      		for( int r=0; r<modelHt/shrink; r++ ) {
	        		cids[m++] = z*width*height + c*height + r;
	        		// std::cout << "cids[m] " << cids[m-1] << std::endl;
	        	}
	        }
	    }

		/*
		// debug: prints values of several variables, all of these return correct results
		// shrink=4, modelHt=128, modelWd=64, stride=4, cascThr=-1.000000, treeDepth=2
		// height=152, width=186, nChns=10, nTreeNodes=7, nTrees=2048, height1=121, width1=171, nFtrs=5120
		std::cout << "shrink=" << shrink << ", modelHt=" << modelHt << ", modelWd=" << modelWd << ", stride=" << stride << ", cascThr=" << cascThr << ", treeDepth=" << treeDepth <<  ", modelDs=(" <<
			opts.modelDs[0] << "," << opts.modelDs[1] << ")" << std::endl;
		std::cout << "height=" << height << ", width=" << width << ", nChns=" << nChns <<  ", nTreeNodes=" << nTreeNodes << ", nTrees=" << nTrees << ", height1=" << height1 << 
			", width1=" << width1 << ", nFtrs=" << nFtrs << std::endl;
		// debug */

		n_candidates += width1*height1;

		// apply classifier to each patch
  		std::vector<int> rs, cs; std::vector<float> hs1;
  		for( int c=0; c<width1; c++ ) 
  		{
  			for( int r=0; r<height1; r++ ) 
  			{
  				bool debug = false;
  				if (r==37 && c==121) {
  					//debug = true;
  				}

			    float h=0, *chns1=chns+(r*stride/shrink) + (c*stride/shrink)*height;
			    if (debug) {
			    	std::cout << "*chns1 not smart " << *chns1 << std::endl;
			    }
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
			        h += hs[k]; if( h<=config.classifierThreshold ) break;
			      }
		    }

		    if (h>max_h)
		    	max_h = h;

		    if(h>config.classifierThreshold) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }

		    if (debug) {
  					double shift[2];
					shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
					shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

					BoundingBox bb;
					bb.topLeftPoint.x = c*stride;
					bb.topLeftPoint.x = (bb.topLeftPoint.x+shift[1])/opts.pPyramid.scales_w[i];
					bb.topLeftPoint.y = r*stride;
					bb.topLeftPoint.y = (bb.topLeftPoint.y+shift[0])/opts.pPyramid.scales_h[i];
					bb.height = opts.modelDs[0]/opts.pPyramid.scales[i];
					bb.width = opts.modelDs[1]/opts.pPyramid.scales[i];
					bb.scale = i;


  				}
		  }
		}
		delete [] cids;
		free(chns);
		m=cs.size();

		// shift=(modelDsPad-modelDs)/2-pad;
		// double shift[2];
		// shift[0] = (modelHt-double(opts.modelDs[0]))/2-opts.pPyramid.pad[0];
		// shift[1] = (modelWd-double(opts.modelDs[1]))/2-opts.pPyramid.pad[1];

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
	}

	std::cout << "max_h " << max_h << std::endl;
	std::cout << "Number of candidates (not smart) : " << n_candidates << std::endl; 

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


	if (config.resizeImage != 1.0) {
		float s = config.resizeImage;
		float scale_matrix[9] = {s, 0.0, 0.0, 0.0, s, 0.0, 0.0, 0.0, 1.0};

		cv::Mat_<float> S(3, 3, scale_matrix);
		*(config.calibrationP) = S * *(config.calibrationP);
	}

	std::ofstream out_file;
	if (config.saveLog) {
		std::cout << "It's going to save a log file in " << config.logFilename << std::endl;
		out_file.open(config.logFilename);
		if (out_file.fail()) {
        	std::cerr << "open failure: " << strerror(errno) << '\n';
        }
        
	}

	for (int i = firstFrame; i < imageNames.size() && i < lastFrame; i++)
	{
		clock_t frameStart = clock();

		// this conversion is necessary, so we don't apply this transformation multiple times, which would break the image inside chnsPyramid
		cv::Mat image = cv::imread(dataSetDirectoryName + '/' + imageNames[i]);

		// if resize image is set different to 1.0, resize before computing the pyramid
		std::cout << "Resizing the image, if required..." << std::endl;
		if (config.resizeImage != 1.0) {
			cv::resize(image, image, cv::Size(), config.resizeImage, config.resizeImage);
		}

		cv::Mat I;
		// which one of these conversions is best?
		//image.convertTo(I, CV_32FC3, 1.0/255.0);
		cv::normalize(image, I, 0.0, 1.0, cv::NORM_MINMAX, CV_32FC3);

		// compute feature pyramid
		std::vector<Info> framePyramid;
		//framePyramid = opts.pPyramid.computeMultiScaleChannelFeaturePyramid(I);
	
		clock_t detectionStart = clock();

		BB_Array frameDetections;
		if (config.useCalibration) {
			double maxHeight = 0;
			
			clock_t candidatesStart = clock();
			BB_Array *bbox_candidates = generateCandidatesFaster(image.rows, image.cols, 4, *(config.calibrationP), &maxHeight, I);
			clock_t candidatesEnd = clock();

			double elapsed_secsc = double(candidatesEnd - candidatesStart) / CLOCKS_PER_SEC;
			std::cout << "Time to create candidates: " << elapsed_secsc << std::endl;
			std::cout << "Max candidate height in the image: " << maxHeight <<  std::endl;
			framePyramid = opts.pPyramid.computeMultiScaleChannelFeaturePyramidSmart(I, opts.modelDs, maxHeight);
			//exit(42);
			frameDetections = applyDetectorToFrameSmart(framePyramid, bbox_candidates, shrink, modelHt, modelWd, stride, cascThr, thrs, hs, fids, child, nTreeNodes, nTrees, treeDepth, nChns, image.cols, image.rows, *(config.calibrationP), image);

			free(bbox_candidates);
		}
		else {
			framePyramid = opts.pPyramid.computeMultiScaleChannelFeaturePyramid(I);
			frameDetections = applyDetectorToFrame(framePyramid, shrink, modelHt, modelWd, stride, cascThr, thrs, hs, fids, child, nTreeNodes, nTrees, treeDepth, nChns);		
		}
			
		
		//BB_Array frameDetections = applyDetectorToFrame(framePyramid, shrink, modelHt, modelWd, stride, cascThr, thrs, hs, fids, child, nTreeNodes, nTrees, treeDepth, nChns);
		detections.push_back(frameDetections);
		frameDetections.clear(); //doesn't seem to make a difference

		clock_t detectionEnd = clock();
		timeSpentInDetection = timeSpentInDetection + (double(detectionEnd - detectionStart) / CLOCKS_PER_SEC);

		
		// debug: shows detections before suppression
		// cv::imshow("source image", I);
		// showDetections(I, detections[i], "detections before suppression");
		// // debug 
		// cv::waitKey();


		if (config.useCalibration)
			detections[i] = nonMaximalSuppressionSmart(detections[i], 1800, 100);
		else
			detections[i] = nonMaximalSuppression(detections[i]);


		if (config.saveLog) {
			for (int x = 0; x < detections[i].size(); ++x) {
				std::cout << i+1 << " " << detections[i][x].topLeftPoint.x << " " 
						 << detections[i][x].topLeftPoint.y << " "  
						 << detections[i][x].height << " "
						 << detections[i][x].width << " "
						 << detections[i][x].score << std::endl;
				out_file << i+1 << " " << detections[i][x].topLeftPoint.x << " " 
						 << detections[i][x].topLeftPoint.y << " "  
						 << detections[i][x].height << " "
						 << detections[i][x].width << " "
						 << detections[i][x].score << std::endl;
			}
		}
		

		
		// debug: shows detections after suppression
		showDetections(I, detections[i], "detections after suppression");
		//printDetections(detections[i], i);
		cv::waitKey(30);
		// debug */

		
		// debug: saves image with embedded detections
		if (config.saveFrames) {
			for (int j = 0; j<detections[i].size(); j++) 
				detections[i][j].plot(image, cv::Scalar(0,255,0));

			std::string outfilename = config.outputFolder + imageNames[i];
			cv::imwrite(outfilename, image);
		}
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

	out_file.close();
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

	std::cout << "opts.suppressionThreshold: " << opts.suppressionThreshold << std::endl;

	//keep just the bounding boxes with scores higher than the threshold
	for (int i=0; i < bbs.size(); i++)
		if (bbs[i].score > config.supressionThreshold)
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

double gaussianFunction(double mean, double std, double x) {
	return exp(-pow(x-mean, 2)/(2*pow(std,2)));
}


BB_Array Detector::nonMaximalSuppressionSmart(BB_Array bbs, double meanHeight, double stdHeight)
{
	BB_Array result;

	//keep just the bounding boxes with scores higher than the threshold
	// for (int i=0; i < bbs.size(); i++)
	// 	if (bbs[i].score > opts.suppressionThreshold)
	// 		result.push_back(bbs[i]);

	for (int i=0; i < bbs.size(); ++i) {
		bbs[i].score = bbs[i].score*gaussianFunction(meanHeight, stdHeight, bbs[i].world_height);
		if (bbs[i].score > config.supressionThreshold)
			result.push_back(bbs[i]);
	}


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