#include "Pyramid.h"
#include <highgui.h>

void Pyramid::readPyramid(cv::FileNode pyramidNode)
{
	pChns.readChannelFeatures(pyramidNode["pChns"]);

	scalesPerOctave = pyramidNode["nPerOct"];
	upsampledOctaves = pyramidNode["nOctUp"];
	approximatedScales = pyramidNode["nApprox"];

	providedLambdas = pyramidNode["providedLambdas"];
	lambdas[0] = pyramidNode["lambdas"][0];
	lambdas[1] = pyramidNode["lambdas"][1];
	lambdas[2] = pyramidNode["lambdas"][2];
	
	padSize = pyramidNode["padSize"];
	for (int i=0; i < padSize; i++)
		pad.push_back(pyramidNode["pad"][i]);

	minImgSize[0] = pyramidNode["minDs"][0];
	minImgSize[1] = pyramidNode["minDs"][1];

	smoothRadius = pyramidNode["smooth"];
	concatenateChannels = pyramidNode["concat"];
	completeInput = pyramidNode["complete"];

	totalTimeForRealScales = 0;
	totalTimeForApproxScales = 0;
}

// translation of the chnsPyramid.m file
std::vector<Info> Pyramid::computeFeaturePyramid(cv::Mat I, bool useCalibration)
{
	int colorChannels = pChns.pColor.nChannels;
	int histogramChannels = pChns.pGradHist.nChannels;
	clock_t start, end;

	// p.pad=round(p.pad/shrink)*shrink;
	for (int i=0; i < padSize; i++)
		pad[i] = round(pad[i]/pChns.shrink)*pChns.shrink;

	// p.minDs=max(p.minDs,shrink*4);
	minImgSize[0] = std::max(minImgSize[0], pChns.shrink);
	minImgSize[1] = std::max(minImgSize[1], pChns.shrink);

	// if(p.nApprox<0), p.nApprox=p.nPerOct-1; end
	if (approximatedScales < 0)
		approximatedScales = scalesPerOctave-1;

	// converting image from cv::Mat to float*
	float* floatImg = (float*)malloc(I.rows*I.cols*colorChannels*sizeof(float));
	cvMat2floatArray(I, floatImg, colorChannels);
	
	// convert I to appropriate color space (or simply normalize)
	// I=rgbConvert(I,cs); pChns.pColor.colorSpace='orig';
	float* convertedImage = (float*)malloc(I.rows*I.cols*colorChannels*sizeof(float));
	int previousColorSpaceType = pChns.pColor.colorSpaceType; // saves the value to be reloaded afterwards
	rgbConvert(floatImg, convertedImage, I.rows*I.cols, colorChannels, pChns.pColor.colorSpaceType, 1.0f);
	pChns.pColor.colorSpaceType = RGB;
	free(floatImg);

	// get scales at which to compute features and list of real/approx scales
	// [scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	if (!useCalibration)
		getScales(I.rows, I.cols, pChns.shrink);

	std::vector<Info> computedChannels(computedScales); 

	int new_h, new_w;
	float* I1;
	int numberOfRealScales=0;
	int i;

	// compute image pyramid [real scales]
	start = clock();
	for (i=0; i < computedScales; i = i+approximatedScales+1)
	{
		// sz=[size(I,1) size(I,2)];
		// sz1=round(sz*s/shrink)*shrink;
		new_h = round(I.rows*scales[i]/pChns.shrink)*pChns.shrink;
		new_w = round(I.cols*scales[i]/pChns.shrink)*pChns.shrink;

		//printf("Shrink: %d\n", pChns.shrink);
		//printf("Compute real scale: %f\n", scales[i]);
		//printf("Size of the image: %d x %d\n", new_w, new_h);

		if (new_h == I.rows && new_w == I.cols)
			I1 = convertedImage; // does this work?
		else // I1=imResampleMex(I,sz1(1),sz1(2),1);
		{
			I1 = (float*)malloc(new_h*new_w*colorChannels*sizeof(float));
			resample(convertedImage, I1, I.rows, new_h, I.cols, new_w, colorChannels, 1.0);
		}

		// if(s==.5 && (nApprox>0 || nPerOct==1)), I=I1;
		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
		{
			free(convertedImage);
			convertedImage = I1; // is this causing the memory leaks?
		}

		//computedChannels.insert(computedChannels.begin()+i, computeSingleScaleChannelFeatures(I1, new_h, new_w));
		computedChannels[i] = computeSingleScaleChannelFeatures(I1, new_h, new_w);

		if (I1 != convertedImage)
			free(I1);
		numberOfRealScales++;
	}
	free(convertedImage);

	end = clock();
	totalTimeForRealScales = totalTimeForRealScales + (double(end - start) / CLOCKS_PER_SEC);

	/* compute the approximated scales */
	start = clock();
	for (int i=0; i<computedScales; i++)
	{
		if (i % (approximatedScales+1) != 0)
		{
			// sz1=round(sz*scales(i)/shrink);
			int new_h = round(I.rows*scales[i]/pChns.shrink);
			int new_w = round(I.cols*scales[i]/pChns.shrink);

			//printf("Compute approx scale: %f\n", scales[i]);
			//printf("Size of the image: %d x %d\n", new_w, new_h);
			
			double ratio[3];
			int iR = 0;

			// real scale changes in i=5 and i=13 to iR=8 and iR=16  
			for (int j=0; j < numberOfRealScales-1; j++)
			{
				if (i > floor((j*(approximatedScales+1)+(j+1)*(approximatedScales+1))/2))
					iR = (j+1)*(approximatedScales+1);
			}

			int realScaleRows = computedChannels[iR].image.rows;
			int realScaleCols = computedChannels[iR].image.cols;

			// resample color channels
			ratio[0] = pow(scales[i]/scales[iR],-lambdas[0]);

			/*
			// old conversion
			int imgIndex2=0;
			float* floatImg2 = (float*)malloc(realScaleRows*realScaleCols*colorChannels*sizeof(float));
			std::vector<cv::Mat> rgb2;
		  	cv::split(computedChannels[iR].image, rgb2);
		  	for (int j=0; j < rgb2.size(); j++)
		  	{
		    	cv::Mat tempMat2;
		    	cv::transpose(rgb2[j], tempMat2);
		    	float* tempFloat2 = (float*)tempMat2.data;
		    	for (int k=0; k < realScaleRows*realScaleCols; k++)
		      		floatImg2[imgIndex2++] = tempFloat2[k];
		  	}
		  	float* tempOutput = (float*)malloc(new_h*new_w*colorChannels*sizeof(float));
		  	// old conversion */

		  	// improved conversion
		  	float* floatImg2 = (float*)malloc(realScaleRows*realScaleCols*colorChannels*sizeof(float));
		  	cvMat2floatArray(computedChannels[iR].image, floatImg2, colorChannels);
		  	float* tempOutput = (float*)malloc(new_h*new_w*colorChannels*sizeof(float));
		  	// improved conversion */
		  	
			resample(floatImg2, tempOutput, realScaleRows, new_h, realScaleCols, new_w, colorChannels, ratio[0]);

			/*
			// old conversion
			float* tempIdata = (float*)malloc(new_h*new_w*colorChannels*sizeof(float));
			floatArray2cvData(tempOutput, tempIdata, new_h, new_w, colorChannels);
			cv::Mat tempImg2(new_h, new_w, CV_32FC3);
			tempImg2.data = (uchar*)tempIdata;
			tempImg2.copyTo(computedChannels[i].image);
			//computedChannels[i].image = tempImg2.clone();
			tempImg2.release();
			free(tempIdata);
			// old conversion */

			// improved conversion
			computedChannels[i].image = floatArray2cvMat(tempOutput, new_h, new_w, colorChannels);
			// improved conversion */

			free(floatImg2); 
			free(tempOutput); 
		
			// resample gradMag channel
			ratio[1] = pow(scales[i]/scales[iR],-lambdas[1]);
			cv::Mat tempMag;
			cv::transpose(computedChannels[iR].gradientMagnitude, tempMag);
			float *floatMag = (float*)tempMag.data;
			float* tempOutput1 = (float*)malloc(new_h*new_w*1*sizeof(float));
			resample(floatMag, tempOutput1, realScaleRows, new_h, realScaleCols, new_w, 1, ratio[1]);

			/*
			// old conversion
			float* tempMdata = (float*)malloc(new_h*new_w*sizeof(float));
			floatArray2cvData(tempOutput1, tempMdata, new_h, new_w, 1);
			cv::Mat tempMag2(new_h, new_w, CV_32FC1);
			tempMag2.data = (uchar*)tempMdata;
			tempMag2.copyTo(computedChannels[i].gradientMagnitude);
			tempMag2.release();
			free(tempMdata);
			// old conversion */

			// improved conversion
			computedChannels[i].gradientMagnitude = floatArray2cvMat(tempOutput1, new_h, new_w, 1);

			free(tempOutput1);
			
			// resample histogram channels
			ratio[2] = pow(scales[i]/scales[iR],-lambdas[2]);
			computedChannels[i].gradientHistogram.reserve(histogramChannels);
			for (int k=0; k < histogramChannels; k++)
			{
				cv::Mat tempHist;
				cv::transpose(computedChannels[iR].gradientHistogram[k], tempHist);
				float *floatHist = (float*)tempHist.data;
				float* tempOutput2 = (float*)malloc(new_h*new_w*1*sizeof(float));
				resample(floatHist, tempOutput2, realScaleRows, new_h, realScaleCols, new_w, 1, ratio[2]);

				/*				
				// old conversion
				float* tempHdata = (float*)malloc(new_h*new_w*sizeof(float));
				floatArray2cvData(tempOutput2, tempHdata, new_h, new_w, 1);
				cv::Mat tempHist2(new_h, new_w, CV_32FC1);
				tempHist2.data = (uchar*)tempHdata;
				computedChannels[i].gradientHistogram.push_back(tempHist2.clone());
				tempHist2.release();
				free(tempHdata);
				// old conversion */

				// improved conversion
				computedChannels[i].gradientHistogram.push_back(floatArray2cvMat(tempOutput2, new_h, new_w, 1));

				free(tempOutput2);
			}
		}
	}
	end = clock();
	totalTimeForApproxScales = totalTimeForApproxScales + (double(end - start) / CLOCKS_PER_SEC);

	int tempPad[padSize];
	for (int i=0; i < padSize; i++)
		tempPad[i] = pad[i]/pChns.shrink;

	int smoothingRadius = pChns.pColor.smoothingRadius;

	/*
	% smooth channels, optionally pad and concatenate channels
	for i=1:nScales*nTypes, data{i}=convTri(data{i},smooth); end
	if(any(pad)), for i=1:nScales, for j=1:nTypes
	      data{i,j}=imPad(data{i,j},pad/shrink,info(j).padWith); end; end; end
	if(concat && nTypes), data0=data; data=cell(nScales,1); end
	if(concat && nTypes), for i=1:nScales, data{i}=cat(3,data0{i,:}); end; end
	*/
	//smooth channels, optionally pad and concatenate channels
	for (int i=0; i < computedScales; i++)
	{
		/*
		// debug: images before convolution
		cv::imshow("image", computedChannels[i].image);
		cv::imshow("gradMag", computedChannels[i].gradientMagnitude);
		cv::imshow("gradHist", computedChannels[i].gradientHistogram[1]);
		cv::waitKey();
		// debug */

		int h = computedChannels[i].image.rows;
		int w = computedChannels[i].image.cols;

		/*
		// old conversion 
		//this is the only substitution that changes some of the results, maybe there's an error here
		int imgIndex3=0;
		float* floatImg3 = (float*)malloc(h*w*colorChannels*sizeof(float));
		std::vector<cv::Mat> rgb3;
	  	cv::split(computedChannels[i].image, rgb3);
	  	for (int j=0; j < rgb3.size(); j++)
	  	{
	    	cv::Mat tempMat3;
	    	cv::transpose(rgb3[j], tempMat3);

	    	//floatArray2cvImage(float* source, int rows, int cols, int channels)
	    
	    	float* tempFloat3 = (float*)tempMat3.data;
	    	for (int k=0; k < h*w; k++)
	      		floatImg3[imgIndex3++] = tempFloat3[k];
	  	}
	  	// old conversion */

	  	// improved conversion
	  	float* floatImg3 = (float*)malloc(h*w*colorChannels*sizeof(float));
	  	cvMat2floatArray(computedChannels[i].image, floatImg3, colorChannels);

	  	float* tempOutput = (float*)malloc(h*w*colorChannels*sizeof(float));
		convolution(floatImg3, tempOutput, h, w, colorChannels, smoothingRadius, 1);	

		/*
		// old conversion
		float* tempIdata = (float*)malloc(h*w*colorChannels*sizeof(float));
		floatArray2cvData(tempOutput, tempIdata, h, w, 3);
		cv::Mat tempImg2(h, w, CV_32FC3);
		tempImg2.data = (uchar*)tempIdata;
		tempImg2.copyTo(computedChannels[i].image);
		tempImg2.release();
		free(tempIdata);
		// old conversion */

		// improved conversion
		computedChannels[i].image = floatArray2cvMat(tempOutput, h, w, 3);

		free(tempOutput); 
		free(floatImg3);


		cv::Mat tempMag;
		cv::transpose(computedChannels[i].gradientMagnitude, tempMag);
		float *floatMag = (float*)tempMag.data;
		float* tempOutput1 = (float*)malloc(h*w*sizeof(float));
		convolution(floatMag, tempOutput1, h, w, 1, smoothingRadius, 1);	

		/*
		// old conversion
		float* tempMdata = (float*)malloc(h*w*sizeof(float));
		floatArray2cvData(tempOutput1, tempMdata, h, w, 1);
		cv::Mat tempMag2(h, w, CV_32FC1);
		tempMag2.data = (uchar*)tempMdata;
		tempMag2.copyTo(computedChannels[i].gradientMagnitude);
		tempMag2.release();
		free(tempMdata);
		// old conversion */

		// improved conversion
		computedChannels[i].gradientMagnitude = floatArray2cvMat(tempOutput1, h, w, 1);
		
		free(tempOutput1);

		for (int j=0; j < pChns.pGradHist.nChannels; j++)
		{
			cv::Mat tempHist;
			cv::transpose(computedChannels[i].gradientHistogram[j], tempHist);
			float *floatHist = (float*)tempHist.data;
			float* tempOutput2 = (float*)malloc(h*w*sizeof(float));
			convolution(floatHist, tempOutput2, h, w, 1, smoothingRadius, 1);

			/*
			// old conversion
			float* tempHdata = (float*)malloc(h*w*sizeof(float));
			floatArray2cvData(tempOutput2, tempHdata, h, w, 1);
			cv::Mat tempHist2(h, w, CV_32FC1);
			tempHist2.data = (uchar*)tempHdata;
			tempHist2.copyTo(computedChannels[i].gradientHistogram[j]);
			tempHist2.release();
			free(tempHdata);
			// old conversion */

			// improved conversion
			computedChannels[i].gradientHistogram[j] = floatArray2cvMat(tempOutput2, h, w, 1);

			free(tempOutput2);
		}

		if (pad[0]!=0 || pad[1]!=0)
		{
			// data{i,j}=imPad(data{i,j},pad/shrink,info(j).padWith);
			// J = imPadMex( I, pad, type );
			// cv::Mat padImage(cv::Mat source, int channels, int *pad, int padSize, int type)
			computedChannels[i].image = padImage(computedChannels[i].image, 3, tempPad, padSize, REPLICATE);
			computedChannels[i].gradientMagnitude = padImage(computedChannels[i].gradientMagnitude, 1, tempPad, padSize, 0);

			for (int j=0; j < pChns.pGradHist.nChannels; j++)
				computedChannels[i].gradientHistogram[j] = padImage(computedChannels[i].gradientHistogram[j], 1, tempPad, padSize, 0);
		}
	}

	pChns.pColor.colorSpaceType = previousColorSpaceType;
	return computedChannels;
}

// translation of the chnsCompute.m file
// this procedure is currently between five to eight times slower than the original
Info Pyramid::computeSingleScaleChannelFeatures(float* source, int rows, int cols)
{
	cv::Mat gradOrientation;
	Info result;
	float *M;
	float *O;
	float *I;
	std::vector<float*> H; 

	int colorChannels = pChns.pColor.nChannels;

	/*
	[h,w,~]=size(I); cr=mod([h w],shrink);
	if(any(cr)), h=h-cr(1); w=w-cr(2); I=I(1:h,1:w,:); end
	*/
	//crop I so it becomes divisible by shrink
	int height = rows - (rows % pChns.shrink);
	int width =  cols - (cols % pChns.shrink);

	/*
	// compute color channels, old version
	float* tempI = (float*)malloc(height*width*colorChannels*sizeof(float));
	tempI = rgbConvert(source, height*width, colorChannels, pChns.pColor.colorSpaceType, 1.0f);
	I = (float*)malloc(height/pChns.pColor.smoothingRadius*width/pChns.pColor.smoothingRadius*3*sizeof(float));
	convolution(tempI, I, height, width, colorChannels, pChns.pColor.smoothingRadius, 1, CONV_TRI);
	free(tempI);
	// old version */

	// compute color channels, new version
	I = (float*)malloc(height*width*colorChannels*sizeof(float));
	convolution(source, I, height, width, colorChannels, pChns.pColor.smoothingRadius, 1);

	if (pChns.pGradHist.enabled)
	{
		std::vector<float*> tempResult = pChns.pGradMag.mGradMag(I, height, width, 0);

		if (tempResult.size() > 0)
			M = tempResult[0];
		if (tempResult.size() > 1)
			O = tempResult[1];

		if (pChns.pGradMag.normalizationRadius != 0)
		{
			float* S = (float*)malloc(height*width*sizeof(float));
			convolution(M, S, height, width, 1, pChns.pGradMag.normalizationRadius, 1);

			// normalization constant is read inside the procedure
			pChns.pGradMag.gradMagNorm(M, S, height, width);

			free(S);
		}	
	}		
	else
	{
		if (pChns.pGradMag.enabled)
		{
			// since there is no histogram channel, we just need the first return of mGradMag
			M = pChns.pGradMag.mGradMag(I, height, width, 0)[0];		

			if (pChns.pGradMag.normalizationRadius != 0)
			{
				float* S = (float*)malloc(height*width*sizeof(float));
				convolution(M, S, height, width, 1, pChns.pGradMag.normalizationRadius, 1);			
				pChns.pGradMag.gradMagNorm(M, S, height, width);
				free(S);
			}
		}	
	}

	// h=h/shrink; w=w/shrink;
	int shrinkedHeight = height/pChns.shrink;
	int shrinkedWidth = width/pChns.shrink;

	//compute gradient histogram channels
	if (pChns.pGradHist.enabled)
	{
		H = pChns.pGradHist.mGradHist(M, O, height, width, pChns.pGradMag.full);

		// this is needed because images returned from mGradHist have dimensions height/binSize and width/binSize
		int binSize = pChns.pGradHist.binSize;

		result.gradientHistogram.reserve(pChns.pGradHist.nChannels);

		for (int i=0; i < pChns.pGradHist.nChannels; i++)
		{
			float* tempH = (float*)malloc(shrinkedHeight*shrinkedWidth*1*sizeof(float));
			resample(H[i], tempH, height/binSize, shrinkedHeight, width/binSize, shrinkedWidth, 1, 1.0);

			/*
			// old conversion
			float* tempHdata = (float*)malloc(shrinkedHeight*shrinkedWidth*1*sizeof(float));
			floatArray2cvData(tempH, tempHdata, shrinkedHeight, shrinkedWidth, 1);
			cv::Mat tempHist(shrinkedHeight, shrinkedWidth, CV_32FC1);
			tempHist.data = (uchar*)tempHdata;
			result.gradientHistogram.push_back(tempHist.clone());
			tempHist.release();
			free(tempHdata); 
			// old conversion */

			// improved conversion
			result.gradientHistogram.push_back(floatArray2cvMat(tempH, shrinkedHeight, shrinkedWidth, 1));

			free(tempH);
			free(H[i]);
		}

		free(O); // alters results
	}	

	// if(p.enabled), chns=addChn(chns,I,nm,p,'replicate',h,w); end
	if (pChns.pColor.enabled)
	{
		// data=imResampleMex(data,h,w,1);
		float* tempI2 = (float*)malloc(shrinkedHeight*shrinkedWidth*colorChannels*sizeof(float));
		resample(I, tempI2, height, shrinkedHeight, width, shrinkedWidth, colorChannels, 1.0);

		/*
		// old conversion
		float* tempI2data = (float*)malloc(shrinkedHeight*shrinkedWidth*colorChannels*sizeof(float));
		floatArray2cvData(tempI2, tempI2data, shrinkedHeight, shrinkedWidth, colorChannels);
		cv::Mat tempImg(shrinkedHeight, shrinkedWidth, CV_32FC3);
		tempImg.data = (uchar*)tempI2data;
		tempImg.copyTo(result.image);
		tempImg.release();
		free(tempI2data);
		// old conversion */

		// improved conversion
		result.image = floatArray2cvMat(tempI2, shrinkedHeight, shrinkedWidth, colorChannels);

		free(tempI2);
		free(I);
	}

	// if(p.enabled), chns=addChn(chns,M,nm,p,0,h,w); end
	if (pChns.pGradMag.enabled)
	{
		float* tempM = (float*)malloc(shrinkedHeight*shrinkedWidth*pChns.pGradMag.nChannels*sizeof(float));
		resample(M, tempM, height, shrinkedHeight, width, shrinkedWidth, pChns.pGradMag.nChannels, 1.0);

		/*
		// old conversion
		float* tempMdata = (float*)malloc(shrinkedHeight*shrinkedWidth*1*sizeof(float));
		floatArray2cvData(tempM, tempMdata, shrinkedHeight, shrinkedWidth, pChns.pGradMag.nChannels);
		cv::Mat tempMag(shrinkedHeight, shrinkedWidth, CV_32FC1);
		tempMag.data = (uchar*)tempMdata;
		tempMag.copyTo(result.gradientMagnitude);
		tempMag.release();
		free(tempMdata);
		// old conversion */

		// improved conversion
		result.gradientMagnitude = floatArray2cvMat(tempM, shrinkedHeight, shrinkedWidth, pChns.pGradMag.nChannels);

		free(tempM);
		free(M);
	}

	return result;
}


void Pyramid::calibratedGetScales(int h, int w, int shrink, int boundingBoxImageWidth, int boundingBoxImageHeight, double maxPedestrianWorldHeight, 
	cv::Mat_<float> &projection, cv::Mat_<float> &homography)
{
	double curScale, lastScale, lastOctave=1.0, step;

	if (h!=0 && w!=0)
	{
		computedScales = 0;
		lastScale = findLastNecessaryScaleInAPoint(0, h, h, boundingBoxImageHeight, maxPedestrianWorldHeight, projection, homography);
		curScale = findLastNecessaryScaleInAPoint(w-boundingBoxImageWidth, h, h, boundingBoxImageHeight, maxPedestrianWorldHeight, projection, homography);
		if (curScale < lastScale)
			lastScale = curScale;

		curScale = 1.0;
		step = 0.5/scalesPerOctave;
		while (curScale >= lastScale)
		{
			scales.push_back(curScale);

			if (curScale == lastOctave/2)
			{
				lastOctave = curScale;
				step = step/2;
			}

			curScale = curScale - step;
			computedScales++;
		}

		scales_w.reserve(computedScales);
		scales_h.reserve(computedScales);
		for (int i=0; i<computedScales; i++)
		{
			scales_w.push_back(round(w*scales[i]/shrink)*shrink/w);
			scales_h.push_back(round(h*scales[i]/shrink)*shrink/h);
		}
	}
	else
		std::cout << " # getScales error: both height and width of an image need to be greater than 0!";
}


// set each scale s such that max(abs(round(sz*s/shrink)*shrink-sz*s)) is minimized 
// without changing the smaller dim of sz (tricky algebra)
void Pyramid::getScales(int h, int w, int shrink)
{
	int minSize, bgDim, smDim;
	double minSizeRatio;
	double *tempScales;
	
	if (h!=0 && w!=0)
	{
		if (h/minImgSize[0] < w/minImgSize[1])
			minSizeRatio = float(h) / minImgSize[0];
		else
			minSizeRatio = float(w) / minImgSize[1];

		// nScales = floor(nPerOct*(nOctUp+log2(min(sz./minDs)))+1);
		computedScales = floor(scalesPerOctave*(upsampledOctaves+log2(minSizeRatio))+1);		

		double s0, s1;
		double epsilon = std::numeric_limits<double>::epsilon();	

		// scales = (double*)malloc(computedScales * sizeof(double));
		tempScales = (double*)malloc(computedScales * sizeof(double));

		// if(sz(1)<sz(2)), d0=sz(1); d1=sz(2); else d0=sz(2); d1=sz(1); end
		// d0 is the small dimension and d1 is the big dimension
		if (h < w)
		{
			bgDim = w;
			smDim = h;
		}
		else
		{
			bgDim = h;
			smDim = w;
		}

		for (int i=0; i < computedScales; i++)
		{
			// scales = 2.^(-(0:nScales-1)/nPerOct+nOctUp);
			tempScales[i]=pow(2, -(float(i)/scalesPerOctave+upsampledOctaves));

			// s0=(round(d0*s/shrink)*shrink-.25*shrink)./d0;
			// s1=(round(d0*s/shrink)*shrink+.25*shrink)./d0;
			s0=(round(smDim*tempScales[i]/shrink)*shrink-.25*shrink)/smDim;
			s1=(round(smDim*tempScales[i]/shrink)*shrink+.25*shrink)/smDim;

			//what follows will substitute ss=(0:.01:1-epsilon())*(s1-s0)+s0;
			double ss[(int)round((1-epsilon)/0.01)], es0[(int)round((1-epsilon)/0.01)], es1[(int)round((1-epsilon)/0.01)];		
			int ssIndex = 0;			
			for (double j=0; j < 1-epsilon; j = j + 0.01)
			{
				ss[ssIndex] = j*(s1-s0)+s0;
				es0[ssIndex]=smDim*ss[ssIndex]; 
				es0[ssIndex]=abs(es0[ssIndex]-round(es0[ssIndex]/shrink)*shrink);
				es1[ssIndex]=bgDim*ss[ssIndex]; 
				es1[ssIndex]=abs(es1[ssIndex]-round(es1[ssIndex]/shrink)*shrink);
				ssIndex++;
			}

			double tempMaxScale[ssIndex];				
		
			// this is the max part of [~,x]=min(max(es0,es1)); 
			for (int k=0; k < ssIndex; k++)
				if (es0[k] > es1[k])
					tempMaxScale[k] = es0[k];
				else
					tempMaxScale[k] = es1[k];

			// this is the min part of [~,x]=min(max(es0,es1));
			double minScaleValue = tempMaxScale[0];
			int minScaleIndex = 0; 
			for (int j=1; j < ssIndex; j++)
				if (tempMaxScale[j] < minScaleValue)
				{
					minScaleValue = tempMaxScale[j];
					minScaleIndex = j;
				}
			
			// scales(i)=ss(x);
			// scales[i] = ss[minScaleIndex];
			tempScales[i] = ss[minScaleIndex];	
		}
		
		// just keep the values of scales[i] which are different from their neighbours
		int scalesIndex=1;
		scales.insert(scales.begin(),tempScales[0]);
		for (int i=1; i < computedScales; i++)
			if(tempScales[i] != tempScales[i-1])
			{
				scales.insert(scales.begin()+scalesIndex,tempScales[i]);
				scalesIndex++;
			}

		// this updates the value of computedScales, since some of them have been suppressed
		computedScales = scalesIndex;
		
		// scaleshw = 	[round(sz(1)*scales/shrink)*shrink/sz(1);
  		//				round(sz(2)*scales/shrink)*shrink/sz(2)]';
		for (int i=0; i<computedScales; i++)
		{
			scales_w.insert(scales_w.begin()+i, round(w*scales[i]/shrink)*shrink/w);
			scales_h.insert(scales_h.begin()+i, round(h*scales[i]/shrink)*shrink/h);
		}

		free(tempScales);
	}
	else //error, height or width of the image are wrong
		std::cout << " # getScales error: both height and width of an image need to be greater than 0!";
}


