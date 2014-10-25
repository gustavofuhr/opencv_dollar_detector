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
std::vector<Info> Pyramid::computeMultiScaleChannelFeaturePyramid(cv::Mat I)
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

	// convert I to appropriate color space (or simply normalize)
	// I=rgbConvert(I,cs); pChns.pColor.colorSpace='orig';
	float* convertedImage;

	int imgIndex=0;
	float* floatImg = (float*)malloc(I.rows*I.cols*colorChannels*sizeof(float));
	std::vector<cv::Mat> rgb;
  	cv::split(I, rgb);
  	for (int j=0; j < rgb.size(); j++)
  	{
    	cv::Mat tempMat;
    	cv::transpose(rgb[j], tempMat);
    	float* tempFloat = (float*)tempMat.data;
    	for (int i=0; i < I.rows*I.cols; i++)
      		floatImg[imgIndex++] = tempFloat[i];
  	}

	//float* floatImg = cvImage2floatArray(I, colorChannels);
	int previousColorSpaceType = pChns.pColor.colorSpaceType; // saves the value to be reloaded afterwards
	convertedImage = rgbConvert(floatImg, I.rows*I.cols, colorChannels, pChns.pColor.colorSpaceType, 1.0f);
	// if(flag==4), flag=1; end;
	pChns.pColor.colorSpaceType = RGB;

	free(floatImg);

	/*
	// debug
	std::cout << "nPerOct = " << scalesPerOctave << ", nOctUp = " << upsampledOctaves << ", minDs = (" << minImgSize[0] << "," << minImgSize[1] << 
	"), shrink = " << pChns.shrink << ", sz = (" << convertedImage.rows << "," << convertedImage.cols << ")" <<std::endl;
	// debug */
	
	/*
	nPerOct = 8, nOctUp = 0, minDs = (100,41), shrink = 4, sz = (576,720)
	*/

	// get scales at which to compute features and list of real/approx scales
	// [scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
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
			convertedImage = I1; 
		}

		computedChannels.insert(computedChannels.begin()+i, computeSingleScaleChannelFeatures(I1, new_h, new_w));
		if (I1 != convertedImage)
			free(I1);
		numberOfRealScales++;
	}
	end = clock();
	totalTimeForRealScales = totalTimeForRealScales + (double(end - start) / CLOCKS_PER_SEC);

	/*
	// calculates lambdas, maybe this will not be implemented
	int is[computedScales/(approximatedScales+1)]; //needs a better name
	int isIndex = 0;
	bool isError = false;
	// if lambdas not specified compute image specific lambdas
	if (!providedLambdas && computedScales>0 && approximatedScales>0)
	{
		for (int j=1+upsampledOctaves*scalesPerOctave; j < computedScales; j = j+approximatedScales+1)
		{
			//since j is being used as the index in this for statement,
			//the next line is equivalent to is=1+nOctUp*nPerOct:nApprox+1:nScales;
			is[isIndex] = j;	

			//the next if statement substitutes assert(length(is)>=2)
			if (is[isIndex] < 2)
				isError = true;

			isIndex++;	
		}

		//next line substitutes if(length(is)>2), is=is(2:3); end	
		//because isIndex is going to carry the number of iteractions
		//made in the previous for statement
		if (isIndex > 1);
			//is=is(2:3)

		//channelTypes is the substitute for the nTypes value
		double *f0 = (double*)memset(f0,0,channelTypes*sizeof(double));		
		double *f1 = (double*)memset(f1, 0, channelTypes*sizeof(double));

		// only two of the channel types seem to be considered here
		// since chnsCompute computes the histogram channel after
		// the other two, its ommited
		for (int j=0; j<channelTypes; j++)
		{
			double sum=0;
			int i;
			for (i=0; i<computedChannels[j].image.dims; i++)
				sum = sum + computedChannels[j].image.at<double>(j,i);			
			f0[j] = sum/i;

			sum=0;
			int k;
			for (k=0; k<computedChannels[j].gradientMagnitude.dims; k++)
			{
				sum = sum + computedChannels[j].gradientMagnitude.at<double>(j,k);	
			}
			f1[j] = sum/k;
		}
		
		//next, we use the results of f0 and f1 to compute lambdas
		//lambdas = - log2(f0./f1) / log2(scales(is(1))/scales(is(2)));
		//this will be revisited at a later time
	}
	*/

	/*
	% compute image pyramid [approximated scales]
	for i=isA
	  iR=isN(i); sz1=round(sz*scales(i)/shrink);
	  for j=1:nTypes, 
	  	ratio=(scales(i)/scales(iR)).^-lambdas(j);
	    data{i,j}=imResampleMex(data{iR,j},sz1(1),sz1(2),ratio); 
	   end
	end
	*/
	start = clock();
	for (int i=0; i<computedScales; i++)
	{
		if (i % (approximatedScales+1) != 0)
		{
			// sz1=round(sz*scales(i)/shrink);
			int new_h = round(I.rows*scales[i]/pChns.shrink);
			int new_w = round(I.cols*scales[i]/pChns.shrink);
			
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
			resample(floatImg2, tempOutput, realScaleRows, new_h, realScaleCols, new_w, colorChannels, ratio[0]);
			computedChannels[i].image = floatArray2cvImage(tempOutput, new_h, new_w, colorChannels);
			//free(tempOutput); // this is probably necessary, but slightly changes the results so is suppressed for now
		
			// resample gradMag channel
			ratio[1] = pow(scales[i]/scales[iR],-lambdas[1]);
			cv::Mat tempMag;
			cv::transpose(computedChannels[iR].gradientMagnitude, tempMag);
			float *floatMag = (float*)tempMag.data;
			float* tempOutput1 = (float*)malloc(new_h*new_w*1*sizeof(float));
			resample(floatMag, tempOutput1, realScaleRows, new_h, realScaleCols, new_w, 1, ratio[1]);
			computedChannels[i].gradientMagnitude = floatArray2cvImage(tempOutput1, new_h, new_w, 1);
			free(tempOutput1);
			
			// resample histogram channels
			ratio[2] = pow(scales[i]/scales[iR],-lambdas[2]);
			for (int k=0; k < histogramChannels; k++)
			{
				cv::Mat tempHist;
				cv::transpose(computedChannels[iR].gradientHistogram[k], tempHist);
				float *floatHist = (float*)tempHist.data;
				float* tempOutput2 = (float*)malloc(new_h*new_w*1*sizeof(float));
				resample(floatHist, tempOutput2, realScaleRows, new_h, realScaleCols, new_w, 1, ratio[2]);
				computedChannels[i].gradientHistogram.push_back(floatArray2cvImage(tempOutput2, new_h, new_w, 1));
				free(tempOutput2);
			}
		}
	}
	end = clock();
	totalTimeForApproxScales = totalTimeForApproxScales + (double(end - start) / CLOCKS_PER_SEC);

	int tempPad[padSize];
	for (int i=0; i < padSize; i++)
		tempPad[i] = pad[i]/pChns.shrink;

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
		// debug */

		int h = computedChannels[i].image.rows;
		int w = computedChannels[i].image.cols;

		//float* floatImg3 = cvImage2floatArray(computedChannels[i].image, 3);
		//this is the only substitution that changes some of the results, maybe there's an error here
		int imgIndex3=0;
		float* floatImg3 = (float*)malloc(h*w*colorChannels*sizeof(float));
		std::vector<cv::Mat> rgb3;
	  	cv::split(computedChannels[i].image, rgb3);
	  	for (int j=0; j < rgb3.size(); j++)
	  	{
	    	cv::Mat tempMat3;
	    	cv::transpose(rgb3[j], tempMat3);
	    	float* tempFloat3 = (float*)tempMat3.data;
	    	for (int k=0; k < h*w; k++)
	      		floatImg3[imgIndex3++] = tempFloat3[k];
	  	}
	  	// */
	  	float* tempOutput = (float*)malloc(h/pChns.pColor.smoothingRadius*w/pChns.pColor.smoothingRadius*3*sizeof(float));
		convolution(floatImg3, tempOutput, h, w, 3, pChns.pColor.smoothingRadius, 1, CONV_TRI);	
		computedChannels[i].image = floatArray2cvImage(tempOutput, h, w, 3);

		cv::Mat tempMag;
		cv::transpose(computedChannels[i].gradientMagnitude, tempMag);
		float *floatMag = (float*)tempMag.data;
		float* tempOutput1 = (float*)malloc(h/pChns.pColor.smoothingRadius*w/pChns.pColor.smoothingRadius*sizeof(float));
		convolution(floatMag, tempOutput1, h, w, 1, pChns.pColor.smoothingRadius, 1, CONV_TRI);	
		computedChannels[i].gradientMagnitude = floatArray2cvImage(tempOutput1, h, w, 1);

		for (int j=0; j < pChns.pGradHist.nChannels; j++)
		{
			cv::Mat tempHist;
			cv::transpose(computedChannels[i].gradientHistogram[j], tempHist);
			float *floatHist = (float*)tempHist.data;
			float* tempOutput2 = (float*)malloc(h/pChns.pColor.smoothingRadius*w/pChns.pColor.smoothingRadius*sizeof(float));
			convolution(floatHist, tempOutput2, h, w, 1, pChns.pColor.smoothingRadius, 1, CONV_TRI);
			computedChannels[i].gradientHistogram[j] = floatArray2cvImage(tempOutput2, h, w, 1);
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

//translation of the chnsCompute.m file
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

	// compute color channels
	float* tempI = (float*)malloc(height*width*colorChannels*sizeof(float));
	tempI = rgbConvert(source, height*width, colorChannels, pChns.pColor.colorSpaceType, 1.0f);
	I = (float*)malloc(height/pChns.pColor.smoothingRadius*width/pChns.pColor.smoothingRadius*3*sizeof(float));
	convolution(tempI, I, height, width, colorChannels, pChns.pColor.smoothingRadius, 1, CONV_TRI);
	free(tempI);

	if (pChns.pGradHist.enabled)
	{
		std::vector<float*> tempResult = pChns.pGradMag.mGradMag(I, height, width, 0);

		if (tempResult.size() > 0)
			M = tempResult[0];
		if (tempResult.size() > 1)
			O = tempResult[1];

		if (pChns.pGradMag.normalizationRadius != 0)
		{
			float* S = (float*)malloc(height/pChns.pColor.smoothingRadius*width/pChns.pColor.smoothingRadius*sizeof(float));
			convolution(M, S, height, width, 1, pChns.pGradMag.normalizationRadius, 1, CONV_TRI);

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
				float* S = (float*)malloc(height/pChns.pColor.smoothingRadius*width/pChns.pColor.smoothingRadius*sizeof(float));
				convolution(M, S, height, width, 1, pChns.pGradMag.normalizationRadius, 1, CONV_TRI);			
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

		for (int i=0; i < pChns.pGradHist.nChannels; i++)
		{
			float* tempH = (float*)malloc(shrinkedHeight*shrinkedWidth*1*sizeof(float));
			resample(H[i], tempH, height/binSize, shrinkedHeight, width/binSize, shrinkedWidth, 1, 1.0);
			result.gradientHistogram.push_back(floatArray2cvImage(tempH, shrinkedHeight, shrinkedWidth, 1));
			free(tempH);
			free(H[i]);
		}
	}	

	// if(p.enabled), chns=addChn(chns,I,nm,p,'replicate',h,w); end
	if (pChns.pColor.enabled)
	{
		// data=imResampleMex(data,h,w,1);
		float* tempI = (float*)malloc(shrinkedHeight*shrinkedWidth*pChns.pColor.nChannels*sizeof(float));
		resample(I, tempI, height, shrinkedHeight, width, shrinkedWidth, pChns.pGradMag.nChannels, 1.0);
		result.image = floatArray2cvImage(I, shrinkedHeight, shrinkedWidth, pChns.pColor.nChannels);
		free(tempI);
		free(I);
	}

	// if(p.enabled), chns=addChn(chns,M,nm,p,0,h,w); end
	if (pChns.pGradMag.enabled)
	{
		float* tempM = (float*)malloc(shrinkedHeight*shrinkedWidth*1*sizeof(float));
		resample(M, tempM, height, shrinkedHeight, width, shrinkedWidth, pChns.pGradMag.nChannels, 1.0);
		result.gradientMagnitude = floatArray2cvImage(tempM, shrinkedHeight, shrinkedWidth, pChns.pGradMag.nChannels);
		free(tempM);
		free(M);
	}

	//for now, i wont add computation of custom channels

	return result;
}

// set each scale s such that max(abs(round(sz*s/shrink)*shrink-sz*s)) is minimized 
// without changing the smaller dim of sz (tricky algebra)
// getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);
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
	}
	else //error, height or width of the image are wrong
		std::cout << " # getScales error: both height and width need to be greater than 0!";
}


