#include "Pyramid.h"

void Pyramid::debugWindow (cv::String name, double value)
{
	cv::Mat debugImg = cv::imread("../opencv_dollar_detector/frame0254.png");
	char str[200];
	sprintf(str,"%f",value);
	putText(debugImg, str, cv::Point2f(100,100), cv::FONT_HERSHEY_PLAIN, 2,  cv::Scalar(0,0,255,255));
	cv::imshow(name,debugImg);
	cv::waitKey();
	cv::destroyWindow(name);
}

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
	pad[0] = pyramidNode["pad"][0];
	pad[1] = pyramidNode["pad"][1];
	minImgSize[0] = pyramidNode["minDs"][0];
	minImgSize[1] = pyramidNode["minDs"][1];
	smoothRadius = pyramidNode["smooth"];
	concatenateChannels = pyramidNode["concat"];
	completeInput = pyramidNode["complete"];
}

//translation of the chnsPyramid.m file
void Pyramid::computeMultiScaleChannelFeaturePyramid(cv::Mat I)
{
	cv::Mat convertedImage;

	//if we are to allow incomplete Pyramids, we need to set some values to default.
	//for now, it wont be implemented. (lines 115-128 of chnsPyramid.m)

	//convert I to appropriate color space (or simply normalize)
	convertedImage = pChns.pColor.rgbConvert(I);
	pChns.pColor.colorSpaceType = "orig";

	// get scales at which to compute features and list of real/approx scales
	// [scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);

	computedChannels = new Info[computedScales];
	
	//next, the values os isA, isR and isN need to be set
	//isA has all the values from 1 to nScales that isR doesn't
	//j has the averages of all adjacent values of isR,
	//but starts with 0 and ends with nScales
	//isN is the full size array that has just the approximated
	//values, so each value in isR occupies nApprox cells in isN 
	
	//compute image pyramid, from chnsPyramid.m, line 144
	//there is a for statement where the index i can only assume values of the array isR
	//I guess isR might be useless, if we do it like this instead:
	int h1, w1;
	cv::Mat I1;
	int ccIndex=0;
	int isN[computedScales];

	for (int i=0; i < computedScales; i = i+approximatedScales+1)
	{
		for (int j = 0; j < approximatedScales; j++)
			isN[i+j] = i;
		h1 = round(I.rows*scales[i]/pChns.shrink)*pChns.shrink;
		w1 = round(I.cols*scales[i]/pChns.shrink)*pChns.shrink;

		//debug, while we dont have imResampleMex
		I1 = I;

		if(h1 == I.rows && w1 == I.cols)
			I1 = I;
		else //need to implement the imResampleMex functions
			//I1 = imResampleMex(I.data,h1,w1,1);
			;
		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
			convertedImage = I1; //is this correct?

		computedChannels[ccIndex] = computeSingleScaleChannelFeatures(I1);
		ccIndex++;

		//why is this here?
		//couldn't it be outside the for?
		//or is it really needed now that computedChannels=data?
		//if(i==isR(1)), nTypes=chns.nTypes; data=cell(nScales,nTypes); end

	}
		
	int is[computedScales/approximatedScales+1]; //needs a better name
	int isIndex = 0;
	bool isError = false;

	//if lambdas not specified compute image specific lambdas
	//chnsPyramid.m, line 154
	//computedScales is the substitute for nScales, might need to 	
	//set it in the computeSingleScaleChannelFeatures function
	if (computedScales>0 && approximatedScales>0 && !providedLambdas)
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
		double *f1 = (double*)memset(f1, 0, channelTypes*sizeof(double));;

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
	
	//compute image pyramid [approximated scales]
	//the next for statement is controlled by the isA array
	//h1, w1 and I1 were declared earlier
	//this will use imResampleMex too, so it will be erased for now
	for (int i=0; i<computedScales; i++)
	{
		h1 = round(I.rows*scales[i]/pChns.shrink);
		w1 = round(I.cols*scales[i]/pChns.shrink);
		
		double ratio[3];
		int iR = isN[i];
		ratio[0] = pow(scales[i]/scales[iR],-lambdas[0]);
		ratio[1] = pow(scales[i]/scales[iR],-lambdas[1]);
		ratio[2] = pow(scales[i]/scales[iR],-lambdas[2]);
		
	}
	
	//smooth channels, optionally pad and concatenate channels
	for (int i=0; i < computedScales; i++)
	{
		// this is kinda of what needs to happen
		// but the actual convolution needs to be reworked
		/*computedChannels[i].image = convTri(computedChannels[i].image, smooth);
		computedChannels[i].gradientMagnitude = convTri(computedChannels[i].gradientMagnitude, smooth);
		computedChannels[i].gradientHistogram = convTri(computedChannels[i].gradientHistogram, smooth);*/
	}

}

//translation of the chnsCompute.m file
Info Pyramid::computeSingleScaleChannelFeatures(cv::Mat I)
{
	cv::Mat gradOrientation;
	Info result;

	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);

	cv::imshow("image before conversion", I);
	cv::waitKey();				
	cv::destroyAllWindows();

	//compute color channels
	result.image = pChns.pColor.rgbConvert(I);

	cv::imshow("image after conversion, before convolution", result.image);
	cv::waitKey();				
	cv::destroyAllWindows();

	result.image = pChns.pColor.convolution(result.image, pChns.pColor.smooth, 1, CONV_TRI);

	if (result.image.cols < 2)
	{
		cv::imshow("image is empty after convolution", I);
		cv::waitKey();				
		cv::destroyAllWindows();
	}
	else
	{
		cv::imshow("image after convolution", result.image);
		cv::waitKey();				
		cv::destroyAllWindows();
	}

	if (pChns.pColor.enabled)
		result.colorCh = pChns.pColor;

	if (pChns.pGradHist.enabled)
	{
		cv::imshow("before mGradMag", I);
		cv::waitKey();				
		cv::destroyAllWindows();

		if (result.image.cols < 2)
		{
			cv::imshow("empty result.image", I);
			cv::waitKey();				
			cv::destroyAllWindows();
		}

		//i need to identify which is the color channel to know
		//how to represent it in integer, or change mGradMag
		std::vector<cv::Mat> tempResult = pChns.pGradMag.mGradMag(result.image,COLOR_CHANNEL);

		if (tempResult.size() > 1)
			result.gradientMagnitude = tempResult[0];
		if (tempResult.size() > 1)
			gradOrientation = tempResult[1];

		cv::imshow("after assignments", I);
		cv::waitKey();				
		cv::destroyAllWindows();

		//still need to understand this next part:
		if (pChns.pGradMag.normalizationRadius != 0)
		{
			float *S = (float*)(pChns.pColor.convolution(result.gradientMagnitude, pChns.pGradMag.normalizationRadius, 1, CONV_TRI)).data;
			float *M = (float*)result.gradientMagnitude.data;
			int h = result.gradientMagnitude.rows;
			int w = result.gradientMagnitude.cols;
			result.gradientMagnitude = pChns.pGradMag.gradMagNorm(M, S, h, w); 
		}
	}		
	else
	{
		if (pChns.pGradMag.enabled)
		{
			result.gradientMagnitude = (pChns.pGradMag.mGradMag(result.image, COLOR_CHANNEL))[0];			

			if (pChns.pGradMag.normalizationRadius != 0)
			{
				float *S = (float*)(pChns.pColor.convolution(result.gradientMagnitude, pChns.pGradMag.normalizationRadius, 1, CONV_TRI)).data;
				float *M = (float*)result.gradientMagnitude.data;
				int h = result.gradientMagnitude.rows;
				int w = result.gradientMagnitude.cols;
				result.gradientMagnitude = pChns.pGradMag.gradMagNorm(M, S, h, w);
			}
		}	
	}
	
	//compute gradient histogram channels
	if (pChns.pGradHist.enabled)
	{
		result.gradientHistogram = pChns.pGradHist.mGradHist(result.gradientMagnitude, gradOrientation, pChns.pGradMag.full);
	}	

	//for now, i wont add computation of custom channels
	
	cv::imshow("results of chnsCompute", I);
	cv::waitKey();				
	cv::destroyAllWindows();	

	return result;
}

// set each scale s such that max(abs(round(sz*s/shrink)*shrink-sz*s)) is minimized 
// without changing the smaller dim of sz (tricky algebra)
// getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);
void Pyramid::getScales(int h, int w, int shrink)
{
	int minSize, bgDim, smDim;
	double minSizeRatio;

	
	if (h!=0 && w!=0)
	{
		if (h <= w)
		{
			bgDim = w;
			smDim = h;
			minSizeRatio = h / minImgSize[0];
		}
		else
		{
			bgDim = h;
			smDim = w;
			minSizeRatio = w / minImgSize[1];
		}

		//nScales = floor(nPerOct*(nOctUp+log2(min(sz./minDs)))+1);
		computedScales = floor(scalesPerOctave*(upsampledOctaves+log2(minSizeRatio))+1);
		int s0, s1;
		double epsilon = std::numeric_limits<double>::epsilon();	

		scales = (double*)malloc(computedScales * sizeof(double));

		//this next for statement will substitute the following:
		//scales = 2.^(-(0:nScales-1)/nPerOct+nOctUp);
		for (int i=0; i < computedScales; i++)
		{
			scales[i]=pow(-(i/scalesPerOctave+upsampledOctaves),2);
			s0=(round(smDim*scales[i]/shrink)*shrink-.25*shrink)/smDim;
			s1=(round(smDim*scales[i]/shrink)*shrink+.25*shrink)/smDim;

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

			// should tempMaxScale be the size of ssIndex or is it 
			// related to computedScales, or both?
			double tempMaxScale[ssIndex];				
		
			//this is the max part of [~,x]=min(max(es0,es1)); 
			for (int k=0; k < ssIndex; k++)
				if (es0[k] > es1[k])
					tempMaxScale[k] = es0[k];
				else
					tempMaxScale[k] = es1[k];

			//this is the min part of [~,x]=min(max(es0,es1));
			double minScaleValue = tempMaxScale[0];
			int minScaleIndex = 0; 
			for (int j=1; j < ssIndex; j++)
				if (tempMaxScale[j] < minScaleValue)
				{
					minScaleValue = tempMaxScale[j];
					minScaleIndex = j;
				}
			
			//scales(i)=ss(x);
			scales[i] = ss[minScaleIndex];	
		}
		
		//in here, we would have a text to just keep the values of 
		//scales[i] which are different from their neighbours
		//for now, this will be suppressed
		
		scaleshw = (cv::Point*)malloc(computedScales * sizeof(cv::Point));

		for (int i=0; i<computedScales; i++)
		{
			scaleshw[i].x = round(w*scales[i]/shrink)*shrink/w;
			scaleshw[i].y = round(h*scales[i]/shrink)*shrink/h;
		}
	}
	else //error, height or width of the image are wrong
		;
}


