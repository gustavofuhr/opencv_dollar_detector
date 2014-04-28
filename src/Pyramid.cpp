#include "Pyramid.h"

//translation of the chnsPyramid.m file
Pyramid Pyramid::computeMultiScaleChannelFeaturePyramid(Mat I)
{
	Mat convertedImage;

	//if we are to allow incomplete Pyramids, we need to set some values to default.
	//for now, it wont be implemented. (lines 115-128 of chnsPyramid.m)

	//convert I to appropriate color space (or simply normalize)
	convertedImage = pChns.pColor.rgbConvert(I);
	pChns.pColor.colorSpaceType = "orig";

	//get scales at which to compute features and list of real/approx scales
	//i dont think scalehw is ever used
	//[scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	double* scales = getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);

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
	Mat I1;
	int ccIndex=0;

	for (int i=0; i < computedScales; i = i + approximatedScales + 1)
	{
		h1 = round(I.rows*scales[i]/pChns.shrink)*pChns.shrink;
		w1 = round(I.cols*scales[i]/pChns.shrink)*pChns.shrink;

		if(h1 == I.rows && w1 == I.cols)
			I1 = I;
		else //calling the mex.cpp file didnt work, need to implement those functions
			//I1 = imResampleMex(I.data,h1,w1,1);
	
		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
			convertedImage = I1;

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

	//suppliedLambdas must determine if the lambdas where supplied
	//this value needs to be set properly in the future
	bool suppliedLambdas = false;

	//if lambdas not specified compute image specific lambdas
	//chnsPyramid.m, line 154
	//computedScales is the substitute for nScales, might need to 	
	//set it in the computeSingleScaleChannelFeatures function
	if (computedScales>0 && approximatedScales>0 && !suppliedLambdas)
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
		double *f0 = (double*)memset(f0, 0, channelTypes*sizeof(double));		
		double *f1 = (double*)memset(f1, 0, channelTypes*sizeof(double));;

		//only two of the channel types seem to be considered here
		//since chnsCompute computes the histogram channel after
		//the other two, its ommited
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
	for (int i=0; i<computedScales; i++)
	{
		h1 = round(I.rows*scales[i]/pChns.shrink);
		w1 = round(I.cols*scales[i]/pChns.shrink);
		
		//to know which elements of scales are accessed, i need to
		//finish get scales part, so this will be completed later
		double ratio=0;
		for (int j=0; j<channelTypes; j++);
			//ratio = (scales[i]/scales)
	}
	
	//smooth channels, optionally pad and concatenate channels
	for (int i=0; i < computedScales*channelTypes; i++)
	{
		;
	}

	//create output struct
	//maybe i'll return a new Pyramid or i'll just update the current
}

//translation of the chnsCompute.m file
Info Pyramid::computeSingleScaleChannelFeatures(Mat I)
{
	Mat gradOrientation;
	Info result;

	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);

	//compute color channels
	result.image = pChns.pColor.rgbConvert(I);
	result.image = pChns.pColor.convolution(result.image, pChns.pColor.smooth, 1, CONV_TRI);
	if (pChns.pColor.enabled)
		result.colorCh = pChns.pColor;

	//compute gradient magnitude channel
	if (pChns.pGradHist.enabled)
	{
		//i need to identify which is the color channel to know
		//how to represent it in integer, or change mGradMag
		/*Mat *tempResult = pChns.pGradMag.mGradMag(result.image,result.colorCh);
		result.gradientMagnitude = tempResult[0];
		gradOrientation = tempResult[1];*/

		//still need to understand this next part:
		if (pChns.pGradMag.normalizationRadius != 0)
		{
			//where does this M come from?
			//this part depends on the last one i took out
			/*float* S = (pChns.pColor.convolution(M, pChns.pGradMag.normalizationRadius, 1, CONV_TRI)).data;
			result.gradientMagnitude = pChns.pGradMag.gradMagNorm(result.gradientMagnitude,S);*/ 
		}
	}		
	else
	{
		if (pChns.pGradMag.enabled)
		{
			//came problem with channel selection as before
			//but maybe this is always done in the color channel
			//result.gradientMagnitude = (pChns.pGradMag.mGradMag(result.image,pChns.pColor))[0];			

			if (pChns.pGradMag.normalizationRadius != 0)
			{
				//also depends on previous problems
				/*float* S = (pChns.pColor.convolution(M, pChns.pGradMag.normalizationRadius, 1, CONV_TRI)).data;
			result.gradientMagnitude = pChns.pGradMag.gradMagNorm(result.gradientMagnitude,S);*/ 
			}
		}	
	}

	//compute gradient histogram channels
	if (pChns.pGradHist.enabled)
	{
		//decided to send full as one
		//this needs to be checked later, because of the
		//gradOrientation Matrix
		//result.gradientHistogram = pChns.pGradHist.mGradHist(result.gradientMagnitude, gradOrientation, 1);
	}	

	//for now, i wont add computation of custom channels
	
	return result;
}

//set each scale s such that max(abs(round(sz*s/shrink)*shrink-sz*s)) is minimized 
//without changing the smaller dim of sz (tricky algebra)
double* Pyramid::getScales(int h, int w, int shrink)
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

		computedScales = floor(scalesPerOctave*(upsampledOctaves+log2(minSizeRatio))+1);
		double scales[computedScales];
		double tempMaxScale[computedScales];
		int s0, s1;
		double epsilon = std::numeric_limits<double>::epsilon();

		//this next for statement will substitute the following:
		//scales = 2.^(-(0:nScales-1)/nPerOct+nOctUp);
		for (int i=0; i < computedScales; i++)
		{
			scales[i]=pow(-(i/scalesPerOctave+upsampledOctaves),2);
			s0=(round(smDim*scales[i]/shrink)*shrink-.25*shrink)/smDim;
			s1=(round(smDim*scales[i]/shrink)*shrink+.25*shrink)/smDim;

			//what follows will substitute ss=(0:.01:1-epsilon())*(s1-s0)+s0;
			double ss[(int)round(1-epsilon/0.01)], es0[(int)round(1-epsilon/0.01)], es1[(int)round(1-epsilon/0.01)];						
			int ssIndex = 0;			
			for (int j=0; j < 1-epsilon; j = j + 0.01)
			{
				ss[ssIndex] = j*(s1-s0)+s0;
				es0[ssIndex]=smDim*ss[ssIndex]; 
				//es0[ssIndex]=abs(es0-round(es0/shrink)*shrink);
				es1[ssIndex]=bgDim*ss[ssIndex]; 
				//es1[ssIndex]=abs(es1-round(es1/shrink)*shrink);
				ssIndex++;
			}

			//this is the max part of [~,x]=min(max(es0,es1)); 
			for (int j=0; j < ssIndex; j++)
				if (es0[j] > es1[j])
					tempMaxScale[j] = es0[j];
				else
					tempMaxScale[j] = es1[j];

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
		
		//still have to decide whether i need to return scaleshw
		/*scaleshw = [round(h*scales/shrink)*shrink/h;
  round(w*scales/shrink)*shrink/w]';*/
		
		//both scales and scalehw are returned, but i dont see scaleshw being used ever
		return scales;
	}
	else //error, height or width are wrong
		;
}


