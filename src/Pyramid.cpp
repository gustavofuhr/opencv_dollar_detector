#include "Pyramid.h"

void Pyramid::debugWindow (string name, double value)
{
	Mat debugImg = imread("../opencv_dollar_detector/frame0254.png");
	char str[200];
	sprintf(str,"%f",value);
	putText(debugImg, str, Point2f(100,100), FONT_HERSHEY_PLAIN, 2,  Scalar(0,0,255,255));
	imshow(name,debugImg);
	waitKey();
	destroyWindow(name);
}

void Pyramid::readPyramid(FileNode pyramidNode)
{
	pChns.readChannelFeatures(pyramidNode["pChns"]);

	scalesPerOctave = pyramidNode["nPerOct"];
	upsampledOctaves = pyramidNode["nOctUp"];
	approximatedScales = pyramidNode["nApprox"];
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
void Pyramid::computeMultiScaleChannelFeaturePyramid(Mat I)
{
	Mat convertedImage;

	//if we are to allow incomplete Pyramids, we need to set some values to default.
	//for now, it wont be implemented. (lines 115-128 of chnsPyramid.m)

	//convert I to appropriate color space (or simply normalize)
	convertedImage = pChns.pColor.rgbConvert(I);
	pChns.pColor.colorSpaceType = "orig";

	// this in here is just for debug purposes
	debugWindow("conv",0.0);

	// get scales at which to compute features and list of real/approx scales
	// [scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);

	// debug
	debugWindow("after scales",0.0);
	
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

		//debug
		I1 = I;

		if(h1 == I.rows && w1 == I.cols)
			I1 = I;
		else //need to implement the imResampleMex functions
			//I1 = imResampleMex(I.data,h1,w1,1);
			;
		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
			convertedImage = I1;

		//debug
		debugWindow("before chnsCompute",0.0);

		computedChannels[ccIndex] = computeSingleScaleChannelFeatures(I1);

		//debug
		debugWindow("after chnsCompute",0.0);

		ccIndex++;

		//why is this here?
		//couldn't it be outside the for?
		//or is it really needed now that computedChannels=data?
		//if(i==isR(1)), nTypes=chns.nTypes; data=cell(nScales,nTypes); end

	}
	
	// debug
	debugWindow("after computedChannels loop",0.0);
	
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

}

//translation of the chnsCompute.m file
Info Pyramid::computeSingleScaleChannelFeatures(Mat I)
{
	Mat gradOrientation;
	Info result;

	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);

		// debug
	debugWindow("before rgbConvert",0.0);

	//compute color channels
	result.image = pChns.pColor.rgbConvert(I);
	result.image = pChns.pColor.convolution(result.image, pChns.pColor.smooth, 1, CONV_TRI);
	if (pChns.pColor.enabled)
		result.colorCh = pChns.pColor;

	// debug
	debugWindow("after rgbConvert",0.0);

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
	
	debugWindow("after big if",0.0);
	
	//compute gradient histogram channels
	if (pChns.pGradHist.enabled)
	{
		//decided to send full as one
		//this needs to be checked later, because of the
		//gradOrientation Matrix
		//result.gradientHistogram = pChns.pGradHist.mGradHist(result.gradientMagnitude, gradOrientation, 1);
	}	

	//for now, i wont add computation of custom channels

	debugWindow("before return",0.0);
	
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
		
		// debug
		debugWindow("before scale loop",0.0);

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

		// debug
		debugWindow("after scale loop",0.0);
		
		//in here, we would have a text to just keep the values of 
		//scales[i] which are different from their neighbours
		//for now, this will be suppressed
		
		scaleshw = (Point*)malloc(computedScales * sizeof(Point));

		for (int i=0; i<computedScales; i++)
		{
			scaleshw[i].x = round(w*scales[i]/shrink)*shrink/w;
			scaleshw[i].y = round(h*scales[i]/shrink)*shrink/h;
		}
		// debug
		debugWindow("after scalehw",0.0);
	}
	else //error, height or width of the image are wrong
		;
}


