#include "Pyramid.h"

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

	// convert I to appropriate color space (or simply normalize)
	// I=rgbConvert(I,cs); pChns.pColor.colorSpace='orig';
	convertedImage = pChns.pColor.rgbConvert(I);
	pChns.pColor.colorSpaceType = "orig";

	std::cout << "nPerOct = " << scalesPerOctave << ", nOctUp = " << upsampledOctaves << ", minDs = (" << minImgSize[0] << "," << minImgSize[1] << 
	"), shrink = " << pChns.shrink << ", sz = (" << convertedImage.rows << "," << convertedImage.cols << ")" <<std::endl;

	/*
	nPerOct = 8, nOctUp = 0, minDs = (100,41), shrink = 4, sz = (576,720)
	scales[0]=0
	scales[1]=3
	*/

	// get scales at which to compute features and list of real/approx scales
	// [scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);

	computedChannels = new Info[computedScales];
	
	//compute image pyramid, from chnsPyramid.m, line 144
	//there is a for statement where the index i can only assume values of the array isR
	//I guess isR might be useless, if we do it like this instead:
	int h1, w1;
	cv::Mat I1;
	int ccIndex=0;
	int numberOfRealScales;
	int i;

	// real scales are multiples of approximatedScales+1 !
	// approximated scales are the others
	// what are upsampled scales?
	for (i=0; i < computedScales; i = i+approximatedScales+1)
	{

		// sz=[size(I,1) size(I,2)];
		// sz1=round(sz*s/shrink)*shrink;
		h1 = round(I.rows*scales[i]/pChns.shrink)*pChns.shrink;
		w1 = round(I.cols*scales[i]/pChns.shrink)*pChns.shrink;
		std::cout << std::endl << "now computing scales[" << i << "] = " << scales[i] << ", shrink=" << pChns.shrink << std::endl;

		if (h1 == I.rows && w1 == I.cols)
			I1 = convertedImage;
		else // I1=imResampleMex(I,sz1(1),sz1(2),1);
		{
			std::cout << "inside chnsPyramid, before resample" << std::endl;
			I1 = resample(convertedImage,h1,w1,1.0, 3);
			std::cout << "inside chnsPyramid, after resample" << std::endl;
		}

		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
			convertedImage = I1; //is this correct?

		std::cout << "inside chnsPyramid, before chnsCompute" << std::endl;
		computedChannels[ccIndex] = computeSingleScaleChannelFeatures(I1);
		ccIndex++;
		std::cout << "inside chnsPyramid, after chnsCompute" << std::endl;

		//why is this here?
		//couldn't it be outside the for?
		//or is it really needed now that computedChannels=data?
		//if(i==isR(1)), nTypes=chns.nTypes; data=cell(nScales,nTypes); end
	}
	numberOfRealScales = i;
		
	int is[computedScales/approximatedScales+1]; //needs a better name
	int isIndex = 0;
	bool isError = false;

	//if lambdas not specified compute image specific lambdas
	//chnsPyramid.m, line 154
	//computedScales is the substitute for nScales
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
	
	// isR=isR:nApprox+1:nScales; 
	// isR = 1, 9, 17 in Matlab, in here, it becomes isR = 0, 8, 16
	// isA=1:nScales; isA(isR)=[];
	// isA holds the indices that are not in isR
	// isA = 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20
	// j as the points at which we change the real scale used to approximate the current one
	// j=[0 floor((isR(1:end-1)+isR(2:end))/2) nScales];
	// j=[0,5,13,21] 
	// isN holds the scales to be used as a basis for approximations
	// isN=1:nScales; for i=1:length(isR), isN(j(i)+1:j(i+1))=isR(i); end
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
	for (int i=0; i<computedScales; i++)
	{
		if (i % (approximatedScales+1) != 0)
		{
			h1 = round(I.rows*scales[i]/pChns.shrink);
			w1 = round(I.cols*scales[i]/pChns.shrink);
			
			double ratio[3];
			int iR = 0;

			for (int j=0; j < numberOfRealScales-1; j++)
			{
				if (i > floor((j*(approximatedScales+1)+(j+1)*(approximatedScales+1))/2))
					iR = (j+1)*(approximatedScales+1);
			}

			std::cout << "approximating, i=" << i << ", scales[iR]=" << scales[iR] << ", iR=" << iR << ", h1=" << h1 << ", w1=" << w1 << std::endl;
			ratio[0] = pow(scales[i]/scales[iR],-lambdas[0]);
			computedChannels[i].image = resample(computedChannels[iR].image, h1, w1, ratio[0], 3);
			ratio[1] = pow(scales[i]/scales[iR],-lambdas[1]);
			computedChannels[i].gradientMagnitude = resample(computedChannels[iR].gradientMagnitude, h1, w1, ratio[1], 1);
			ratio[2] = pow(scales[i]/scales[iR],-lambdas[2]);
			computedChannels[i].gradientHistogram = resample(computedChannels[iR].gradientHistogram, h1, w1, ratio[2], 1);
		}
	}

	//smooth channels, optionally pad and concatenate channels
	for (int i=0; i < computedScales; i++)
	{
		computedChannels[i].image = convolution(computedChannels[i].image, 3, pChns.pColor.smoothingRadius, 1, CONV_TRI);
		computedChannels[i].gradientMagnitude = convolution(computedChannels[i].gradientMagnitude, 1, pChns.pColor.smoothingRadius, 1, CONV_TRI);
		computedChannels[i].gradientHistogram = convolution(computedChannels[i].gradientHistogram, 1, pChns.pColor.smoothingRadius, 1, CONV_TRI);
	}

	std::cout << "end of chnsPyramid" << std::endl;
}

//translation of the chnsCompute.m file
Info Pyramid::computeSingleScaleChannelFeatures(cv::Mat I)
{
	cv::Mat gradOrientation;
	Info result;

	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);
	result.image = cv::Mat(height, width, I.type());

	cv::imshow("testing source image", I);

	//compute color channels
	result.image = pChns.pColor.rgbConvert(I);

	cv::imshow("testing rgbConvert", result.image);

	// I = convTri(I,p.smooth);
	// J = convConst('convTri',I,r,s);
	// A = (float*) mxGetData(prhs[1]);
	// p = (float) mxGetScalar(prhs[2]);
 	// r = (int) mxGetScalar(prhs[2]);
  	// s = (int) mxGetScalar(prhs[3]);
  	// ms[0]=ns[0]/s; ms[1]=ns[1]/s; ms[2]=d;
  	// B = (float*) mxMalloc(ms[0]*ms[1]*d*sizeof(float));
  	// nDims = mxGetNumberOfDimensions(prhs[1]);
  	// ns = (int*) mxGetDimensions(prhs[1]);
  	// d = (nDims == 3) ? ns[2] : 1;
  	// convTri( A, B, ns[0], ns[1], d, r, s );
	result.image = convolution(result.image, 3, pChns.pColor.smoothingRadius, 1, CONV_TRI);

	cv::imshow("testing convolution", result.image);
	cv::waitKey();

	if (pChns.pColor.enabled)
		result.colorCh = pChns.pColor;

	if (pChns.pGradHist.enabled)
	{
		// debug
		std::cout << "before mGradMag, rows=" << result.image.rows << ", cols=" << result.image.cols << std::endl;

		std::vector<cv::Mat> tempResult = pChns.pGradMag.mGradMag(result.image,COLOR_CHANNEL);

		// debug
		std::cout << "chnsCompute, after mGradMag" << std::endl;

		if (tempResult.size() > 0)
			result.gradientMagnitude = tempResult[0];
		if (tempResult.size() > 1)
			gradOrientation = tempResult[1];

		if (pChns.pGradMag.normalizationRadius != 0)
		{

			// debug
			std::cout << "inside chnsCompute, inside normalization if" << std::endl;

			// this convolution happens in a single channel matrix, causing a segmentation fault!
			cv::Mat convRes = convolution(result.gradientMagnitude, 1, pChns.pGradMag.normalizationRadius, 1, CONV_TRI);

			// debug
			std::cout << "inside chnsCompute, after convolution" << std::endl;

			float *S = cvMat2floatArray(convRes, 1);

			// debug
			std::cout << "inside chnsCompute, after S conversion" << std::endl;

			float *M = cvMat2floatArray(result.gradientMagnitude, 1);

			// debug
			std::cout << "inside chnsCompute, after M conversion" << std::endl;

			int h = result.gradientMagnitude.rows;
			int w = result.gradientMagnitude.cols;

			// debug
			std::cout << "inside chnsCompute, before gradMagNorm" << std::endl;

			pChns.pGradMag.gradMagNorm(M, S, h, w);

			// debug
			std::cout << "inside chnsCompute, after gradMagNorm" << std::endl;


			// debug
			std::cout << "inside chnsCompute, before conversion 0" << std::endl;

			result.gradientMagnitude = floatArray2cvMat(M, h, w, 1); // only one channel

			// debug
			std::cout << "inside chnsCompute, after conversion 0" << std::endl;
		}
	}		
	else
	{
		if (pChns.pGradMag.enabled)
		{
			result.gradientMagnitude = (pChns.pGradMag.mGradMag(result.image, COLOR_CHANNEL))[0];			

			if (pChns.pGradMag.normalizationRadius != 0)
			{
				float *S = cvMat2floatArray(convolution(result.gradientMagnitude, 1, pChns.pGradMag.normalizationRadius, 1, CONV_TRI), 1);
				float *M = cvMat2floatArray(result.gradientMagnitude, 1);
				int h = result.gradientMagnitude.rows;
				int w = result.gradientMagnitude.cols;
				pChns.pGradMag.gradMagNorm(M, S, h, w);
				result.gradientMagnitude = floatArray2cvMat(M, h, w, 1); // only one channel
			}
		}	
	}

	//compute gradient histogram channels
	if (pChns.pGradHist.enabled)
	{

		// debug
		std::cout << "inside chnsCompute, before mGradHist" << std::endl;

		result.gradientHistogram = pChns.pGradHist.mGradHist(result.gradientMagnitude, gradOrientation, pChns.pGradMag.full);

		// debug
		std::cout << "inside chnsCompute, after mGradHist" << std::endl;
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

	/*
	nPerOct = 8, nOctUp = 0, minDs = (100,41), shrink = 4, sz = (576,720)
	scales[0]=0
	scales[1]=3
	*/
	
	if (h!=0 && w!=0)
	{
		if (h/minImgSize[0] < w/minImgSize[1])
			minSizeRatio = float(h) / minImgSize[0];
		else
			minSizeRatio = float(w) / minImgSize[1];

		//nScales = floor(nPerOct*(nOctUp+log2(min(sz./minDs)))+1);
		// nScales = floor(8*(0+log2(min(576/100, 720/41)))+1);
		// in the current tests, returns 21 inside MATLAB. In here, returns 21 too
		computedScales = floor(scalesPerOctave*(upsampledOctaves+log2(minSizeRatio))+1);

		// prints for testing computedScales
		// std::cout << "computedScales = " << computedScales << std::endl;
		// std::cout << "h = " << h << ", w = " << w << ", minImgSize[0] = " << minImgSize[0] << ", minImgSize[1] = " << minImgSize[1] << ", minSizeRatio = " << minSizeRatio  << std::endl;
		

		double s0, s1;
		double epsilon = std::numeric_limits<double>::epsilon();	

		//scales = (double*)malloc(computedScales * sizeof(double));
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
			// this print is now returning the correct values!
			// std::cout << "tempScales[" << i << "] = " << tempScales[i] << std::endl;

			// s0=(round(d0*s/shrink)*shrink-.25*shrink)./d0;
			// s1=(round(d0*s/shrink)*shrink+.25*shrink)./d0;
			s0=(round(smDim*tempScales[i]/shrink)*shrink-.25*shrink)/smDim;
			s1=(round(smDim*tempScales[i]/shrink)*shrink+.25*shrink)/smDim;
			// prints for testing s0 and s1 now returning correct values
			// std::cout << "s0 = " << s0 << ", s1 = " << s1 << std::endl;

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
			//scales[i] = ss[minScaleIndex];
			tempScales[i] = ss[minScaleIndex];	
		}
		
		//just keep the values of scales[i] which are different from their neighbours
		scales = (double*)malloc(computedScales*sizeof(double));
		scales[0] = tempScales[0];
		int scalesIndex=1;
		for (int i=1; i < computedScales; i++)
			if(tempScales[i] != tempScales[i-1])
			{
				scales[scalesIndex++] = tempScales[i];
			}

		computedScales = scalesIndex;
		
		scaleshw = (cv::Point*)malloc(computedScales * sizeof(cv::Point));

		for (int i=0; i<computedScales; i++)
		{
			scaleshw[i].x = round(w*scales[i]/shrink)*shrink/w;
			scaleshw[i].y = round(h*scales[i]/shrink)*shrink/h;
		}
	}
	else //error, height or width of the image are wrong
		std::cout << " # getScales error: both height and width need to be greater than 0!";
}


