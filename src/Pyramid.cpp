#include "Pyramid.h"
#include <highgui.h>

void Pyramid::readScalesFromXML(cv::FileNode pyramid)
{
	cv::Mat scale01;
	pyramid["scale01"] >> scale01;
	float* scale1F = cvImage2floatArray(scale01, 1);
	print_100_elements(scale1F, scale01.rows, "scale 1 read from xml file");

	cv::Mat scale02;
	pyramid["scale02"] >> scale02;
	float* scale2F = cvImage2floatArray(scale02, 1);
	print_100_elements(scale2F, scale02.rows, "scale 2 read from xml file");

	cv::Mat scale03;
	pyramid["scale03"] >> scale03;
	float* scale3F = cvImage2floatArray(scale03, 1);
	print_100_elements(scale3F, scale03.rows, "scale 3 read from xml file");

	std::cin.get();
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
	
	padSize = pyramidNode["padSize"];
	pad = (int*)malloc(padSize*sizeof(int));
	for (int i=0; i < padSize; i++)
		pad[i] = pyramidNode["pad"][i];

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
	convertedImage = rgbConvert(I, pChns.pColor.colorSpaceType);
	pChns.pColor.colorSpaceType = ORIG;

	/*
	// debug: loads image converted inside matlab
	cv::FileStorage xml;
	xml.open("../opencv_dollar_detector/luvImage.xml", cv::FileStorage::READ);
	if (!xml.isOpened())
	{
		std::cerr << "Failed to open luvImage.xml" << std::endl;
		std::cin.get();
	}
	cv::Mat ch1, ch2, ch3;
	std::vector<cv::Mat> channels;
	xml["luvImage"]["ch1"] >> ch1;
	xml["luvImage"]["ch2"] >> ch2;
	xml["luvImage"]["ch3"] >> ch3;

	channels.push_back(ch1);
	channels.push_back(ch2);
	channels.push_back(ch3);

	cv::merge(channels, convertedImage);

	//cv::imshow("image loaded from xml file", convertedImage);
	//cv::waitKey();
	// debug */

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
	getScales(convertedImage.rows, convertedImage.cols, pChns.shrink);

	computedChannels = new Info[computedScales];
	
	int h1, w1;
	cv::Mat I1;
	int numberOfRealScales;
	int i;

	// compute image pyramid [real scales]
	for (i=0; i < computedScales; i = i+approximatedScales+1)
	{
		// sz=[size(I,1) size(I,2)];
		// sz1=round(sz*s/shrink)*shrink;
		h1 = round(I.rows*scales[i]/pChns.shrink)*pChns.shrink;
		w1 = round(I.cols*scales[i]/pChns.shrink)*pChns.shrink;

		// debug
		// std::cout << std::endl << "now computing scales[" << i << "] = " << scales[i] << ", shrink=" << pChns.shrink << ", h1=" << h1 << ", w1=" << w1 << std::endl;

		if (h1 == I.rows && w1 == I.cols)
			I1 = convertedImage;
		else // I1=imResampleMex(I,sz1(1),sz1(2),1);
			I1 = resample(convertedImage,convertedImage.rows,convertedImage.cols,h1,w1,1.0, 3);

		// if(s==.5 && (nApprox>0 || nPerOct==1)), I=I1;
		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
			convertedImage = I1; 

		computedChannels[i] = computeSingleScaleChannelFeatures(I1);
	}
	numberOfRealScales = i;

	// debug
	std::cout << "after real scales" << std::endl;

	int is[computedScales/(approximatedScales+1)]; //needs a better name
	int isIndex = 0;
	bool isError = false;
	// if lambdas not specified compute image specific lambdas
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
			// sz1=round(sz*scales(i)/shrink);
			h1 = round(I.rows*scales[i]/pChns.shrink);
			w1 = round(I.cols*scales[i]/pChns.shrink);
			
			double ratio[3];
			int iR = 0;

			// real scale changes in i=5 and i=13 to iR=8 and iR=16  
			for (int j=0; j < numberOfRealScales-1; j++)
			{
				if (i > floor((j*(approximatedScales+1)+(j+1)*(approximatedScales+1))/2))
					iR = (j+1)*(approximatedScales+1);
			}

			// debug
			// std::cout << std::endl << "approximating, i=" << i << ", scales[iR]=" << scales[iR] << ", iR=" << iR << ", h1=" << h1 << ", w1=" << w1 << std::endl;

			ratio[0] = pow(scales[i]/scales[iR],-lambdas[0]);
			computedChannels[i].image = resample(computedChannels[iR].image, computedChannels[iR].image.rows, computedChannels[iR].image.cols, h1, w1, ratio[0], 3);
		
			ratio[1] = pow(scales[i]/scales[iR],-lambdas[1]);
			computedChannels[i].gradientMagnitude = resample(computedChannels[iR].gradientMagnitude, computedChannels[iR].gradientMagnitude.rows, computedChannels[iR].gradientMagnitude.cols, h1, w1, ratio[1], 1);
			
			ratio[2] = pow(scales[i]/scales[iR],-lambdas[2]);

			for (int k=0; k < pChns.pGradHist.nChannels; k++)
			{
				int h = computedChannels[iR].gradientHistogram[k].rows;
				int w = computedChannels[iR].gradientHistogram[k].cols;
				computedChannels[i].gradientHistogram.push_back(resample(computedChannels[iR].gradientHistogram[k], h, w, h1, w1, ratio[2], 1));

				/*
				// debug 
				std::cout << "approx gradHist for scale[" << i << "], rows=" << h << ", cols=" 
				<< w << ", type=" << computedChannels[i].gradientHistogram[k].type() << ", ratio=" << ratio[2] << std::endl;
				// debug */
			}

			/*
			// debug: test results of image and gradMag channel approximations
			cv::imshow("iR image", computedChannels[iR].image);
			cv::imshow("iR gradMag", computedChannels[iR].gradientMagnitude);

			std::cout << "approx image for scale[" << i << "], rows=" << computedChannels[i].image.rows << ", cols=" 
			<< computedChannels[i].image.cols << ", type=" << computedChannels[i].image.type() << ", ratio=" << ratio[0] << std::endl;
			cv::imshow("approx image", computedChannels[i].image);

			std::cout << "approx gradMag for scale[" << i << "], rows=" << computedChannels[i].gradientMagnitude.rows << ", cols=" 
			<< computedChannels[i].gradientMagnitude.cols << ", type=" << computedChannels[i].gradientMagnitude.type() << ", ratio=" << ratio[1] << std::endl;
			cv::imshow("approx gradMag", computedChannels[i].gradientMagnitude);
			// debug */

			/*
			// debug: test results of gradHist approximations
			cv::imshow("iR gradHist0", computedChannels[iR].gradientHistogram[0]);
			cv::imshow("app gradHist0", computedChannels[i].gradientHistogram[0]);
			cv::imshow("iR gradHist1", computedChannels[iR].gradientHistogram[1]);
			cv::imshow("app gradHist1", computedChannels[i].gradientHistogram[1]);
			cv::imshow("iR gradHist2", computedChannels[iR].gradientHistogram[2]);
			cv::imshow("app gradHist2", computedChannels[i].gradientHistogram[2]);
			cv::imshow("iR gradHist3", computedChannels[iR].gradientHistogram[3]);
			cv::imshow("app gradHist3", computedChannels[i].gradientHistogram[3]);
			cv::imshow("iR gradHist4", computedChannels[iR].gradientHistogram[4]);
			cv::imshow("app gradHist4", computedChannels[i].gradientHistogram[4]);
			cv::imshow("iR gradHist5", computedChannels[iR].gradientHistogram[5]);
			cv::imshow("app gradHist5", computedChannels[i].gradientHistogram[5]);
			// */

			// debug
			// cv::waitKey();
			
			// debug
			// std::cout << "end of i=" << i << ", scales[iR]=" << scales[iR] << ", iR=" << iR << ", h1=" << h1 << ", w1=" << w1 << std::endl;
		}
	}

	// debug
	std::cout << "after approx scales" << std::endl;

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


		computedChannels[i].image = convolution(computedChannels[i].image, 3, pChns.pColor.smoothingRadius, 1, CONV_TRI);		
		computedChannels[i].gradientMagnitude = convolution(computedChannels[i].gradientMagnitude, 1, pChns.pColor.smoothingRadius, 1, CONV_TRI);

		for (int j=0; j < pChns.pGradHist.nChannels; j++)
			computedChannels[i].gradientHistogram[j] = convolution(computedChannels[i].gradientHistogram[j], 1, pChns.pColor.smoothingRadius, 1, CONV_TRI);


		// debug: test values inside computedChannels
		// color channels are still correct here
		// testFeatures(computedChannels[0], "after convolution");
		//std::cin.get();
		// debug*/

		/*
		// debug: images after convolution
		cv::imshow("conv image", computedChannels[i].image);
		cv::imshow("conv gradMag", computedChannels[i].gradientMagnitude);
		cv::imshow("conv gradHist", computedChannels[i].gradientHistogram[1]);
		cv::waitKey();
		// debug */	

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

		// debug: test values inside computedChannels
		// testFeatures(computedChannels[0], "after padding");
		//std::cin.get();
		// debug*/

		/*
		// debug 
		std::cout << "scale[" << i << "], imgSize=(" << computedChannels[i].image.rows << "," << computedChannels[i].image.cols << 
			"), magSize=(" << computedChannels[i].gradientMagnitude.rows << "," << computedChannels[i].gradientMagnitude.cols << 
			"), hisSize=(" << computedChannels[i].gradientHistogram[0].rows << "," << computedChannels[i].gradientHistogram[0].cols << ")" << std::endl;
		// debug */

	}
}

//translation of the chnsCompute.m file
Info Pyramid::computeSingleScaleChannelFeatures(cv::Mat I)
{
	cv::Mat gradOrientation;
	Info result;

	/*
	[h,w,~]=size(I); cr=mod([h w],shrink);
	if(any(cr)), h=h-cr(1); w=w-cr(2); I=I(1:h,1:w,:); end
	*/
	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);

	
	// debug
	float *If = cvImage2floatArray(I, 3);
	print_100_elements(If, I.rows, "image inside chnsCompute before rgbConvert");
	// debug */

	// compute color channels
	result.image = rgbConvert(I, pChns.pColor.colorSpaceType);
	result.image = convolution(result.image, 3, pChns.pColor.smoothingRadius, 1, CONV_TRI);

	
	// debug
	float *If2 = cvImage2floatArray(result.image, 3);
	print_100_elements(If2, I.rows, "image after convolution");
	// debug */

	// h=h/shrink; w=w/shrink;
	height = height/pChns.shrink;
	width = width/pChns.shrink;

	if (pChns.pGradHist.enabled)
	{
		std::vector<cv::Mat> tempResult = pChns.pGradMag.mGradMag(result.image,0);

		if (tempResult.size() > 0)
			result.gradientMagnitude = tempResult[0];
		if (tempResult.size() > 1)
			gradOrientation = tempResult[1];

		if (pChns.pGradMag.normalizationRadius != 0)
		{
			cv::Mat convRes = convolution(result.gradientMagnitude, 1, pChns.pGradMag.normalizationRadius, 1, CONV_TRI);

			float *S = cvImage2floatArray(convRes, 1);
			float *M = cvImage2floatArray(result.gradientMagnitude, 1);

			int h = result.gradientMagnitude.rows;
			int w = result.gradientMagnitude.cols;

			// normalization constant is read inside the procedure
			pChns.pGradMag.gradMagNorm(M, S, h, w);

			
			// debug: after normalization, S matrix is correct but gradMag is wrong (but equal to compiled mex result, which is the wrong one)
			print_100_elements(M, h, "gradMag after normalization");
			print_100_elements(S, h, "S matrix");
			//std::cin.get();
			// debug */

			result.gradientMagnitude = floatArray2cvImage(M, h, w, 1); 		}	
	}		
	else
	{
		if (pChns.pGradMag.enabled)
		{
			result.gradientMagnitude = (pChns.pGradMag.mGradMag(result.image, 0))[0];			

			if (pChns.pGradMag.normalizationRadius != 0)
			{
				float *S = cvImage2floatArray(convolution(result.gradientMagnitude, 1, pChns.pGradMag.normalizationRadius, 1, CONV_TRI), 1);
				float *M = cvImage2floatArray(result.gradientMagnitude, 1);

				int h = result.gradientMagnitude.rows;
				int w = result.gradientMagnitude.cols;

				pChns.pGradMag.gradMagNorm(M, S, h, w);

				result.gradientMagnitude = floatArray2cvImage(M, h, w, 1); // only one channel
			}
		}	
	}

	//compute gradient histogram channels
	if (pChns.pGradHist.enabled)
	{
		result.gradientHistogram = pChns.pGradHist.mGradHist(result.gradientMagnitude, gradOrientation, pChns.pGradMag.full);

		for (int i=0; i < pChns.pGradHist.nChannels; i++)
			result.gradientHistogram[i] = resample(result.gradientHistogram[i], result.gradientHistogram[i].rows, result.gradientHistogram[i].cols, height, width, 1.0, 1);
	}	

	// if(p.enabled), chns=addChn(chns,I,nm,p,'replicate',h,w); end
	if (pChns.pColor.enabled)
	{
		if (result.image.rows!=height || result.image.cols!=width)
		{
			// data=imResampleMex(data,h,w,1);
			result.image = resample(result.image, result.image.rows, result.image.cols, height, width, 1.0, pChns.pColor.nChannels);
		}
		// result.colorCh = pChns.pColor;
	}

	// if(p.enabled), chns=addChn(chns,M,nm,p,0,h,w); end
	if (pChns.pGradMag.enabled)
		result.gradientMagnitude = resample(result.gradientMagnitude, result.gradientMagnitude.rows, result.gradientMagnitude.cols, height, width, 1.0, pChns.pGradMag.nChannels);

	/*
	// debug: print results for every channel
	cv::imshow("0 - testing rgbConvert", result.image);
	cv::imshow("1 - testing convolution", result.image);
	cv::imshow("2 - testing gradientMagnitude", result.gradientMagnitude);
	cv::imshow("3 - testing gradientOrientation", gradOrientation);
	cv::imshow("4 - testing gradientHistogram", result.gradientHistogram[0]);
	cv::waitKey();		
	// debug */

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

		// nScales = floor(nPerOct*(nOctUp+log2(min(sz./minDs)))+1);
		// nScales = floor(8*(0+log2(min(576/100, 720/41)))+1);
		// in the current tests, returns 21 inside MATLAB. In here, returns 21 too
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
		scales = (double*)malloc(computedScales*sizeof(double));
		scales[0] = tempScales[0];
		int scalesIndex=1;
		for (int i=1; i < computedScales; i++)
			if(tempScales[i] != tempScales[i-1])
			{
				scales[scalesIndex++] = tempScales[i];
			}

		// this updates the value of computedScales, since some of them have been suppressed
		computedScales = scalesIndex;
		
		scales_h = (double*)malloc(computedScales * sizeof(double));
		scales_w = (double*)malloc(computedScales * sizeof(double));

		// scaleshw = 	[round(sz(1)*scales/shrink)*shrink/sz(1);
  		//				round(sz(2)*scales/shrink)*shrink/sz(2)]';
		for (int i=0; i<computedScales; i++)
		{
			scales_w[i] = round(w*scales[i]/shrink)*shrink/w;
			scales_h[i] = round(h*scales[i]/shrink)*shrink/h;
		}
	}
	else //error, height or width of the image are wrong
		std::cout << " # getScales error: both height and width need to be greater than 0!";
}


