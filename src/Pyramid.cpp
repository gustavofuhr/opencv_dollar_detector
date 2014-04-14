#include "Pyramid.h"

//translation of the chnsPyramid.m file
Pyramid computeMultiScaleChannelFeaturePyramid(Mat I)
{
	Mat convertedImage;
	Info computedChannels;

	//if we are to allow incomplete Pyramids, we need to set some values to default.
	//for now, it wont be implemented. (lines 115-128 of chnsPyramid.m)

	//convert I to appropriate color space (or simply normalize)
    convertedImage = this->pChns.pColor.rgbConvert(I);
	this->pChns.pColor.colorSpaceType = "orig";

	//get scales at which to compute features and list of real/approx scales
	/*
	[scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	//next, there's a bunch of weirdness happening:	
	nScales=length(scales); if(1), isR=1; else isR=1+nOctUp*nPerOct; end
isR=isR:nApprox+1:nScales; isA=1:nScales; isA(isR)=[];
j=[0 floor((isR(1:end-1)+isR(2:end))/2) nScales];
isN=1:nScales; for i=1:length(isR), isN(j(i)+1:j(i+1))=isR(i); end
nTypes=0; data=cell(nScales,nTypes); info=struct([]);
	*/
	
	//compute image pyramid, from chnsPyramid.m, line 144
	//there is a for statement where the index i can only assume values of the array isR
	//I guess isR might be useless, if we do it like this instead:
	int h1, w1;
	Mat I1;
	for (int i=0; i < nScales; i = i + approximatedScales + 1)
	{
		h1 = round(I.rows*scales[i]/shrink)*shrink;
		w1 = round(I.cols*scales[i]/shrink)*shrink;

		if(h1 == I.rows && w1 == I.cols)
			I1 = I;
		else
			//I1 = imResampleMex(I,h1,w1,1);
	
		if (scales[i] == 0.5 && (approximatedScales>0 || scalesPerOctave == 1))
			convertedImage = I1;

		computedChannels = computeSingleScaleChannelFeatures(I1);

		//why is this here?
		//couldn't it be outside the for?
		//if(i==isR(1)), nTypes=chns.nTypes; data=cell(nScales,nTypes); end

		//data is a structure initilized as data=cell(nScales,nTypes);
		//in here data receives the content of chns.data, which...
		//is the return of the chnsCompute function
		data[i] = computedChannels;
	}

	int is[nScales/approximatedScales+1]; //needs a better name
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
		for (int j=1+upsampledOctaves*scalesPerOctave; j < nScales; j = j+approximatedScales+1)
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
		if (isIndex > 1)
			//is=is(2:3)

		//channelTypes is the substitute for the nTypes value
		double f0[channelTypes] = {0};
		double f1[channelTypes] = {0};

		//these matrices are receiving the value of chns.data
		//that would be the computed matrices of the Info type
		//returned on the chnsCompute.m function
		//the question is: why aren't the three channels read?
		Mat d0 = computedChannels.image;
		Mat d1 = computedChannels.gradMag;

		for (int j=0; j<channelTypes; j++)
		{
			//f0 and f1 seem to be the average value of each cell 
			//of d0 and d1
			double sum=0;
			int i;
			for (i=0; i<computedChannels.image.dims; i++)
				sum = sum + computedChannels.image.data[j][i];			
			f0[j] = sum/i;

			sum=0;
			int k;
			for (k=0; k<computedChannels.gradMag.dims; k++)
			{
				sum = sum + computedChannels.gradMag.data[j][k];	
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
		h1 = round(I.rows*scales[i]/shrink);
		w1 = round(I.cols*scales[i]/shrink);
		
		//to know which elements of scales are accessed, i need to
		//finished get scales part, so this will be completed later
		double ratio=0;
		for (int j=0; j<channelTypes; j++)
			ratio = (scales[i]/scales)
	}
	
	//smooth channels, optionally pad and concatenate channels

	//create output struct
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
	result.image = this->pChns.pColor.rgbConvert(I);
	result.image = this->pChns.pColor.convolution(result.image, this->pChns.pColor.smooth, 1, CONV_TRI);
	if (this->pChns.pColor.enabled)
		result.colorCh = this->pChns.pColor;

	//compute gradient magnitude channel
	if (this->pChns.pGradHist.enabled)
	{
		Mat *tempResult = this->pChns.pGradMag.mGradMag(result.image,result.colorChn,full);
		result.gradientMagnitude = tempResult[0];
		gradOrientation = tempResult[1];

		//still need to understand this next part:
		if (this->pChns.pGradMag.normalizationRadius != 0)
		{
			float* S = (this->pChns.pColor.convolution(M, this->pChns.pGradMag.normalizationRadius, 1, CONV_TRI)).data;
			result.gradientMagnitude = this->pChns.pGradMag.gradMagNorm(result.gradientMagnitude,S); 
		}
	}		
	else
	{
		if (this->pChns.pGradMag.enabled)
		{
			result.gradMagnitude = (this->pChns.pGradMag.mGradMag(result.image,p.colorChn,full))[0];			

			if (this->pChns.pGradMag.normalizationRadius != 0)
			{
				float* S = (this->pChns.pColor.convolution(M, this->pChns.pGradMag.normalizationRadius, 1, CONV_TRI)).data;
			result.gradientMagnitude = this->pChns.pGradMag.gradMagNorm(result.gradientMagnitude,S); 
			}
		}	
	}

	//compute gradient histogram channels
	if (this->pChns.pGradHist.enabled)
	{
		result.gradientHistogram = this->pChns.pGradHist.mGradHist(result.gradientMagnitude, gradOrientation, full);
	}	

	//for now, i wont add computation of custom channels
	
	return result;
}

/*Mat Pyramid::TriangleFilterConvolution(Mat I, int r, int s, int nomex)
{
    Mat result = I;
    if(!I.empty() && !(r==0 && s==1))
    {
        int m = min(I.rows, I.cols);
        if(m>=4 && m>2*r+1)
        {
            if(nomex==0)
            {
                if(r>0 && r<=1 && s<=2)
                    result = this->pChns.pColor.convConst(I,CONV_TRI1);
                else
                    result = this->pChns.pColor.convConst(I,CONV_TRI);
            }
            else
            {
                double f[3];
                if (r <= 1)
                {
                    double p=12/r/(r+2)-2;
                    f[0] = 1/(2+p);
                    f[2] = f[0];
                    f[3] = p / (2+p);
                    r=1;
                }
                else
                {
                    for (int i=1; i <= r+1; i++)
                        f[i] = i/(pow(r+1,2));
                    for (int i=2*r+1; i > r+1; i--)
                        f[i] = i/(pow(r+1,2));
                }
                //result = padarray(I,[r r],"symmetric","both");
                copyMakeBorder(I, result, 0, 0, r, r, BORDER_REPLICATE);
                //result = convn(convn(result,f,"valid"),f',"valid");
                if (s>1)
                {
                    int t=floor(s/2)+1;
                    Mat temp;
                    for (int i=t; i < result.rows-s+t+1; i = i + s)
                        for (int j=t; j < result.cols-s+t+1; j = j + s)
                            temp
                    //result = result(t:s:end-s+t,t:s:end-s+t,:);
                }
            }
        }
    }
    return result;
}*/

//set each scale s such that max(abs(round(sz*s/shrink)*shrink-sz*s)) is minimized 
//without changing the smaller dim of sz (tricky algebra)
void Pyramid::getScales(int h, int w, int shrink)
{
	int nScales; 
	int minSize, bgDim, smDim;
	
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

		nScales = floor(scalesPerOctave*(upsampledOctaves+log2(minSizeRatio))+1);
		double scales[nScales];
		double tempMaxScale[nscales];

		//this next for statement will substitute the following:
		//scales = 2.^(-(0:nScales-1)/nPerOct+nOctUp);
		for (int i=0; i < nScales; i++)
		{
			scales[i]=pow(-(i/scalesPerOctave+upsampledOctaves),2);
			s0=(round(smDim*scales[i]/shrink)*shrink-.25*shrink)./smDim;
			s1=(round(smDim*scales[i]/shrink)*shrink+.25*shrink)./smDim;

			//what follows will substitute ss=(0:.01:1-epsilon())*(s1-s0)+s0;
			double ss[round(1-epsilon()/0.01)], es0[round(1-epsilon()/0.01)], es1[round(1-epsilon()/0.01)];						
			int ssIndex = 0;			
			for (int j=0; j < 1-epsilon(); j = j + 0.01)
			{
				ss[ssIndex] = j*(s1-s0)+s0;
				es0[ssIndex]=smDim*ss[ssIndex]; 
				es0[ssIndex]=abs(es0-round(es0/shrink)*shrink);
				es1[ssIndex]=bgDim*ss[ssIndex]; 
				es1[ssIndex]=abs(es1-round(es1/shrink)*shrink);
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
		//the kp variable needs to be set up
		//i'm yet to understand the purpose of this:		
		//kp=[scales(1:end-1)~=scales(2:end) true]; 
		
		//next, the scales array is changed using kp:		
		//scales=scales(kp);
		
		scaleshw = [round(h*scales/shrink)*shrink/h;
  round(w*scales/shrink)*shrink/w]';
		
		//both scales and scalehw are returned
	}
	else //error, height or width are wrong
		;
}


