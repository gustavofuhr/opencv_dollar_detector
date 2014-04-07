#include "Pyramid.h"

//translation of the chnsPyramid.m file
Pyramid computeMultiScaleChannelFeaturePyramid(Mat I)
{
	Mat convertedImage;
	//if we are to allow incomplete Pyramids, we need to set what has to value to the default.
	//for now, it wont be implemented. (lines 115-128 of chnsPyramid.m)

	//convert I to appropriate color space (or simply normalize)
    convertedImage = this->pChns.pColor.rgbConvert(I);
	this->pChns.pColor.colorSpaceType = "orig";

	//get scales at which to compute features and list of real/approx scales
	/*
	[scales,scaleshw]=getScales(nPerOct,nOctUp,minDs,shrink,sz);
	*/
}

//translation of the chnsCompute.m file
Info Pyramid::computeSingleScaleChannelFeaturePyramid(Mat I)
{
	Mat gradOrientation;
	Info result;

	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);


	//compute color channels
	result.image = this->pChns.pColor.rgbConvert(I);
	result.image = this->pChns.pColor.convConst(colorChResult, CONV_TRI);
	if (this->pChns.pColor.enabled)
		result.colorCh = this->pChns.pColor;

	//compute gradient magnitude channel
	if (this->pChns.pGradHist.enabled)
	{
		Mat tempResult = this->pChns.pGradMag.mGradMag(result.image,result.colorChn,full);
		result.gradientMagnitude = tempResult[0];
		gradOrientation = tempResult[1];
		//still need to understand this next part:
		/*if (this->pChns.pGradMag.normalizationRadius != 0)
		{
			S = convTri(M, normRad);
			//gradientMex('gradientMagNorm',M,S,normConst); 
		}*/
	}		
	else
	{
		if (this->pChns.pGradMag.enabled)
		{
			result.gradMagnitude = (this->pChns.pGradMag.mGradMag(result.image,p.colorChn,full))[0];			
			/*if (this->pChns.pGradMag.normalizationRadius != 0)
			{
				S = convTri(M, normRad);
				//gradientMex('gradientMagNorm',M,S,normConst); 
			}*/
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

Mat Pyramid::TriangleFilterConvolution(Mat I, int r, int s, int nomex)
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
                    /*for (int i=t; i < result.rows-s+t+1; i = i + s)
                        for (int j=t; j < result.cols-s+t+1; j = j + s)
                            temp*/
                    //result = result(t:s:end-s+t,t:s:end-s+t,:);
                }
            }
        }
    }
    return result;
}

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
		//this next for is made to substitute the following:
		//scales = 2.^(-(0:nScales-1)/nPerOct+nOctUp);
		for (int i=0; i < nScales; i++)
		{
			scales[i]=pow(-(i/scalesPerOctave+upsampledOctaves),2);
			s0=(round(smDim*scales[i]/shrink)*shrink-.25*shrink)./smDim;
			s1=(round(smDim*scales[i]/shrink)*shrink+.25*shrink)./smDim;
			//the folowing will substitute ss=(0:.01:1-epsilon())*(s1-s0)+s0;
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
			//gotta test to see exactly what this is doing.
			/*[~,x]=min(max(es0,es1)); 
			scales(i)=ss(x);*/
		}
		/*
		kp=[scales(1:end-1)~=scales(2:end) true]; scales=scales(kp);
		scaleshw = [round(sz(1)*scales/shrink)*shrink/sz(1);
  round(sz(2)*scales/shrink)*shrink/sz(2)]';
		*/
	}
	else //error
		;
}


