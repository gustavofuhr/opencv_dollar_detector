#include "Pyramid.h"

//translation of the chnsPyramid.m file
Pyramid computeMultiScaleChannelFeaturePyramid(Mat I)
{
    ;
}

//translation of the chnsCompute.m file
//is the addChn function really needed?
Pyramid Pyramid::computeSingleScaleChannelFeaturePyramid(Mat I)
{
	Mat colorChResult, gradMagResult, graHistResult, gradOrientation;

	//crop I so it becomes divisible by shrink
	int height = I.rows - (I.rows % pChns.shrink);
	int width =  I.cols - (I.cols % pChns.shrink);


	//compute color channels
	colorChResult = this->pChns.pColor.rgbConvert(I);
	colorChResult = this->pChns.pColor.convConst(colorChResult, CONV_TRI);
	if (this->pChns.pColor.enabled)
		//chns = addChn(chns,I,nm,p,"replicate",h,w);

	//compute gradient magnitude channel
	if (this->pChns.pGradHist.enabled)
	{
		Mat tempResult = this->pChns.pGradMag.mGradMag(I,p.colorChn,full);
		gradMagResult 	= tempResult[0];
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
			gradMagResult = (this->pChns.pGradMag.mGradMag(I,p.colorChn,full))[0];			
			/*if (this->pChns.pGradMag.normalizationRadius != 0)
			{
				S = convTri(M, normRad);
				//gradientMex('gradientMagNorm',M,S,normConst); 
			}*/
		}	
	}
	if (this->pChns.pGradMag.enabled)
		//chns=addChn(chns,gradMagResult,nm,p,0,h,w);

	//compute gradient histogram channels
	if (this->pChns.pGradHist.enabled)
	{
		gradHistResult = this->pChns.pGradHist.mGradHist(gradMagResult, gradOrientation, full);
		//chns=addChn(chns,gradHistResult,nm,pChns.pGradHist,0,h,w);
	}	

	//for now, i wont add computation of custom channels
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
