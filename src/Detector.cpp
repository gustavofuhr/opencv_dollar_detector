#include "Detector.h"

//reads the detector model from the xml model
void Detector::readDetectorModel(String fileName)
{
    FileStorage xml;
    FileNode currentRootNode;

    xml.open(fileName, FileStorage::READ);

	if (!xml.isOpened())
	{
		std::cerr << "Failed to open " << fileName << std::endl;
	}
	else
	{
        currentRootNode = xml["detector"]["opts"]["pPyramid"]["pChns"];
        opts.pPyramid.pChns.shrink = currentRootNode["shrink"];
        opts.pPyramid.pChns.pColor.enabled = currentRootNode["pColor"]["enabled"];
        opts.pPyramid.pChns.pColor.smooth = currentRootNode["pColor"]["smooth"];
        opts.pPyramid.pChns.pColor.colorSpaceType = (string)currentRootNode["pColor"]["colorSpace"];
        opts.pPyramid.pChns.pGradMag.enabled = currentRootNode["pColor"]["pGradMag"]["enabled"];
        opts.pPyramid.pChns.pGradMag.colorChannelIndex = currentRootNode["pColor"]["pGradMag"]["colorChn"];
        opts.pPyramid.pChns.pGradMag.normalizationRadius = currentRootNode["pColor"]["pGradMag"]["normRad"];
        opts.pPyramid.pChns.pGradMag.normalizationConstant = currentRootNode["pColor"]["pGradMag"]["normConst"];
        opts.pPyramid.pChns.pGradMag.full = currentRootNode["pColor"]["pGradMag"]["full"];
        opts.pPyramid.pChns.pGradHist.enable = currentRootNode["pColor"]["pGradHist"]["enabled"];
        opts.pPyramid.pChns.pGradHist.orientationChannels = currentRootNode["pColor"]["pGradHist"]["nOrients"];
        opts.pPyramid.pChns.pGradHist.useSoftBinning = currentRootNode["pColor"]["pGradHist"]["softBin"];
        opts.pPyramid.pChns.pGradHist.useHogNormalization = currentRootNode["pColor"]["pGradHist"]["useHog"];
        opts.pPyramid.pChns.pGradHist.clipHog = currentRootNode["pColor"]["pGradHist"]["clipHog"];
        opts.pPyramid.pChns.complete = currentRootNode["complete"];

        currentRootNode = xml["detector"]["opts"]["pPyramid"];
        opts.pPyramid.scalesPerOctave = currentRootNode["nPerOct"];
        opts.pPyramid.upsampledOctaves = currentRootNode["nOctUp"];
        opts.pPyramid.approximatedScales = currentRootNode["nApprox"];
        opts.pPyramid.lambdas[0] = currentRootNode["lambdas"][0];
        opts.pPyramid.lambdas[1] = currentRootNode["lambdas"][1];
        opts.pPyramid.lambdas[2] = currentRootNode["lambdas"][2];
        opts.pPyramid.pad[0] = currentRootNode["pad"][0];
        opts.pPyramid.pad[1] = currentRootNode["pad"][1];
        opts.pPyramid.minImgSize[0] = currentRootNode["minDs"][0];
        opts.pPyramid.minImgSize[1] = currentRootNode["minDs"][1];
        opts.pPyramid.smoothRadius = currentRootNode["smooth"];
        opts.pPyramid.concatenateChannels = currentRootNode["concat"];
        opts.pPyramid.completeInput = currentRootNode["complete"];

        currentRootNode = xml["detector"]["opts"];
        opts.modelDs[0] = currentRootNode["modelDs"][0];
        opts.modelDs[1] = currentRootNode["modelDs"][1];
        opts.modelDsPad[0] = currentRootNode["modelDsPad"][0];
        opts.modelDsPad[1] = currentRootNode["modelDsPad"][1];
        opts.stride = currentRootNode["stride"];
        opts.cascadeThreshold = currentRootNode["cascThr"];
        opts.cascadeCalibration = currentRootNode["cascCal"];
        opts.nWeak[0] = currentRootNode["nWeak"][0];
        opts.nWeak[1] = currentRootNode["nWeak"][1];
        opts.nWeak[2] = currentRootNode["nWeak"][2];
        opts.nWeak[3] = currentRootNode["nWeak"][3];
        opts.seed = currentRootNode["seed"];
        opts.name = (string)currentRootNode["name"];
        opts.posGtDir = (string)currentRootNode["posGtDir"];
        opts.posImgDir = (string)currentRootNode["posImgDir"];
        opts.negImgDir = (string)currentRootNode["negImgDir"];
        opts.nPos = currentRootNode["nPos"];
        opts.nNeg = currentRootNode["nNeg"];
        opts.nPerNeg = currentRootNode["nPerNeg"];
        opts.nAccNeg = currentRootNode["nAccNeg"];
        opts.winsSave = currentRootNode["winsSave"];

		  //there's a problem to read these big matrices because the numbers are not in scientific notation
        currentRootNode = xml["detector"]["clf"];
        currentRootNode["fids"] >> clf.fids;

        //all of the elements in fids are zero
        qDebug() << clf.fids.data[0];
        qDebug() << clf.fids.data[1];
        qDebug() << clf.fids.data[2];
        currentRootNode["thrs"] >> clf.thrs;

        //some more matrices would need to be read here...

        clf.treeDepth = currentRootNode["treeDepth"];

        currentRootNode = xml["detector"]["info"];
        info.colorCh.enabled = currentRootNode["colorCh"]["enabled"];
        info.colorCh.smooth = currentRootNode["colorCh"]["smooth"];
        info.colorCh.colorSpaceType = (string) currentRootNode["colorCh"]["colorSpace"];
        info.colorCh.nChannels = currentRootNode["colorCh"]["nChns"];
        info.colorCh.padWith = (string) currentRootNode["colorCh"]["padWith"];

        info.gradMag.enabled = currentRootNode["pGradMag"]["enabled"];
        info.gradMag.colorChannelIndex = currentRootNode["pGradMag"]["colorChn"];
        info.gradMag.normalizationRadius = currentRootNode["pGradMag"]["normRad"];
        info.gradMag.normalizationConstant = currentRootNode["pGradMag"]["normConst"];
        info.gradMag.full = currentRootNode["pGradMag"]["full"];
        info.gradMag.nChannels = currentRootNode["pGradMag"]["nChns"];
        info.gradMag.padWith = (string) currentRootNode["pGradMag"]["padWith"];

        info.gradHist.enable = currentRootNode["pGradHist"]["enabled"];
        info.gradHist.orientationChannels = currentRootNode["pGradHist"]["nOrients"];
        info.gradHist.useSoftBinning = currentRootNode["pGradHist"]["softBin"];
        info.gradHist.useHogNormalization = currentRootNode["pGradHist"]["useHog"];
        info.gradHist.clipHog = currentRootNode["pGradHist"]["clipHog"];
        info.gradHist.nChannels = currentRootNode["pGradHist"]["nChns"];
        info.gradHist.padWith = (string) currentRootNode["pGradHist"]["padWith"];
    }
}

//this procedure was just copied verbatim
inline void getChild(float *chns1, uint32_t *cids, uint32_t *fids,
    float *thrs, uint32_t offset, uint32_t &k0, uint32_t &k)
{
	float ftr = chns1[cids[fids[k]]];
	k = (ftr<thrs[k]) ? 1 : 2;
	k0 = k += k0 * 2; k += offset;
}

//bb = acfDetect1(P.data{i},Ds{j}.clf,shrink,modelDsPad(1),modelDsPad(2),opts.stride,opts.cascThr);
BoundingBox* Detector::acfDetect(Mat image)
{
    BoundingBox bb;
    int detectorLength = 100; //what is the length of the detector model?

    //teste para ver se o conteudo da imagem eh char, se for aplica a funcao imreadf

    //calcular a piramide da imagem

    //criar uma matriz de bounding boxes bbs[nmero de escalas da pirmida][nmero de elementos do detector]

    //this loop was just copied from the original file (except some comments)
    for (int i = 0; i < this->opts.pPyramid.computedScales; i++)
    {
        for (int j = 0; j < detectorLength; j++)
        {
            float* chns = (float*)image.data;
            const int shrink = this->opts.pPyramid.pChns.shrink;
            const int modelHt = this->opts.modelDsPad[0];
            const int modelWd = this->opts.modelDsPad[1];
            const int stride = this->opts.stride;
            const float cascThr = this->opts.cascadeThreshold;

            float *thrs = (float*) this->clf.thrs.data;
            float *hs = (float*) this->clf.hs;
            uint32_t *fids = (uint32_t*) this->clf.fids.data;
            uint32_t *child = (uint32_t*) this->clf.child;
            const int treeDepth = this->clf.treeDepth;

            const int height = image.rows;
            const int width = image.cols;
            const int nChns = image.channels();
            //nTreeNodes  number of elements in the first dimension of fids
            const int nTreeNodes = this->clf.fids.rows;
            //nTrees size of the second dimension of fids
            const int nTrees = this->clf.fids.cols;
            const int height1 = (int)ceil(float(height*shrink - modelHt + 1 / stride));
            const int width1 = (int)ceil(float(width*shrink - modelWd + 1 / stride));

            //construct cids array
            int nFtrs = modelHt / shrink*modelWd / shrink*this->info.colorCh.nChannels;
            uint32_t *cids = new uint32_t[nFtrs];
            int m = 0;
            for (int z = 0; z<this->info.colorCh.nChannels; z++)
            for (int c = 0; c<modelWd / shrink; c++)
            for (int r = 0; r<modelHt / shrink; r++)
                cids[m++] = z*width*height + c*height + r;


            // apply classifier to each patch
            vector<int> rs, cs; vector<float> hs1;
            for (int c = 0; c<width1; c++) for (int r = 0; r<height1; r++) {
                float h = 0, *chns1 = chns + (r*stride / shrink) + (c*stride / shrink)*height;
                if (treeDepth == 1) {
                    // specialized case for treeDepth==1
                    for (int t = 0; t < nTrees; t++) {
                        uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
                        this->getChild(chns1, cids, fids, thrs, offset, k0, k);
                        h += hs[k]; if (h <= cascThr) break;
                    }
                }
                else if (treeDepth == 2) {
                    // specialized case for treeDepth==2
                    for (int t = 0; t < nTrees; t++) {
                        uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
                        this->getChild(chns1, cids, fids, thrs, offset, k0, k);
                        this->getChild(chns1, cids, fids, thrs, offset, k0, k);
                        h += hs[k]; if (h <= cascThr) break;
                    }
                }
                else if (treeDepth>2) {
                    // specialized case for treeDepth>2
                    for (int t = 0; t < nTrees; t++) {
                        uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
                        for (int i = 0; i<treeDepth; i++)
                            this->getChild(chns1, cids, fids, thrs, offset, k0, k);
                        h += hs[k]; if (h <= cascThr) break;
                    }
                }
                else {
                    // general case (variable tree depth)
                    for (int t = 0; t < nTrees; t++) {
                        uint32_t offset = t*nTreeNodes, k = offset, k0 = k;
                        while (child[k]) {
                            float ftr = chns1[cids[fids[k]]];
                            k = (ftr<thrs[k]) ? 1 : 0;
                            k0 = k = child[k0] - k + offset;
                        }
                        h += hs[k]; if (h <= cascThr) break;
                    }
                }
                if (h>cascThr) { cs.push_back(c); rs.push_back(r); hs1.push_back(h); }
            }
            delete [] cids;
            m=cs.size();

            // convert to bbs
            double *bbs;
            for( int i=0; i<m; i++ )
            {
                bbs[i+0*m]=cs[i]*stride; bbs[i+2*m]=modelWd;
                bbs[i+1*m]=rs[i]*stride; bbs[i+3*m]=modelHt;
                bbs[i+4*m]=hs1[i];
            }
        }
    }
    return NULL;
}
