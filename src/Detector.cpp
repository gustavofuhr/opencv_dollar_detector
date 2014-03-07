#include "Detector.h"

//Procedimento para leitura do arquivo xml que contém o modelo de detector
void Detector::readDetectorModel(String fileName)
{
    FileStorage xml;
    FileNode currentRootNode;

    //instrução para abertura do arquivo
    xml.open(fileName, FileStorage::READ);

	if (!xml.isOpened())
	{
		std::cerr << "Failed to open " << fileName << std::endl;
	}
	else
	{
        //Série de atribuições dos valores ao objeto a ser criado
        currentRootNode = xml["detector"]["opts"]["pPyramid"]["pChns"];
        opts.pPyramid.pChns.shrink = currentRootNode["shrink"];
        opts.pPyramid.pChns.pColor.enabled = currentRootNode["pColor"]["enabled"];
        opts.pPyramid.pChns.pColor.smooth = currentRootNode["pColor"]["smooth"];
        opts.pPyramid.pChns.pColor.colorSpace = (string)currentRootNode["pColor"]["colorSpace"];
        opts.pPyramid.pChns.pGradMag.enabled = currentRootNode["pColor"]["pGradMag"]["enabled"];
        opts.pPyramid.pChns.pGradMag.colorChn = currentRootNode["pColor"]["pGradMag"]["colorChn"];
        opts.pPyramid.pChns.pGradMag.normRad = currentRootNode["pColor"]["pGradMag"]["normRad"];
        opts.pPyramid.pChns.pGradMag.normConst = currentRootNode["pColor"]["pGradMag"]["normConst"];
        opts.pPyramid.pChns.pGradMag.full = currentRootNode["pColor"]["pGradMag"]["full"];
        opts.pPyramid.pChns.pGradHist.enabled = currentRootNode["pColor"]["pGradHist"]["enabled"];
        opts.pPyramid.pChns.pGradHist.nOrients = currentRootNode["pColor"]["pGradHist"]["nOrients"];
        opts.pPyramid.pChns.pGradHist.softBin = currentRootNode["pColor"]["pGradHist"]["softBin"];
        opts.pPyramid.pChns.pGradHist.useHog = currentRootNode["pColor"]["pGradHist"]["useHog"];
        opts.pPyramid.pChns.pGradHist.clipHog = currentRootNode["pColor"]["pGradHist"]["clipHog"];
        opts.pPyramid.pChns.complete = currentRootNode["complete"];

        currentRootNode = xml["detector"]["opts"]["pPyramid"];
        opts.pPyramid.nPerOct = currentRootNode["nPerOct"];
        opts.pPyramid.nOctUp = currentRootNode["nOctUp"];
        opts.pPyramid.nApprox = currentRootNode["nApprox"];
        opts.pPyramid.lambdas[0] = currentRootNode["lambdas"][0];
        opts.pPyramid.lambdas[1] = currentRootNode["lambdas"][1];
        opts.pPyramid.lambdas[2] = currentRootNode["lambdas"][2];
        opts.pPyramid.pad[0] = currentRootNode["pad"][0];
        opts.pPyramid.pad[1] = currentRootNode["pad"][1];
        opts.pPyramid.minDs[0] = currentRootNode["minDs"][0];
        opts.pPyramid.minDs[1] = currentRootNode["minDs"][1];
        opts.pPyramid.smooth = currentRootNode["smooth"];
        opts.pPyramid.concat = currentRootNode["concat"];
        opts.pPyramid.complete = currentRootNode["complete"];

        currentRootNode = xml["detector"]["opts"];
        opts.modelDs[0] = currentRootNode["modelDs"][0];
        opts.modelDs[1] = currentRootNode["modelDs"][1];
        opts.modelDsPad[0] = currentRootNode["modelDsPad"][0];
        opts.modelDsPad[1] = currentRootNode["modelDsPad"][1];
        opts.stride = currentRootNode["stride"];
        opts.cascThr = currentRootNode["cascThr"];
        opts.cascCal = currentRootNode["cascCal"];
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

        //problema para ler essas matrizes gigantes porque os números não estão em notação científica
        //tentei fazer a inicialização abaixo para testar se o problema é o tipo esperado, não funcionou
        //clf.fids = Mat(7,2048, CV_32S, 0.0);
        currentRootNode = xml["detector"]["clf"];
        currentRootNode["fids"] >> clf.fids;

        //todos os elementos lidos ficam zerados
        qDebug() << clf.fids.data[0];
        qDebug() << clf.fids.data[1];
        qDebug() << clf.fids.data[2];
        currentRootNode["thrs"] >> clf.thrs;

        //mais algumas matrizes precisam ser lidas aqui...

        clf.treeDepth = currentRootNode["treeDepth"];

        currentRootNode = xml["detector"]["info"];
        info.colorCh.enabled = currentRootNode["colorCh"]["enabled"];
        info.colorCh.smooth = currentRootNode["colorCh"]["smooth"];
        info.colorCh.colorSpace = (string) currentRootNode["colorCh"]["colorSpace"];
        info.colorCh.nChannels = currentRootNode["colorCh"]["nChns"];
        info.colorCh.padWith = (string) currentRootNode["colorCh"]["padWith"];

        info.gradMag.enabled = currentRootNode["pGradMag"]["enabled"];
        info.gradMag.colorChn = currentRootNode["pGradMag"]["colorChn"];
        info.gradMag.normRad = currentRootNode["pGradMag"]["normRad"];
        info.gradMag.normConst = currentRootNode["pGradMag"]["normConst"];
        info.gradMag.full = currentRootNode["pGradMag"]["full"];
        info.gradMag.nChannels = currentRootNode["pGradMag"]["nChns"];
        info.gradMag.padWith = (string) currentRootNode["pGradMag"]["padWith"];

        info.gradHist.enabled = currentRootNode["pGradHist"]["enabled"];
        info.gradHist.nOrients = currentRootNode["pGradHist"]["nOrients"];
        info.gradHist.softBin = currentRootNode["pGradHist"]["softBin"];
        info.gradHist.useHog = currentRootNode["pGradHist"]["useHog"];
        info.gradHist.clipHog = currentRootNode["pGradHist"]["clipHog"];
        info.gradHist.nChannels = currentRootNode["pGradHist"]["nChns"];
        info.gradHist.padWith = (string) currentRootNode["pGradHist"]["padWith"];
    }
}

//Procedimento copiado apenas
inline void getChild(float *chns1, uint32_t *cids, uint32_t *fids,
	float *thrs, uint32_t offset, uint32_t &k0, uint32_t &k)
{
	float ftr = chns1[cids[fids[k]]];
	k = (ftr<thrs[k]) ? 1 : 2;
	k0 = k += k0 * 2; k += offset;
}
