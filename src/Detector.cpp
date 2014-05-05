#include "Detector.h"

//i dont know if its gonna be needed but this is start
void Detector::exportDetectorModel(String fileName)
{
	FileStorage xml;
	
	xml.open(fileName, FileStorage::WRITE);

	xml << "opts" << "{";
		xml << "pPyramid" << "{";
			xml << "pChns" << "{";
				xml << "shrink" << opts.pPyramid.pChns.shrink;
				xml << "pColor" << "{";
					xml << "enabled" << opts.pPyramid.pChns.pColor.enabled;
				xml << "}";
			xml << "}";
		xml << "}";
	xml << "stride" << opts.stride;
	xml << "}";
	
	//xml << "clf" << this->clf;

	xml.release();
}

//reads the detector model from the xml model
//for now, it must be like this since the current model
//was not written by this program this will change after we are 
//set on a class structure
void Detector::importDetectorModel(String fileName)
{
	FileStorage xml;

	xml.open(fileName, FileStorage::READ);

	if (!xml.isOpened())
	{
		std::cerr << "Failed to open " << fileName << std::endl;
	}
	else
	{
		opts.readOptions(xml["detector"]["opts"]);

		clf.readClf(xml["detector"]["clf"]);		

		xml.release();
	}
}

//this procedure was just copied verbatim
void Detector::getChild(float *chns1, uint32_t *cids, uint32_t *fids,
    float *thrs, uint32_t offset, uint32_t &k0, uint32_t &k)
{
	float ftr = chns1[cids[fids[k]]];
	k = (ftr<thrs[k]) ? 1 : 2;
	k0 = k += k0 * 2; k += offset;
}

//bb = acfDetect1(P.data{i},Ds{j}.clf,shrink,modelDsPad(1),modelDsPad(2),opts.stride,opts.cascThr);
BoundingBox** Detector::acfDetect(Mat image)
{
	//teste para ver se o conteudo da imagem eh char, se for aplica a funcao imreadf 

	//compute feature pyramid
	opts.pPyramid.computeMultiScaleChannelFeaturePyramid(image);

	// creates a bounding box matrix
	// needs another dimension, one for the number of scales
	// another for the number of detections in each scale
	BoundingBox* detections[opts.pPyramid.computedScales];

	//this became a simple loop because we will apply just one detector here, 
	//to apply multiple detector models you need to create multiple Detector objects. 
	for (int i = 0; i < opts.pPyramid.computedScales; i++)
	{
		float* chns = (float*)image.data;
		const int shrink = opts.pPyramid.pChns.shrink;
		const int modelHt = opts.modelDsPad[0];
		const int modelWd = opts.modelDsPad[1];
		const int stride = opts.stride;
		const float cascThr = opts.cascadeThreshold;

		float *thrs = (float*) clf.thrs.data;
		float *hs = (float*) clf.hs.data;
		uint32_t *fids = (uint32_t*) clf.fids.data;
		uint32_t *child = (uint32_t*) clf.child.data;
		const int treeDepth = clf.treeDepth;

		const int height = opts.pPyramid.channelTypes;
		const int width = opts.pPyramid.computedScales;
		const int nChns = 1;

		//nTreeNodes size of the first dimension of fids
		const int nTreeNodes = clf.fids.rows;
		//nTrees size of the second dimension of fids
		const int nTrees = clf.fids.cols;
		const int height1 = (int)ceil(float(height*shrink-modelHt+1/stride));
		const int width1 = (int)ceil(float(width*shrink-modelWd+1/stride));

		//The number of color channels
		int nChannels = opts.pPyramid.pChns.pColor.nChannels;

		//construct cids array
		int nFtrs = modelHt / shrink*modelWd / shrink*nChannels;
		uint32_t *cids = new uint32_t[nFtrs];
		int m = 0;
		for (int z = 0; z<nChannels; z++)
			for (int c = 0; c<modelWd / shrink; c++)
				for (int r = 0; r<modelHt / shrink; r++)
					cids[m++] = z*width*height + c*height + r;

		// apply classifier to each patch
		vector<int> rs, cs; vector<float> hs1;
		for (int c = 0; c<width1; c++) 
			for (int r = 0; r<height1; r++) 
			{
				float h = 0, *chns1 = chns + (r*stride / shrink) + (c*stride / shrink)*height;
				if (treeDepth == 1) 
				{
					// specialized case for treeDepth==1
					for (int t = 0; t < nTrees; t++) 
					{
						uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
						getChild(chns1, cids, fids, thrs, offset, k0, k);
						h += hs[k]; if (h <= cascThr) break;
					}
				}	
				else if (treeDepth == 2) {
					// specialized case for treeDepth==2
					for (int t = 0; t < nTrees; t++) {
						uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
						getChild(chns1, cids, fids, thrs, offset, k0, k);
						getChild(chns1, cids, fids, thrs, offset, k0, k);
						h += hs[k]; if (h <= cascThr) break;
					}
				}
				else if (treeDepth>2) 
				{
					// specialized case for treeDepth>2
					for (int t = 0; t < nTrees; t++) 
					{
						uint32_t offset = t*nTreeNodes, k = offset, k0 = 0;
						for (int i = 0; i<treeDepth; i++)
							getChild(chns1, cids, fids, thrs, offset, k0, k);
						h += hs[k]; if (h <= cascThr) break;
					}
				}
				else 
				{
					// general case (variable tree depth)
					for (int t = 0; t < nTrees; t++) 
					{
						uint32_t offset = t*nTreeNodes, k = offset, k0 = k;
						while (child[k]) 
						{
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

			double shift[2];
			shift[0] = (modelHt-opts.modelDs[0])/2-opts.pPyramid.pad[0];
			shift[1] = (modelWd-opts.modelDs[1])/2-opts.pPyramid.pad[1];

			BoundingBox bb[m];
			for( int j=0; j<m; j++ )
			{
				bb[j].firstPoint.x = cs[j]*stride;
				bb[j].firstPoint.x = (bb[j].firstPoint.x+shift[1])/opts.pPyramid.scaleshw[j].x;
				bb[j].firstPoint.y = rs[j]*stride;
				bb[j].firstPoint.y = (bb[j].firstPoint.y+shift[0])/opts.pPyramid.scaleshw[j].y;
				bb[j].height = modelHt/opts.pPyramid.scales[j];
				bb[j].width = modelWd/opts.pPyramid.scales[j];
				bb[j].score = hs1[j];
			}

			// does this next part need to be translated?
			/*
			if(separate), bb(:,6)=j; end; 
			*/

			//detections[i] = new BoundingBox();
			//detections[i] = bb;
	} // */
 	// last part before returning the bounding boxes
	// in here, we create the return of this function
	// the bbNms.m function needs to be brought here too,
	// and this one will need some others.
/*
	bbs=cat(1,bbs{:});
	if(~isempty(pNms)), bbs=bbNms(bbs,pNms); end
*/

	return detections;
}
