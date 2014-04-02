#include "Pyramid.h"
#include "PNms.h"
#include "PBoost.h"

using namespace cv;

class Options
{
public:
	Pyramid pPyramid;
	double modelDs[2];
	double modelDsPad[2];
	PNms pNms;
	int stride;
	int cascadeThreshold;
	double cascadeCalibration;
	int nWeak[4];
	PBoost pBoost;
	int seed;
	String name;
	String posGtDir;
	String posImgDir;
	String negImgDir;
	String posWinDir;
	String negWinDir;
	//function imreadf;
	Mat imreadp;
	double pLoadSquarify[2];
	int nPos;
	int nNeg;
	int nPerNeg;
	int nAccNeg;
	int pJitterFlip;
	int winsSave;
};
