#include "PChns.h"

class Pyramid
{
	public:
		PChns pChns;
		int nPerOct;
		int nOctUp;
		int nApprox;
		double lambdas[3];
		int pad[2];
		int minDs[2];
		int smooth;
		int concat;
		int complete;
		int nScales;

		Pyramid channelFeaturePyramid(Mat);
};