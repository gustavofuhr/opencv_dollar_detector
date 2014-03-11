class Clf
{
public:
    Mat fids; //unsigned integers
    Mat thrs; //
	unsigned int child[7][2048];
	float hs[7][2048];
	float weights[7][2048];
	unsigned int depth[7][2048];
    Mat errs;
	double losses[2048];
	int treeDepth;
};
