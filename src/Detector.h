#include "Options.h"
#include "Clf.h"
#include "Info.h"
#include "BoundingBox.h"

using namespace cv;

class Detector
{
public:
	Options opts; //opts contains the Pyramid
	Clf clf;

	void importDetectorModel(String);
	inline void getChild(float*, uint32_t*, uint32_t*,	float*, uint32_t, uint32_t&, uint32_t&);
	BoundingBox* acfDetect(Mat);
};
