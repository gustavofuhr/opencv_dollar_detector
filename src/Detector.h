#include "Options.h"
#include "Clf.h"
#include "Info.h"
#include "BoundingBox.h"

using namespace cv;

class Detector
{
	public:
    Options opts;
    Clf clf;
    Info info;

    void readDetectorModel(String);
    inline void getChild(float*, uint32_t*, uint32_t*,	float*, uint32_t, uint32_t&, uint32_t&);
    BoundingBox* acfDetect(Mat);
};
