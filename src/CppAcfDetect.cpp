#include "CppAcfDetect.h"
#define DETECTOR_FILE "/home/gfuhr/Dropbox/Fuhr_Charles/QtAcfDetect/QtAcfDetect/detector.xml"

int main(int argc, char** argv) 
{
	Mat image;
    Detector d;

    QApplication a(argc, argv);
    /*MainWindow w;
    w.show();*/

    //call to read the detector file
    d.readDetectorModel(DETECTOR_FILE);

    BoundingBox* bbs = d.acfDetect(image);

    //destroyAllWindows();
    return a.exec();
}
