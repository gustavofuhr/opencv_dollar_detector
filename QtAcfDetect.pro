#-------------------------------------------------
#
# Project created by QtCreator 2014-02-27T16:05:44
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = QtAcfDetect
TEMPLATE = app


SOURCES +=\
        src\mainwindow.cpp \
        src\CppAcfDetect.cpp \
        src\Detector.cpp

HEADERS  +=\
    src\mainwindow.h \
    src\Pyramid.h \
    src\PTree.h \
    src\PNms.h \
    src\PChns.h \
    src\PBoost.h \
    src\Options.h \
    src\Info.h \
    src\GradientMagnitude.h \
    src\GradientHistogram.h \
    src\Detector.h \
    src\CppAcfDetect.h \
    src\ColorChannels.h \
    src\Clf.h \
    src\BoundingBox.h

FORMS    += src\mainwindow.ui

INCLUDEPATH += D:\Programs\OpenCV
LIBS += -L/usr/local/lib
LIBS += -lopencv_core
LIBS += -lopencv_imgproc
LIBS += -lopencv_highgui
LIBS += -lopencv_ml
LIBS += -lopencv_video
LIBS += -lopencv_features2d
LIBS += -lopencv_calib3d
LIBS += -lopencv_objdetect
LIBS += -lopencv_contrib
LIBS += -lopencv_legacy
LIBS += -lopencv_flann
LIBS += -lopencv_nonfree
