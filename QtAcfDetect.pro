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

HEADERS  += mainwindow.h \
    Pyramid.h \
    PTree.h \
    PNms.h \
    PChns.h \
    PBoost.h \
    Options.h \
    Info.h \
    GradientMagnitude.h \
    GradientHistogram.h \
    Detector.h \
    CppAcfDetect.h \
    ColorChannels.h \
    Clf.h \
    BoundingBox.h

FORMS    += mainwindow.ui

INCLUDEPATH += /usr/local/include/opencv
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
