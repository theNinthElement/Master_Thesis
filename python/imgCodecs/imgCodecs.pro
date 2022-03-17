#-------------------------------------------------
#
# Project created by QtCreator 2019-10-28T17:33:13
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += console

TARGET = imgCodecs
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
        main.cpp \
        imgcodecs.cpp \
    randommask.cpp \
    supplementary_functions.cpp \
    masks.cpp \
    reggridmaskdialog.cpp \
    eedinpaintingdialog.cpp \
    inpainting.cpp \
    derivatives.cpp \
    fed.cpp \
    xorresdialog.cpp \
    diffresdialog.cpp \
    foeedinpaintingdialog.cpp \
    eedsmoothingdialog.cpp \
    smoothing.cpp \
    timeaxiscombinedialog.cpp \
    chain_coding.cpp \
    slicewiseeeddialog.cpp \
    eed_4_4d_dialog.cpp \
    directneighbordilationdialog.cpp \
    primizedialog.cpp \
    r_ilh_0_dialog.cpp

HEADERS += \
        imgcodecs.h \
    randommask.h \
    supplementary_functions.h \
    masks.h \
    reggridmaskdialog.h \
    eedinpaintingdialog.h \
    inpainting.h \
    derivatives.h \
    fed.h \
    fed_kappa.h \
    xorresdialog.h \
    diffresdialog.h \
    foeedinpaintingdialog.h \
    eedsmoothingdialog.h \
    smoothing.h \
    timeaxiscombinedialog.h \
    chain_coding.h \
    slicewiseeeddialog.h \
    eed_4_4d_dialog.h \
    directneighbordilationdialog.h \
    primizedialog.h \
    r_ilh_0_dialog.h

FORMS += \
        imgcodecs.ui \
    randommask.ui \
    reggridmaskdialog.ui \
    eedinpaintingdialog.ui \
    xorresdialog.ui \
    diffresdialog.ui \
    foeedinpaintingdialog.ui \
    eedsmoothingdialog.ui \
    timeaxiscombinedialog.ui \
    chain_coding.ui \
    slicewiseeeddialog.ui \
    eed_4_4d_dialog.ui \
    directneighbordilationdialog.ui \
    primizedialog.ui \
    r_ilh_0_dialog.ui

include(niftilib/nifti1.pri)
include(niftilib/nifti1_io.pri)
include(niftilib/znzlib.pri)
include(niftilib/ctpl_stl.pri)
include(cannylib/CannyEdgeDetector.pri)


INCLUDEPATH += /home/anirban/miniconda3/include #"/usr/local/include/opencv4" #/home/anirban/miniconda/include/opencv2
INCLUDEPATH += /home/anirban/miniconda3/include/python3.7m

LIBS += -L/home/anirban/miniconda3/lib \
-lopencv_ml \
-lopencv_objdetect \
-lopencv_stitching\
-lopencv_calib3d\
-lopencv_features2d\
-lopencv_highgui\
-lopencv_videoio\
-lopencv_imgcodecs\
-lopencv_video\
-lopencv_photo\
-lopencv_imgproc\
-lopencv_flann\
-lopencv_core\
#-lQt5Widgets\
#-lQt5Gui\
#-lQt5Core\
#-lGL\
#-lpthread

DEPENDPATH += -L/home/anirban/Qt/5.9.9/gcc_64/lib
LIBS += -L$$[QT_INSTALL_LIBS] -lQt5Gui -lQt5Core -lpthread -lGL

LIBS += -lpython3.8
#LIBS += =lopencv

