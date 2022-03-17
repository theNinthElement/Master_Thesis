#include <cassert>
#include <cstdlib>
#include <cmath>

#include <qdebug.h>

//~ Create an array of random different pixels.
double* randDiffPixels(double percent, double upperBound)
{
    int length;
    double n,m;

    length = (int)(upperBound*percent)/100;

    // Decleartion and memory allocation
    double* randPixelsArr = new double[sizeof(double) * length];
//    srand(time(NULL));    //For generating different random scattered images
    srand(1);

    m = 0;
    for(n=0; n<upperBound && m<length; ++n) {
        double rn = upperBound - n;
        double rm = length - m;
        if(rand() % (int)rn < rm)
            randPixelsArr[(int)m++] = n;
    }
    assert(m == length);
    qDebug() << "Random generated once";
    return randPixelsArr;
}

//~ Create an array of regular grid pixels by given ratio size.
double* regGridPixels(int samplingSize, double upperBound)
{
        int length;
        length = (int)(upperBound/samplingSize);

//        qDebug() << length+1;
        // Decleartion and memory allocation
        double* regGridPixelsArr = new double[sizeof(double) * (length+1)];

        int n = 0;
        for(int i=0; i<upperBound; i++) {
            if(i%samplingSize == 0) {
                regGridPixelsArr[n] = (double)i;
//                printf("Grid array: %lf\n", regGridPixelsArr[n]);
                n++;
            }
        }
//        qDebug() << n;
        return regGridPixelsArr;
}

//~ Create an array of regular grid pixels for each slice seperately by given sampling size.
double* regGridPixelsSlicewise(int samplingSize, double imgWid, double imgHei, double imgDep, int& outArrSize)
{
        int length;
        length = (int)(ceil(imgWid/samplingSize)*ceil(imgHei/samplingSize)*imgDep);

        qDebug() << "Expected Mask Size: " << length;
        // Decleartion and memory allocation
        double* regGridPixelsArr = new double[sizeof(double) * (length)];

        int n = 0;
        for(int d=0; d<imgDep; d++) {
            for(int k=0; k<imgHei; k=k+samplingSize) {
                for(int i=0; i<imgWid; i++) {
                    if(i%samplingSize == 0) {
                        regGridPixelsArr[n] = (double)(i+imgWid*k+imgWid*imgHei*d);
        //                printf("Grid array: %lf\n", regGridPixelsArr[n]);
                        n++;
                    }
                }
            }
        }
        qDebug() << "Reg. Mask Size: " << n;
        outArrSize = n;

        return regGridPixelsArr;
}

//~ Create an array of regular grid pixels for volume by given sampling size.
double* regGridPixelsVolum(int samplingSize, double imgWid, double imgHei, double imgDep, int& outArrSize)
{
        int length;
        length = (int)(ceil(imgWid/samplingSize)*ceil(imgHei/samplingSize)*ceil(imgDep/samplingSize));

        qDebug() << "Expected Mask Size: " << length;
        // Decleartion and memory allocation
        double* regGridPixelsArr = new double[sizeof(double) * (length)];

        int n = 0;
        for(int d=0; d<imgDep; d=d+samplingSize) {
            for(int k=0; k<imgHei; k=k+samplingSize) {
                for(int i=0; i<imgWid; i++) {
                    if(i%samplingSize == 0) {
                        regGridPixelsArr[n] = (double)(i+imgWid*k+imgWid*imgHei*d);
        //                printf("Grid array: %lf\n", regGridPixelsArr[n]);
                        n++;
                    }
                }
            }
        }
        qDebug() << "Reg. Mask Size: " << n;
        outArrSize = n;

        return regGridPixelsArr;
}

//~ Create a mask(2D=time-spatial) of regular grid pixels for volumes of 4D image by given sampling size.
double** regGridPixelsVolum_4D(int samplingSize, int imgWid, int imgHei, int imgDep, int imgTime)
{
    // Decleartion and memory allocation
    double** regGridBinMask = new double*[sizeof(double)*imgTime];
    for(int i = 0; i < imgTime; ++i) {
        regGridBinMask[i] = new double[sizeof(double)*imgWid*imgHei*imgDep];
    }
    for(int i = 0; i < imgTime; ++i) {
        for(int j=0; j<imgWid*imgHei*imgDep; j++) {
            regGridBinMask[i][j] = 0;
        }
    }

//        for(int i = 0; i < imgTime; ++i) {
//            for(int j=0; j<imgWid*imgHei*imgDep; j++) {
//                qDebug() << regGridBinMask[i][j] << i << j;
//            }
//        }

    for(int t=0; t<imgTime; t++) {
        for(int d=0; d<imgDep; d=d+samplingSize) {
            for(int k=0; k<imgHei; k=k+samplingSize) {
                for(int i=0; i<imgWid; i=i+samplingSize) {
                    regGridBinMask[t][i+imgWid*k+imgWid*imgHei*d] = 1;
//                            qDebug() << regGridBinMask[t][i+imgWid*k+imgWid*imgHei*d] << t << i+imgWid*k+imgWid*imgHei*d;
                }
            }
        }
    }
    return regGridBinMask;
}

//~ Moore Neighborhood Dilation Slicewise.
double* dilatedMask2D(int imgWid, int imgHei, int imgDep, int* alreadyTaken, int& dilatedMaskLen, bool excludIt)
{
    // Decleartion and memory allocation
    int* binaryMaskArr = new int[sizeof(int) * imgWid*imgHei*imgDep];

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        binaryMaskArr[i] = 0;
    }
    int arrTraceIndex=0, j=0;

    for(int k=0; k<imgDep; k++) {
        for(int i=0; i<imgWid*imgHei; i++) {
            int depth_i = i + imgWid*imgHei*k;
            if(depth_i == alreadyTaken[arrTraceIndex]) {
                //Make immediate neighbors also selected
                if(i == 0) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i+1] = 1;
                    binaryMaskArr[depth_i+imgWid] = 1;
                    binaryMaskArr[depth_i+imgWid+1] = 1;
                } else if(i == imgWid-1) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i+imgWid] = 1;
                    binaryMaskArr[depth_i+imgWid-1] = 1;
                } else if(i == imgWid*(imgHei-1)) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i+1] = 1;
                    binaryMaskArr[depth_i-imgWid+1] = 1;
                } else if(i == imgWid*imgHei-1) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i-imgWid] = 1;
                    binaryMaskArr[depth_i-imgWid-1] = 1;
                } else if(i > 0 && i < imgWid-1) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i+1] = 1;
                    binaryMaskArr[depth_i+imgWid] = 1;
                    binaryMaskArr[depth_i+imgWid-1] = 1;
                    binaryMaskArr[depth_i+imgWid+1] = 1;
                } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i+1] = 1;
                    binaryMaskArr[depth_i-imgWid] = 1;
                    binaryMaskArr[depth_i-imgWid-1] = 1;
                    binaryMaskArr[depth_i-imgWid+1] = 1;
                } else if (i % imgWid == 0) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-imgWid] = 1;
                    binaryMaskArr[depth_i+1] = 1;
                    binaryMaskArr[depth_i+imgWid] = 1;
                    binaryMaskArr[depth_i-imgWid+1] = 1;
                    binaryMaskArr[depth_i+imgWid+1] = 1;
                } else if ((i+1) % imgWid == 0) {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-imgWid] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i+imgWid] = 1;
                    binaryMaskArr[depth_i-imgWid-1] = 1;
                    binaryMaskArr[depth_i+imgWid-1] = 1;
                }  else {
                    binaryMaskArr[depth_i] = 1;
                    binaryMaskArr[depth_i-imgWid] = 1;
                    binaryMaskArr[depth_i-imgWid-1] = 1;
                    binaryMaskArr[depth_i-imgWid+1] = 1;
                    binaryMaskArr[depth_i-1] = 1;
                    binaryMaskArr[depth_i+1] = 1;
                    binaryMaskArr[depth_i+imgWid] = 1;
                    binaryMaskArr[depth_i+imgWid-1] = 1;
                    binaryMaskArr[depth_i+imgWid+1] = 1;
                }
                arrTraceIndex++;
                continue;
            }
        }
    }

    if(excludIt) {
        int randArrTraceIndex = 0;
        for(int i=0; i<imgWid*imgHei*imgDep; i++) {
            if(i == alreadyTaken[randArrTraceIndex]) {
                randArrTraceIndex++;
                binaryMaskArr[i] = 0;
                continue;
            }
        }
        qDebug() << "Prev. Mask Size: " << randArrTraceIndex;
    }

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(binaryMaskArr[i] == 1) {
            dilatedMaskLen++;
        }
//        qDebug() << "Binary Mask: " << binaryMaskArr[i];
    }

    // Decleartion and memory allocation
    double* randPixelsArr = new double[sizeof(double) * dilatedMaskLen];
    int tempIndx=0;
    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(binaryMaskArr[i] == 1) {
            randPixelsArr[tempIndx] = i;
            tempIndx++;
        }
    }
//    qDebug() << "Dilated Mask Len: " << dilatedMaskLen;
//    for(int i=0; i<dilatedMaskLen; i++) {
//        qDebug() << "After All: " << i << randPixelsArr[i];
//    }

    return randPixelsArr;
}

//~ Moore Neighborhood 3D Dilation.
double* dilatedMask3D(int imgWid, int imgHei, int imgDep, int* alreadyTaken, int& dilatedMaskLen, bool excludIt)
{
    // Decleartion and memory allocation
    int* binaryMaskArr = new int[sizeof(int) * imgWid*imgHei*imgDep];

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        binaryMaskArr[i] = 0;
    }
    int arrTraceIndex=0, j=0;

    for(int k=0; k<imgDep; k++) {
        for(int i=0; i<imgWid*imgHei; i++) {
            int depth_i = i + imgWid*imgHei*k;
            if(depth_i == alreadyTaken[arrTraceIndex]) {
                if(k == 0) {
                    //Make immediate neighbors also selected
                    if(i == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    } else if(i == imgWid-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                    } else if(i == imgWid*(imgHei-1)) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                    } else if(i == imgWid*imgHei-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                    } else if(i > 0 && i < imgWid-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                    } else if (i % imgWid == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    } else if ((i+1) % imgWid == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                    }  else {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    }
                } else if(k == imgDep-1) {
                    //Make immediate neighbors also selected
                    if(i == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                    } else if(i == imgWid-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                    } else if(i == imgWid*(imgHei-1)) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                    } else if(i == imgWid*imgHei-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                    } else if(i > 0 && i < imgWid-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                    } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                    } else if (i % imgWid == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                    } else if ((i+1) % imgWid == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                    }  else {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                    }
                } else {
                    //Make immediate neighbors also selected
                    if(i == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    } else if(i == imgWid-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                    } else if(i == imgWid*(imgHei-1)) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                    } else if(i == imgWid*imgHei-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                    } else if(i > 0 && i < imgWid-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                    } else if (i % imgWid == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    } else if ((i+1) % imgWid == 0) {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                    }  else {
                        binaryMaskArr[depth_i] = 1;
                        binaryMaskArr[depth_i-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid+1] = 1;
                        binaryMaskArr[depth_i-1] = 1;
                        binaryMaskArr[depth_i+1] = 1;
                        binaryMaskArr[depth_i+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i-imgWid*imgHei+imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-imgWid+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid-1] = 1;
                        binaryMaskArr[depth_i+imgWid*imgHei+imgWid+1] = 1;
                    }
                }
                arrTraceIndex++;
                continue;
            }
        }
    }

    if(excludIt) {
        int randArrTraceIndex = 0;
        for(int i=0; i<imgWid*imgHei*imgDep; i++) {
            if(i == alreadyTaken[randArrTraceIndex]) {
                randArrTraceIndex++;
                binaryMaskArr[i] = 0;
                continue;
            }
        }
        qDebug() << "Prev. Mask Size: " << randArrTraceIndex;
    }

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(binaryMaskArr[i] == 1) {
            dilatedMaskLen++;
        }
//        qDebug() << "Binary Mask: " << binaryMaskArr[i];
    }

    // Decleartion and memory allocation
    double* randPixelsArr = new double[sizeof(double) * dilatedMaskLen];
    int tempIndx=0;
    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(binaryMaskArr[i] == 1) {
            randPixelsArr[tempIndx] = i;
            tempIndx++;
        }
    }
//    qDebug() << "Dilated Mask Len: " << dilatedMaskLen;
//    for(int i=0; i<dilatedMaskLen; i++) {
//        qDebug() << "After All: " << i << randPixelsArr[i];
//    }

    return randPixelsArr;
}

//~ Partially Moore Neighborhood 3D Dilation by given mode size.
double* dilatedMask3D_partial(int imgWid, int imgHei, int imgDep, int* alreadyTaken, int& dilatedMaskLen, bool excludIt, int mode)
{
    // Decleartion and memory allocation
    int* binaryMaskArr = new int[sizeof(int) * imgWid*imgHei*imgDep];

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        binaryMaskArr[i] = 0;
    }
    int arrTraceIndex=0, j=0;

    if(mode == 1) {
        for(int k=0; k<imgDep; k++) {
            for(int i=0; i<imgWid*imgHei; i++) {
                int depth_i = i + imgWid*imgHei*k;
                if(depth_i == alreadyTaken[arrTraceIndex]) {
                    if(k == 0) {
                        //Make immediate neighbors also selected
                        if(i == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i == imgWid-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i == imgWid*(imgHei-1)) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i == imgWid*imgHei-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i>0 && i<imgWid-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if (i % imgWid == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if ((i+1) % imgWid == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        }  else {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        }
                    } else if(k == imgDep-1) {
                        //Make immediate neighbors also selected
                        if(i == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if(i == imgWid-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if(i == imgWid*(imgHei-1)) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if(i == imgWid*imgHei-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if(i > 0 && i < imgWid-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if (i % imgWid == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        } else if ((i+1) % imgWid == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                        }  else {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                       }
                    } else {
                        //Make immediate neighbors also selected
                        if(i == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i == imgWid-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i == imgWid*(imgHei-1)) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i == imgWid*imgHei-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i > 0 && i < imgWid-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if(i > imgWid*(imgHei-1) && i < imgWid*imgHei-1) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if (i % imgWid == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        } else if ((i+1) % imgWid == 0) {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        }  else {
                            binaryMaskArr[depth_i] = 1;
                            binaryMaskArr[depth_i-imgWid] = 1;
                            binaryMaskArr[depth_i-1] = 1;
                            binaryMaskArr[depth_i+1] = 1;
                            binaryMaskArr[depth_i+imgWid] = 1;
                            binaryMaskArr[depth_i-imgWid*imgHei] = 1;
                            binaryMaskArr[depth_i+imgWid*imgHei] = 1;
                        }
                    }
                    arrTraceIndex++;
                    continue;
                }
            }
        }
    }

    if(excludIt) {
        int randArrTraceIndex = 0;
        for(int i=0; i<imgWid*imgHei*imgDep; i++) {
            if(i == alreadyTaken[randArrTraceIndex]) {
                randArrTraceIndex++;
                binaryMaskArr[i] = 0;
                continue;
            }
        }
        qDebug() << "Prev. Mask Size: " << randArrTraceIndex;
    }

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(binaryMaskArr[i] == 1) {
            dilatedMaskLen++;
        }
//        qDebug() << "Binary Mask: " << binaryMaskArr[i];
    }

    // Decleartion and memory allocation
    double* randPixelsArr = new double[sizeof(double) * dilatedMaskLen];
    int tempIndx=0;
    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(binaryMaskArr[i] == 1) {
            randPixelsArr[tempIndx] = i;
            tempIndx++;
        }
    }
//    qDebug() << "Dilated Mask Len: " << dilatedMaskLen;
//    for(int i=0; i<dilatedMaskLen; i++) {
//        qDebug() << "After All: " << i << randPixelsArr[i];
//    }

    return randPixelsArr;
}
