#ifndef MASKS_H
#define MASKS_H

//~ Create an array of random different pixels.
double* randDiffPixels(double percent, double upperBound);

//~ Create an array of regular grid pixels.
double* regGridPixels(int samplingSize, double upperBound);

//~ Create an array of regular grid pixels for each slice seperately by given sampling size.
double* regGridPixelsSlicewise(int samplingSize, double imgWid, double imgHei, double imgDep, int& outArrSize);

//~ Create an array of regular grid pixels for volume by given sampling size.
double* regGridPixelsVolum(int samplingSize, double imgWid, double imgHei, double imgDep, int& outArrSize);

//~ Create a mask(2D=time-spatial) of regular grid pixels for volumes of 4D image by given sampling size.
double** regGridPixelsVolum_4D(int samplingSize, int imgWid, int imgHei, int imgDep, int imgTime);

//~ Moore Neighborhood Dilation Slicewise.
double* dilatedMask2D(int imgWid, int imgHei, int imgDep, int* alreadyTaken, int& dilatedMaskLen, bool excludIt);

//~ Moore Neighborhood 3D Dilation.
double* dilatedMask3D(int imgWid, int imgHei, int imgDep, int* alreadyTaken, int& dilatedMaskLen, bool excludIt);

//~ Partially Moore Neighborhood 3D Dilation by given mode size.
double* dilatedMask3D_partial(int imgWid, int imgHei, int imgDep, int* alreadyTaken, int& dilatedMaskLen, bool excludIt, int mode);

#endif // MASKS_H
