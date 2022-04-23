#ifndef SUPPLEMENTARY_FUNCTIONS_H
#define SUPPLEMENTARY_FUNCTIONS_H

#include "niftilib/nifti1_io.h"
//#include "/usr/include/x86_64-linux-gnu/qt5/QtGui/QImage"
#include <QImage>
#include <vector>

using namespace std;

//Copy given image array to grayscale PNG.
void array_to_png(double* inputImg, QImage* imgObject , int imgWidth, int imgHeight);

//Copy given short int image array to grayscale PNG.
void si_array_to_png(short int* inputImg, QImage* imgObject , int imgWidth, int imgHeight);

// Copy NIftI format to 1D array.
double* nii_to_array(nifti_image *nim_input);

// Copy NIfTI format to 1D array of 2D images.
double** nii_to_arrayof_2d_arrays(nifti_image *nim_input);

// Copy 4D NIfTI format to 1D time array of 1D volume images.
double** nii4d_to_array(nifti_image *nim_input);

// Copy 4D NIfTI format to 1D time array of 1D volume images.  (for the second image e.g. for the mask)
double** nii4d_to_array_mask(nifti_image *nim_input);

// Copy given 1D(time axis)+1D double array entires to given 4D Nifti format object voxel values.
void array_to_nii4d(nifti_image *nim_input, double** inputImg);

// Copy given 1D(time axis)+1D double array entires to given 4D Nifti format object voxel values.    (for the second image e.g. for the mask)
void array_to_nii4d_mask(nifti_image *nim_input, double** inputImg);

// Copy 1D array of 2D images to NIfTI format.
void arrayof_2d_arrays_to_nii(nifti_image *nim_input, double** inputImgSlices);

// Copy NIfTI format to 3D arry.
double*** nii_to_3d_array(nifti_image *nim_input);

// Copy given 1D array entires to given Nifti format object voxel values.
void array_to_nii(nifti_image *nim_input, double* inputImg);

// Copy given 1D short int array entires to given Nifti format object voxel values.
void si_array_to_nii(nifti_image *nim_input, short int* inputImg);

// Copy NIftI format to 1D uchar array.
unsigned char* nii_to_ucharArray(nifti_image *nim_input);

//Scale any double array to uchar([0,255]) array
void scaleToUchar(double* imgArr, int len);

//Scale any short int array to uchar([0,255]) array
void si_scaleToUchar(short int* imgArr, int len);

// Create 1D Gaussian kernel.
void gaussian1D_kernel(double gKernelX[], int kernelLength, double sigma);

//Horizontal convolution.
double* convolutionX(double* inputAr, int inputArWidth, int inputArHeight, double gKernelX[], int kernelLength);

//Horizontal convolution for 3D image.
double* convolution3DX(double* inputAr, int inputArWidth, int inputArDepth, int inputArHeight, double gKernelX[], int kernelLength);

//Horizontal convolution for 3D image series of 4D data
double** convolution3DX_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArDepth, int inputArHeight, double gKernelX[], int kernelLength);

//Vertical convolution.
double* convolutionY(double* inputAr, int inputArWidth, int inputArHeight, double gKernelY[], int kernelLength);

//Vertical convolution for 3D image.
double* convolution3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength);

//Horizontal convolution for 3D image series of 4D data
double** convolution3DY_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength);

//Z direction convolution for 3D image.
double* convolution3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength);

//Z direction convolution for 3D image.
double** convolution3DZ_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength);

// T direction convolution of 4D image.
double** convolution3DT_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength);

//Mean Squared Error
double mse(double* imageArray1, double* imageArray2, int imageArrayLen);

//Mean Squared Error
double mse_4D(double** imageArray1, double** imageArray2, int imageArrayLen, int imgVolume);

//Mean Squared Error for Dilated Mask Region
double mse_mask_dilatied_region3D(double* imageArray1, double* imageArray2, int imgWidth, int imgHeight, int imgDepth, int* mskPxlLocs);

//Mean Squared Error for Dilated Mask Region
double mse_mask_dilatied_region4D(double** imageArray1, double** imageArray2, int imgWidth, int imgHeight, int imgDepth, int imgVolume, int* mskPxlLocs);


//Mean Squared Error of Sub Array (Region)
double sub_mse(double* imageArray1, double* imageArray2, double* maskArr, int imageArrayLen);

//Mean Squared Error of Sub Array (Region) for 4D (2D(time+data) array) data
double sub_mse_4D(double** imageArray1, double** imageArray2, double** maskArr, int timeLen, int imageArrayLen);

//Absolute Average Error
double aae(double* imageArray1, double* imageArray2, int imageArrayLen);

//Absolute Average Error
double aae_4D(double** imageArray1, double** imageArray2, int imageArrayLen, int imgVolume);

//Absolute Average Error for Dilated Mask Region
double aae_mask_dilatied_region3D(double* imageArray1, double* imageArray2, int imgWidth, int imgHeight, int imgDepth, int* mskPxlLocs);

//Absolute Average Error for Dilated Mask Region
double aae_mask_dilatied_region4D(double** imageArray1, double** imageArray2, int imgWidth, int imgHeight, int imgDepth, int imgVolume, int* mskPxlLocs);

//Absolute Average Error of Sub Array (Region)
double sub_aae(double* imageArray1, double* imageArray2, double* maskArr, int imageArrayLen);

//Absolute Average Error of Sub Array (Region) for 4D (2D(time+data) array) data
double sub_aae_4D(double** imageArray1, double** imageArray2, double** maskArr, int timeLen, int imageArrayLen);

//Discrete l2 norm
double l2Norm(double* imageArray1, double* imageArray2, int imageArrayLen);

//Discrete l2 norm for 4D (2D(time+data) array) data
double l2Norm_4D(double** imageArray1, double** imageArray2, int timeLen, int imageArrayLen);

//Charbonnier diffusivity function
double charbonnier_diff(double s_square, double contrastParam);

//Perona-Malik diffusivity function
double pm_diff(double s_square, double contrastParam);

//Aubert diffusivity function
double aubert_diff(double s_square, double contrastParam);

//Perona-Malik diffusivity function 2
double pm_diff2(double s_square, double contrastParam);

//Green diffusivity function
double green_diff(double s_square, double contrastParam);

//Geman et Reynolds
double gr_diff(double s_square, double contrastParam);

double li1(double s, double contPar);

//Forward sum on x axis.
double* sumForwX(double* inputAr, int inputArWidth, int inputArHeight);

//Backward sum on x axis.
double* sumBackX(double* inputAr, int inputArWidth, int inputArHeight);

//Forward sum on y axis.
double* sumForwY(double* inputAr, int inputArWidth, int inputArHeight);

//Backward derivative on y axis.
double* sumBackY(double* inputAr, int inputArWidth, int inputArHeight);

// 3D Forward sum on x axis.
double* sumForw3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth);

//3D Backward sum on x axis.
double* sumBack3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth);

//3D Forward sum on y axis.
double* sumForw3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth);

//3D Backward derivative on y axis.
double* sumBack3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth);

//3D Forward sum on z axis.
double* sumForw3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth);

//3D Backward derivative on z axis.
double* sumBack3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth);

// 4D Forward sum on x axis.
double* *sumForw4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Backward sum on x axis.
double* *sumBack4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Forward sum on y axis.
double* *sumForw4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Backward derivative on y axis.
double* *sumBack4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Forward sum on z axis.
double* *sumForw4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Backward derivative on z axis.
double* *sumBack4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Backward derivative on t axis.
double* *sumForw4DT(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);

//4D Backward derivative on t axis.
double* *sumBack4DT(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume);


//Number of primes less than given intiger
int countPrimes(int n);

// Function to generate N prime numbers using
// Sieve of Eratosthenes
void SieveOfEratosthenes(vector<short int> &primes, int maxVal);

// Utility function for Sieve of Sundaram
void sieveSundaram(vector<short int> &primes, int maxVal);

//Binary algorithm for the GCD
short int gcd(short int a, short int b);

// Function to perform Goldbach's conjecture
short int* findPrimes(vector<short int> primes, int n);

//Cross product of two vectors
void crossProd(double vec1[3], double vec2[3], double prod[3]);

// Function to convert 1D array to 2D slicewise array of 3D image
int** linear2twoDim(int* arr, int lenArr, int imgWid, int imgHei, int imgDep);

QQueue<int> array2queue(int* arr, int len);

//Nearest Neighbir Mask Initialization
double* nearNeighInit(double* arr, int width, int height, int depth, int* maskInds);

// Method to compare which one is the more close. It assumes that val2 is greater than val1 and target lies between these two.
int getClosest(int val1, int val2, int target);

// Returns element closest to target in array
int findClosest(vector<short int> arr, int n, int target);

//Convert all voxel values of an image to the closest prime numbers
void primize(double* imgArr, int imgWid, int imgHei, int imgDep, int maxVoxelVal);

//Convert all voxel values of an image to the number of primes less than or equal to voxel value
void primeCountImage(double* imgArr, int imgWid, int imgHei, int imgDep);

// Calculate Image "qntl" Quantile given in percentage.
double quantlCriterPM(double qntl, double* imgArr, int imgWidth, int imgHeight, int imgDepth, double sigma, int kernelSize);

// Calculate Image "qntl" Quantile given in percentage.
double quantlCriterPM4Inpainting(double qntl, double* imgArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, double sigma, int kernelSize);

// Calculate Image "qntl" Quantile given in percentage.
double quantlCriterPM4Inpainting4D(double qntl, double* imgArr, double* binMask, int imgWidth, int imgHeight, int imgDepth, double sigma, int kernelSize);

// Calculate Image "qntl" Quantile given in percentage.
double* quantlCriterPM4Inpainting4DNoBin(double qntl, double** imgArr, int** randPxls, int imgWidth, int imgHeight, int imgDepth, int imgVolume, double sigma, int kernelSize);

// Calculate Image "qntl" Quantile given in percentage.
double histoCriterPM4Inpainting(double qntl, double* imgArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth);

//Entropy Computation.
double entropy(double* arr, int len);

//Sum of squared differences (SSD) array of 4D(3D + q space) data in the x direction.
double* ssd_4DX_central(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth);

//Sum of squared differences (SSD) array of 4D(3D + q space) data in the y direction.
double* ssd_4DY_central(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth);

//Sum of squared differences (SSD) array of 4D(3D + q space) data in the z direction.
double* ssd_4DZ_central(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth);

//For non-dropouts: Sum of squared differences (SSD) array of 4D(3D + q space) data in the x direction.
double* ssd_4DX_central_nonDropouts(double** inputAr, double** maskArr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth);

//For non-dropouts: Sum of squared differences (SSD) array of 4D(3D + q space) data in the y direction.
double* ssd_4DY_central_nonDropouts(double** inputAr, double** maskArr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth);

//For non-dropouts: Sum of squared differences (SSD) array of 4D(3D + q space) data in the Z direction.
double* ssd_4DZ_central_nonDropouts(double** inputAr, double** maskArr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth);

//~ Create an array of random n different voxel locations (from the same volume) from the mask image for 4D data.
//double* randVoxLocs(int locLen, double** maskImg, int volNum, int upperBound, int lowerBound);
double* randVoxLocs(int locLen, double** maskImg, int volNum, int upperBound);

// Return corrupted volume numbers a vector of intigers
vector<int> corruptedVolumes(double** binMaskArr, int imgWidth, int imgHeight, int imgDepth, int numVols);

// Given a volume return corrupted slice numbers as a vector of integers
vector<int> corruptedSlices(double* binMaskArr, int imgWidth, int imgHeight, int imgDepth);

// Make a regular grid mask from a slice mask. (Currently, there are 4 modes, TR, TL, BR, BL)
void makeRegGrid_on_slice(double* binMaskrr, int imgWidth, int imgHeight, int currSlice, int gridLen, int mode);

#endif // SUPPLEMENTARY_FUNCTIONS_H
