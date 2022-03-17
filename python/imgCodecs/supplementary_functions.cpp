#include <cmath>
#include <QImage>
#include <vector>
#include <QQueue>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <cassert>
#include "nifti1_io.h"
#include "derivatives.h"
#include "masks.h"

#include <QDebug>
using namespace std;

//Copy given double image array to grayscale PNG.
void array_to_png(double* inputImg, QImage* imgObject , int imgWidth, int imgHeight) {

    for(int i=0; i<imgHeight; i++) {
        uchar* imageCharArr = imgObject->scanLine(i);
        for(int j=0; j<imgWidth; j++) {
            double tempPxlVal;
//            qDebug() << "Before Pixel value at " << "(" << i << "," << j << ")" << inputImg[j+i*imgWidth] << "\n";
//            qDebug("Pixel value at (%d,%d) before: %u",i,j,imageCharArr[j]);
            tempPxlVal = inputImg[j+i*imgWidth];
            if(tempPxlVal < 0) {
                tempPxlVal = 0;
            } else if(tempPxlVal > 255) {
                tempPxlVal = 255;
            } else {
                tempPxlVal = round(tempPxlVal);
            }
            imageCharArr[j] = (unsigned char)tempPxlVal;
//            imageCharArr[j] = tempPxlVal;
//            qDebug("Pixel value at (%d,%d): %u",i,j,imageCharArr[j]);
//            qDebug() << "Pixel value at " << "(" << i << "," << j << ")" << imageCharArr[j] << "\n";
        }
    }
}

//Copy given short int image array to grayscale PNG.
void si_array_to_png(short int* inputImg, QImage* imgObject , int imgWidth, int imgHeight) {

    for(int i=0; i<imgHeight; i++) {
        uchar* imageCharArr = imgObject->scanLine(i);
        for(int j=0; j<imgWidth; j++) {
            short int tempPxlVal;
//            qDebug() << "Before Pixel value at " << "(" << i << "," << j << ")" << inputImg[j+i*imgWidth] << "\n";
//            qDebug("Pixel value at (%d,%d) before: %u",i,j,imageCharArr[j]);
            tempPxlVal = inputImg[j+i*imgWidth];
            if(tempPxlVal < 0) {
                tempPxlVal = 0;
            } else if(tempPxlVal > 255) {
                tempPxlVal = 255;
            } else {
                tempPxlVal = round(tempPxlVal);
            }
            imageCharArr[j] = (unsigned char)tempPxlVal;
//            imageCharArr[j] = tempPxlVal;
//            qDebug("Pixel value at (%d,%d): %u",i,j,imageCharArr[j]);
//            qDebug() << "Pixel value at " << "(" << i << "," << j << ")" << imageCharArr[j] << "\n";
        }
    }
}

// Copy NIfTI format to 1D double array.
double* nii_to_array(nifti_image *nim_input)
{
    double* imageArr;
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3];
    int dataType = nim_input->datatype;

    // Memory allocation for imageArr.
    imageArr = (double*) calloc(imgWid*imgHei*imgDep, sizeof(double));
    if(imageArr == NULL) /* Memory allocation fails */ {
        printf("Couldn't allocate memory for 'nii_to_array' array\n");
    }
//    else /* Memory allocation successful */ {
//        printf("Memory allocation successful for 'nii_to_array' array\n");
//    }

    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

//    char xyztUnits = nim_input->xyz_units;
//    short form = nim_input->qform_code;
//    qDebug() << "xyzt units" << xyztUnits;
//    qDebug() << "Form Code: " << form;

    // Get access to data of nim_input
    if(dataType == 0) {
        qDebug() << "Unknown datatype!";
        return 0;
    } else if(dataType == 1) {
        qDebug() << "Datatype: Bool (1 bit)";
//        bool *nim_input_data = (bool *) nim_input->data;
    } else if(dataType == 2) {
        qDebug() << "Datatype: Unsigned Char (8 bits)";
//        unsigned char *nim_input_data = (unsigned char *) nim_input->data;
    } else if(dataType == 4) {
        qDebug() << "Datatype: Signed Short (16 bits)";
//        signed short *nim_input_data = (signed short *) nim_input->data;
    } else if(dataType == 8) {
        qDebug() << "Datatype: Signed Int (32 bits)";
//        signed int *nim_input_data = (signed int *) nim_input->data;
    } else if(dataType == 16) {
        qDebug() << "Datatype: Float (32 bits)";
//        float *nim_input_data = (float *) nim_input->data;
    } else if(dataType == 64) {
        qDebug() << "Datatype: Double (64 bits)";
//        double *nim_input_data = (double *) nim_input->data;
    } else if(dataType == 256) {
        qDebug() << "Datatype: Signed Char (8 bits)";
//        signed char *nim_input_data = (signed char *) nim_input->data;
    } else if(dataType == 512) {
        qDebug() << "Datatype: Unsigned Short (16 bits)";
//        unsigned short *nim_input_data = (unsigned short *) nim_input->data;
    } else if(dataType == 768) {
        qDebug() << "Datatype: Unsigned Int (32 bits)";
//        unsigned int *nim_input_data = (unsigned int *) nim_input->data;
    }
    //Depending on the nii voxel datatype shoul be chanched.
//    unsigned char *nim_input_data = (unsigned char *) nim_input->data;
//    float  *nim_input_data = (float *) nim_input->data;
//    short int *nim_input_data = (short int *) nim_input->data;
//    signed short *nim_input_data = (signed short *) nim_input->data;
    unsigned short *nim_input_data = (unsigned short *) nim_input->data;

    qDebug() << "So far successfully";

    int counter = 0;
    for(int timestep=0; timestep<nrep; ++timestep){
        for(int islice=0; islice<sizeSlice; ++islice){
            for(int iy=0; iy<sizeRead; ++iy){
                for(int ix=0; ix<sizePhase; ++ix){
//                    float tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix + iy);
//                    short int tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
//                    short int tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    signed short tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];       //Depending on the nii voxel datatype shoul be chanched.
                    unsigned short tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    float tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    unsigned char tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    imageArr[counter] = (double)tempValu;
//                    printf("Voxel values: %f\n", imageArr[counter]);
                    imageArr[counter] = (double)tempValu;
//                    printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
    return imageArr;
}

// Copy NIfTI format to 1D array of 2D images.
double** nii_to_arrayof_2d_arrays(nifti_image *nim_input)
{
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3];

    // Memory allocation for imageArr.
    double** imageArr = new double*[sizeof(double)*imgDep];
    for(int i = 0; i < imgDep; ++i) {
        imageArr[i] = new double[sizeof(double)*imgWid*imgHei];
    }

    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    // Get access to data of nim_input
//    float  *nim_input_data = (float *) nim_input->data;
    short int *nim_input_data = (short int *) nim_input->data;

    for(int timestep=0; timestep<nrep; ++timestep) {
        for(int islice=0; islice<sizeSlice; ++islice) {
            int counter = 0;
            for(int iy=0; iy<sizeRead; ++iy) {
                for(int ix=0; ix<sizePhase; ++ix) {
 //                    float tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix + iy);
 //                    short int tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
                    short int tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
 //                    imageArr[counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    imageArr[islice][counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
    return imageArr;
}

// Copy 4D NIfTI format to 1D time array of 1D volume images.
double** nii4d_to_array(nifti_image *nim_input)
{
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3], imgTime = nim_input->dim[4];

    // Memory allocation for imageArr.
    double** imageArr = new double*[sizeof(double)*imgTime];
    for(int i = 0; i < imgTime; ++i) {
        imageArr[i] = new double[sizeof(double)*imgWid*imgHei*imgDep];
    }

    qDebug() << " nii4d_ declaration done";
    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    // Get access to data of nim_input
    float  *nim_input_data = (float *) nim_input->data;
//    double *nim_input_data = (double *) nim_input->data;
//    short int *nim_input_data = (short int *) nim_input->data;

    qDebug() << " nii4d_ convertion for nrep : " << nrep;
    qDebug() << " nii4d_ convertion for sizeSlice : " << sizeSlice;
    qDebug() << " nii4d_ convertion for sizeRead : " << sizeRead;
    qDebug() << " nii4d_ convertion for sizePhase : " << sizePhase;
    for(int timestep=0; timestep<nrep; ++timestep) {
        int counter = 0;
        for(int islice=0; islice<sizeSlice; ++islice) {
            for(int iy=0; iy<sizeRead; ++iy) {
                for(int ix=0; ix<sizePhase; ++ix) {
 //                    float tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix + iy);
 //                    short int tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
                    float tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    double tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];                        //Depending on the nii voxel datatype shoul be chanched.
//                    short int tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
 //                    imageArr[counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    imageArr[timestep][counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
        qDebug() << " nii4d_ convertion for volume : " << timestep;
    }
    return imageArr;
}

// Copy 4D NIfTI format to 1D time array of 1D volume images.  (for the second image e.g. for the mask)
double** nii4d_to_array_mask(nifti_image *nim_input)
{
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3], imgTime = nim_input->dim[4];

    // Memory allocation for imageArr.
    double** imageArr = new double*[sizeof(double)*imgTime];
    for(int i = 0; i < imgTime; ++i) {
        imageArr[i] = new double[sizeof(double)*imgWid*imgHei*imgDep];
    }

    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    // Get access to data of nim_input
    float  *nim_input_data = (float *) nim_input->data;
//    double *nim_input_data = (double *) nim_input->data;
//    short int *nim_input_data = (short int *) nim_input->data;
//    signed short *nim_input_data = (signed short *) nim_input->data;

    for(int timestep=0; timestep<nrep; ++timestep) {
        int counter = 0;
        for(int islice=0; islice<sizeSlice; ++islice) {
            for(int iy=0; iy<sizeRead; ++iy) {
                for(int ix=0; ix<sizePhase; ++ix) {
 //                    float tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix + iy);
 //                    short int tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
                    float tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    double tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];                        //Depending on the nii voxel datatype shoul be chanched.
//                    short int tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
//                    signed short tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
 //                    imageArr[counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    imageArr[timestep][counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
    return imageArr;
}

// Copy given 1D(time axis)+1D double array entires to given 4D Nifti format object voxel values.
void array_to_nii4d(nifti_image *nim_input, double** inputImg)
{
    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep = nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;
    int dataType = nim_input->datatype;

    // Get access to data of nim_input
    if(dataType == 0) {
        qDebug() << "Unknown datatype!";
        return;
    } else if(dataType == 1) {
        qDebug() << "Datatype: Bool (1 bit)";
//        bool *nim_input_data = (bool *) nim_input->data;
    } else if(dataType == 2) {
        qDebug() << "Datatype: Unsigned Char (8 bits)";
//        unsigned char *nim_input_data = (unsigned char *) nim_input->data;
    } else if(dataType == 4) {
        qDebug() << "Datatype: Signed Short (16 bits)";
//        signed short *nim_input_data = (signed short *) nim_input->data;
    } else if(dataType == 8) {
        qDebug() << "Datatype: Signed Int (32 bits)";
//        signed int *nim_input_data = (signed int *) nim_input->data;
    } else if(dataType == 16) {
        qDebug() << "Datatype: Float (32 bits)";
//        float *nim_input_data = (float *) nim_input->data;
    } else if(dataType == 64) {
        qDebug() << "Datatype: Double (64 bits)";
//        double *nim_input_data = (double *) nim_input->data;
    } else if(dataType == 256) {
        qDebug() << "Datatype: Signed Char (8 bits)";
//        signed char *nim_input_data = (signed char *) nim_input->data;
    } else if(dataType == 512) {
        qDebug() << "Datatype: Unsigned Short (16 bits)";
//        unsigned short *nim_input_data = (unsigned short *) nim_input->data;
    } else if(dataType == 768) {
        qDebug() << "Datatype: Unsigned Int (32 bits)";
//        unsigned int *nim_input_data = (unsigned int *) nim_input->data;
    }
    //Depending on the nii voxel datatype shoul be chanched.
//    unsigned char *nim_input_data = (unsigned char *) nim_input->data;
//    float  *nim_input_data = (float *) nim_input->data;
//    double *nim_input_data = (double *) nim_input->data;
//    short int  *nim_input_data = (short int *) nim_input->data;
    signed short *nim_input_data = (signed short *) nim_input->data;

    for(int timestep=0; timestep<nrep; ++timestep){
        int counter = 0;
        for(int islice=0; islice<sizeSlice; ++islice){
            for(int iy=0; iy<sizeRead; ++iy){
                for(int ix=0; ix<sizePhase; ++ix){
//                    unsigned short tempValu = (unsigned short)inputImg[counter];
//                    float tempValu = (float)inputImg[timestep][counter];                                    //Depending on the nii voxel datatype shoul be chanched.
//                    double tempValu = (double)inputImg[timestep][counter];
                    signed short tempValu = (signed short)inputImg[timestep][counter];
                    *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix) = tempValu;
//                    nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix] = tempValu;
    //              printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
}

// Copy 1D array of 2D images to NIfTI format.
void arrayof_2d_arrays_to_nii(nifti_image *nim_input, double** inputImgSlices)
{
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3];

    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    // Get access to data of nim_input
//    float  *nim_input_data = (float *) nim_input->data;
    short int *nim_input_data = (short int *) nim_input->data;

    for(int timestep=0; timestep<nrep; ++timestep) {
        for(int islice=0; islice<sizeSlice; ++islice) {
            int counter = 0;
            for(int iy=0; iy<sizeRead; ++iy) {
                for(int ix=0; ix<sizePhase; ++ix) {
                    short int tempValu = (short int)inputImgSlices[islice][counter];
                    *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix) = tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
}

// Copy NIfTI format to 3D arry.
double*** nii_to_3d_array(nifti_image *nim_input)
{
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3];

    // Memory allocation for imageArr.
    double*** imageArr = new double**[sizeof(double)*imgDep];
    for(int i = 0; i < imgDep; ++i) {
        imageArr[i] = new double*[sizeof(double)*imgHei];
        for(int j = 0; j < imgHei; ++j) {
            imageArr[i][j] = new double[sizeof(double)*imgWid];
        }

    }

    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    // Get access to data of nim_input
//    float  *nim_input_data = (float *) nim_input->data;
    short int *nim_input_data = (short int *) nim_input->data;

    for(int timestep=0; timestep<nrep; ++timestep) {
        for(int islice=0; islice<sizeSlice; ++islice) {
            for(int iy=0; iy<sizeRead; ++iy) {
                for(int ix=0; ix<sizePhase; ++ix) {
//                    int counter = 0;
 //                    float tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix + iy);
 //                    short int tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
                    short int tempValu = nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix];
 //                    imageArr[counter] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
                    imageArr[islice][iy][ix] = (double)tempValu;
 //                    printf("Voxel values: %f\n", imageArr[counter]);
//                    counter++;
                }
            }
        }
    }
    return imageArr;
}

// Copy given 1D double array entires to given Nifti format object voxel values.
void array_to_nii(nifti_image *nim_input, double* inputImg)
{
    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep = nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;
    int dataType = nim_input->datatype;

    // Get access to data of nim_input
    if(dataType == 0) {
        qDebug() << "Unknown datatype!";
        return;
    } else if(dataType == 1) {
        qDebug() << "Datatype: Bool (1 bit)";
//        bool *nim_input_data = (bool *) nim_input->data;
    } else if(dataType == 2) {
        qDebug() << "Datatype: Unsigned Char (8 bits)";
//        unsigned char *nim_input_data = (unsigned char *) nim_input->data;
    } else if(dataType == 4) {
        qDebug() << "Datatype: Signed Short (16 bits)";
//        signed short *nim_input_data = (signed short *) nim_input->data;
    } else if(dataType == 8) {
        qDebug() << "Datatype: Signed Int (32 bits)";
//        signed int *nim_input_data = (signed int *) nim_input->data;
    } else if(dataType == 16) {
        qDebug() << "Datatype: Float (32 bits)";
//        float *nim_input_data = (float *) nim_input->data;
    } else if(dataType == 64) {
        qDebug() << "Datatype: Double (64 bits)";
//        double *nim_input_data = (double *) nim_input->data;
    } else if(dataType == 256) {
        qDebug() << "Datatype: Signed Char (8 bits)";
//        signed char *nim_input_data = (signed char *) nim_input->data;
    } else if(dataType == 512) {
        qDebug() << "Datatype: Unsigned Short (16 bits)";
//        unsigned short *nim_input_data = (unsigned short *) nim_input->data;
    } else if(dataType == 768) {
        qDebug() << "Datatype: Unsigned Int (32 bits)";
//        unsigned int *nim_input_data = (unsigned int *) nim_input->data;
    }
    //Depending on the nii voxel datatype shoul be chanched.
//    unsigned char *nim_input_data = (unsigned char *) nim_input->data;
//    float  *nim_input_data = (float *) nim_input->data;
//    short int  *nim_input_data = (short int *) nim_input->data;
//    signed short *nim_input_data = (signed short *) nim_input->data;
    unsigned short *nim_input_data = (unsigned short *) nim_input->data;

    int counter = 0;
    for(int timestep=0; timestep<nrep; ++timestep){
        for(int islice=0; islice<sizeSlice; ++islice){
            for(int iy=0; iy<sizeRead; ++iy){
                for(int ix=0; ix<sizePhase; ++ix){
//                    unsigned short tempValu = (unsigned short)inputImg[counter];
//                    float tempValu = (float)inputImg[counter];
//                    signed short tempValu = (signed short)inputImg[counter];         //Depending on the nii voxel datatype shoul be chanched.
                    unsigned short tempValu = (unsigned short)inputImg[counter];         //Depending on the nii voxel datatype shoul be chanched.
//                    unsigned char tempValu = (unsigned char)inputImg[counter];
                    *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix) = tempValu;
//                    nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix] = tempValu;
    //              printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
}

// Copy given 1D short int array entires to given Nifti format object voxel values.
void si_array_to_nii(nifti_image *nim_input, short int* inputImg)
{
    // Get dimensions of input
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    // Get access to data of nim_input
    short int  *nim_input_data = (short int *) nim_input->data;

    int counter = 0;
    for(int timestep=0; timestep<nrep; ++timestep){
        for(int islice=0; islice<sizeSlice; ++islice){
            for(int iy=0; iy<sizeRead; ++iy){
                for(int ix=0; ix<sizePhase; ++ix){
                    short int tempValu = inputImg[counter];
                    *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix) = tempValu;
//                    nim_input_data[nxyz*timestep + nxy*islice + nx*iy + ix] = tempValu;
    //              printf("Voxel values: %f\n", imageArr[counter]);
                    counter++;
                }
            }
        }
    }
}

// Copy NIftI format to 1D uchar array.
unsigned char* nii_to_ucharArray(nifti_image *nim_input)
{
    unsigned char* imageArr;
    double* tempArr;
    int imgWid = nim_input->dim[1], imgHei = nim_input->dim[2], imgDep = nim_input->dim[3], imgTime = nim_input->dim[4];
    short int m=0, M=0;
    int dataType = nim_input->datatype;

    // Memory allocation for imageArr.
    imageArr = (unsigned char*) calloc(imgWid*imgHei*imgDep*imgTime, sizeof(unsigned char));
    if(imageArr == NULL) /* Memory allocation fails */ {
        printf("Couldn't allocate memory for 'imageArr' array\n");
    }
//    else /* Memory allocation successful */ {
//        printf("Memory allocation successful for 'imageArr' array\n");
//    }
    // Memory allocation for tempArr.
    tempArr = (double*) calloc(imgWid*imgHei*imgDep*imgTime, sizeof(double));

    // Get dimensions of input data
    int sizeSlice = nim_input->nz, sizePhase = nim_input->nx, sizeRead = nim_input->ny, nrep =  nim_input->nt;
    int nx =  nim_input->nx, nxy = nim_input->nx * nim_input->ny, nxyz = nim_input->nx * nim_input->ny * nim_input->nz;

    qDebug() << sizePhase << sizeRead << sizeSlice << nrep;
    qDebug() << "Here: " << nim_input->dim[0] << nim_input->dim[1] << nim_input->dim[2] << nim_input->dim[3] << nim_input->dim[4] << nim_input->dim[5];
    qDebug() << "Here1: " << nim_input->datatype << nim_input->intent_code;

    // Get access to data of nim_input
    if(dataType == 0) {
        qDebug() << "Unknown datatype!";
        return 0;
    } else if(dataType == 1) {
        qDebug() << "Datatype: Bool (1 bit)";
//        bool *nim_input_data = (bool *) nim_input->data;
    } else if(dataType == 2) {
        qDebug() << "Datatype: Unsigned Char (8 bits)";
//        unsigned char *nim_input_data = (unsigned char *) nim_input->data;
    } else if(dataType == 4) {
        qDebug() << "Datatype: Signed Short (16 bits)";
//        signed short *nim_input_data = (signed short *) nim_input->data;
    } else if(dataType == 8) {
        qDebug() << "Datatype: Signed Int (32 bits)";
//        signed int *nim_input_data = (signed int *) nim_input->data;
    } else if(dataType == 16) {
        qDebug() << "Datatype: Float (32 bits)";
//        float *nim_input_data = (float *) nim_input->data;
    } else if(dataType == 64) {
        qDebug() << "Datatype: Double (64 bits)";
//        double *nim_input_data = (double *) nim_input->data;
    } else if(dataType == 256) {
        qDebug() << "Datatype: Signed Char (8 bits)";
//        signed char *nim_input_data = (signed char *) nim_input->data;
    } else if(dataType == 512) {
        qDebug() << "Datatype: Unsigned Short (16 bits)";
//        unsigned short *nim_input_data = (unsigned short *) nim_input->data;
    } else if(dataType == 768) {
        qDebug() << "Datatype: Unsigned Int (32 bits)";
//        unsigned int *nim_input_data = (unsigned int *) nim_input->data;
    }
    //Depending on the nii voxel datatype shoul be chanched.
//    unsigned char *nim_input_data = (unsigned char *) nim_input->data;
//        float *nim_input_data = (float *) nim_input->data;
    //    short int *nim_input_data = (short int *) nim_input->data;
    signed short *nim_input_data = (signed short *) nim_input->data;
//    unsigned char *nim_input_data = (unsigned char *) nim_input->data;

    int counter = 0;
    for(int timestep=0; timestep<nrep; ++timestep){
        for(int islice=0; islice<sizeSlice; ++islice){
            for(int iy=0; iy<sizeRead; ++iy){
                for(int ix=0; ix<sizePhase; ++ix){
//                    float tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*ix + iy);
//                    short int tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
                    signed short tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);    //Depending on the nii voxel datatype shoul be chanched.
//                    unsigned short tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);
//                    unsigned char tempValu = *(nim_input_data + nxyz*timestep + nxy*islice + nx*iy + ix);

                    tempArr[counter] = tempValu;
                    if(m > tempValu) {
                        m = tempValu;
                    }
                    if(M < tempValu) {
                        M = tempValu;
                    }
                    counter++;
                }
            }
        }
    }
//    qDebug() << "Hereeeeeeeeeeeeee" << tempArr[33464654];

    for(int i=0; i<counter; i++) {
        M = max(abs(m),abs(M));
        tempArr[i] = 127.5*(tempArr[i]/M + 1);
        imageArr[i] = (unsigned char)tempArr[i];
    }
    delete [] tempArr;
    tempArr = NULL;

    return imageArr;
}

//Scale any double array to uchar([0,255]) array
void scaleToUchar(double* imgArr, int len) {
    double m=0, M=0;

    for(int i=0; i<len; i++) {
        if(m > imgArr[i]) {
            m = imgArr[i];
        }
        if(M < imgArr[i]) {
            M = imgArr[i];
        }
    }
    for(int i=0; i<len; i++) {
        M = max(abs(m),abs(M));
        imgArr[i] = 127.5*(imgArr[i]/M + 1);
//        imgArr[i] = (unsigned char)imgArr[i];
    }
}

//Scale any short int array to uchar([0,255]) array
void si_scaleToUchar(short int* imgArr, int len) {
    double m=0, M=0;

    for(int i=0; i<len; i++) {
        if(m > imgArr[i]) {
            m = imgArr[i];
        }
        if(M < imgArr[i]) {
            M = imgArr[i];
        }
    }
    for(int i=0; i<len; i++) {
        M = max(abs(m),abs(M));
        imgArr[i] = 127.5*(imgArr[i]/M + 1);
//        imgArr[i] = (unsigned char)imgArr[i];
    }
}

// Create 1D Gaussian kernel.
void gaussian1D_kernel(double gKernelX[], int kernelLength, double sigma)
{
    // sum is for normalization
    double sum = 0.0;
    double s = 2.0 * sigma * sigma;
    double r;

    //Initialization input array
    for(int i=0; i<kernelLength; i++) {
        gKernelX[i] = 0;
    }
    // Generate horizontal kernel
    for (int x = -(kernelLength-1)/2; x <= (kernelLength-1)/2; x++)
    {
        r = x*x;
        gKernelX[x + (kernelLength-1)/2] = exp(-r/s)/sqrt(M_PI * s);
        sum += gKernelX[x + (kernelLength-1)/2];
    }
    // Normalize the kernel
    for(int i = 0; i < kernelLength; i++) {
        gKernelX[i] /= sum;
    }
}

//Horizontal convolution.
double* convolutionX(double* inputAr, int inputArWidth, int inputArHeight, double gKernelX[], int kernelLength)
{
//    double* convolvedArr;
    int i, j, k, m, len;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* convolvedArr = new double[sizeof(double) * len];

    // Memory allocation for convolvedArr.
//    convolvedArr = (double*) malloc(sizeof(double) * inputArHeight*inputArWidth);

    // find center position of kernel (half of kernel size)
    kCenter = kernelLength >> 1;                          // center index of kernel array
    endIndex = inputArWidth - kCenter;                    // index for full kernel convolution stop

    // init working indices
    int inputArrIndex = 0;
    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

    // start convolution
    for(i=0; i<inputArHeight; ++i)                        // number of rows
    {
        kOffset = 0;                                      // starting index of partial kernel varies for each sample

        // COLUMN FROM index=0 TO index=kCenter-1
        for(j=0; j < kCenter; ++j)
        {
            convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
            for(k=kCenter+kOffset, m=0; k>=0; --k, ++m)          // convolve with partial of kernel
            {
                convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
            }
            ++convolvedArrIndex;
            ++kOffset;
        }

        // COLUMN FROM index=kCenter TO index=(inputArWidth-kCenter-1)
        for(j = kCenter; j < endIndex; ++j)
        {
            convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate

            for(k = kernelLength-1, m = 0; k >= 0; --k, ++m)  // full kernel
            {
                convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
            }
            ++inputArrIndex;                                // next input
            ++convolvedArrIndex;                               // next output
        }
        kOffset = 1;                                // ending index of partial kernel varies for each sample

        // COLUMN FROM index=(inputArWidth-kCenter) TO index=(inputArWidth-1)
        for(j = endIndex; j < inputArWidth; ++j)
        {
            convolvedArr[convolvedArrIndex] = 0;                            // init to 0 before accumulation

            for(k = kernelLength-1, m=0; k >= kOffset; --k, ++m)   // convolve with partial of kernel
                    {
                        convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                    }
                    ++inputArrIndex;                                // next input
                    ++convolvedArrIndex;                               // next output
                    ++kOffset;                              // increase ending index of partial kernel
                }

                inputArrIndex += kCenter;                           // next row
    }

    return convolvedArr;
}

//Horizontal convolution for 3D image.
double* convolution3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelX[], int kernelLength)
{
//    double* convolvedArr;
    int i, j, k, m, len, z;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* convolvedArr = new double[sizeof(double) * len];

    // find center position of kernel (half of kernel size)
    kCenter = kernelLength >> 1;                          // center index of kernel array
    endIndex = inputArWidth - kCenter;                    // index for full kernel convolution stop

    // init working indices
    int inputArrIndex = 0;
    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

    for(z=0; z<inputArDepth; z++)                             // number of depth, i.e. third z direction
    {
        // start convolution
        for(i=0; i<inputArHeight; ++i)                        // number of rows
        {
            kOffset = 0;                                      // starting index of partial kernel varies for each sample

            // COLUMN FROM index=0 TO index=kCenter-1
            for(j=0; j < kCenter; ++j)
            {
                convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
                for(k=kCenter+kOffset, m=0; k>=0; --k, ++m)          // convolve with partial of kernel
                {
                    convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                }
                ++convolvedArrIndex;
                ++kOffset;
            }

            // COLUMN FROM index=kCenter TO index=(inputArWidth-kCenter-1)
            for(j = kCenter; j < endIndex; ++j)
            {
                convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate

                for(k = kernelLength-1, m = 0; k >= 0; --k, ++m)  // full kernel
                {
                    convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                }
                ++inputArrIndex;                                // next input
                ++convolvedArrIndex;                               // next output
            }
            kOffset = 1;                                // ending index of partial kernel varies for each sample

            // COLUMN FROM index=(inputArWidth-kCenter) TO index=(inputArWidth-1)
            for(j = endIndex; j < inputArWidth; ++j)
            {
                convolvedArr[convolvedArrIndex] = 0;                            // init to 0 before accumulation

                for(k = kernelLength-1, m=0; k >= kOffset; --k, ++m)   // convolve with partial of kernel
                {
                    convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                }
                ++inputArrIndex;                                // next input
                ++convolvedArrIndex;                               // next output
                ++kOffset;                              // increase ending index of partial kernel
            }
            inputArrIndex += kCenter;                           // next row
        }
    }
    return convolvedArr;
}

//Horizontal convolution for 4D Volume.
double* convolution4DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, double gKernelX[], int kernelLength)
{
//    double* convolvedArr;
    int i, j, k, m, len, z;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth*inputArVolume;

    // Decleartion and memory allocation
    double* convolvedArr = new double[sizeof(double) * len];

    // find center position of kernel (half of kernel size)
    kCenter = kernelLength >> 1;                          // center index of kernel array
    endIndex = inputArWidth - kCenter;                    // index for full kernel convolution stop

    // init working indices
    int inputArrIndex = 0;
    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

    for(z=0; z<inputArDepth; z++)                             // number of depth, i.e. third z direction
    {
        // start convolution
        for(i=0; i<inputArHeight; ++i)                        // number of rows
        {
            kOffset = 0;                                      // starting index of partial kernel varies for each sample

            // COLUMN FROM index=0 TO index=kCenter-1
            for(j=0; j < kCenter; ++j)
            {
                convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
                for(k=kCenter+kOffset, m=0; k>=0; --k, ++m)          // convolve with partial of kernel
                {
                    convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                }
                ++convolvedArrIndex;
                ++kOffset;
            }

            // COLUMN FROM index=kCenter TO index=(inputArWidth-kCenter-1)
            for(j = kCenter; j < endIndex; ++j)
            {
                convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate

                for(k = kernelLength-1, m = 0; k >= 0; --k, ++m)  // full kernel
                {
                    convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                }
                ++inputArrIndex;                                // next input
                ++convolvedArrIndex;                               // next output
            }
            kOffset = 1;                                // ending index of partial kernel varies for each sample

            // COLUMN FROM index=(inputArWidth-kCenter) TO index=(inputArWidth-1)
            for(j = endIndex; j < inputArWidth; ++j)
            {
                convolvedArr[convolvedArrIndex] = 0;                            // init to 0 before accumulation

                for(k = kernelLength-1, m=0; k >= kOffset; --k, ++m)   // convolve with partial of kernel
                {
                    convolvedArr[convolvedArrIndex] += inputAr[inputArrIndex+m] * gKernelX[k];
                }
                ++inputArrIndex;                                // next input
                ++convolvedArrIndex;                               // next output
                ++kOffset;                              // increase ending index of partial kernel
            }
            inputArrIndex += kCenter;                           // next row
        }
    }
    return convolvedArr;
}

//Horizontal convolution for 3D image series of 4D data
double** convolution3DX_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelX[], int kernelLength) {
    int i, j, k, m, len, z;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** convolvedArr = new double*[sizeof(double)*timeLen];
    for(int i = 0; i < timeLen; i++) {
        convolvedArr[i] = new double[sizeof(double)*len];
    }

    for(int t=0; t<timeLen; t++) {
        // find center position of kernel (half of kernel size)
        kCenter = kernelLength >> 1;                          // center index of kernel array
        endIndex = inputArWidth - kCenter;                    // index for full kernel convolution stop

        // init working indices
        int inputArrIndex = 0;
        int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

        for(z=0; z<inputArDepth; z++)                             // number of depth, i.e. third z direction
        {
            // start convolution
            for(i=0; i<inputArHeight; ++i)                        // number of rows
            {
                kOffset = 0;                                      // starting index of partial kernel varies for each sample

                // COLUMN FROM index=0 TO index=kCenter-1
                for(j=0; j < kCenter; ++j)
                {
                    convolvedArr[t][convolvedArrIndex] = 0;                              // init to 0 before accumulate
                    for(k=kCenter+kOffset, m=0; k>=0; --k, ++m)          // convolve with partial of kernel
                    {
                        convolvedArr[t][convolvedArrIndex] += inputAr[t][inputArrIndex+m] * gKernelX[k];
                    }
                    ++convolvedArrIndex;
                    ++kOffset;
                }

                // COLUMN FROM index=kCenter TO index=(inputArWidth-kCenter-1)
                for(j = kCenter; j < endIndex; ++j)
                {
                    convolvedArr[t][convolvedArrIndex] = 0;                              // init to 0 before accumulate

                    for(k = kernelLength-1, m = 0; k >= 0; --k, ++m)  // full kernel
                    {
                        convolvedArr[t][convolvedArrIndex] += inputAr[t][inputArrIndex+m] * gKernelX[k];
                    }
                    ++inputArrIndex;                                // next input
                    ++convolvedArrIndex;                               // next output
                }
                kOffset = 1;                                // ending index of partial kernel varies for each sample

                // COLUMN FROM index=(inputArWidth-kCenter) TO index=(inputArWidth-1)
                for(j = endIndex; j < inputArWidth; ++j)
                {
                    convolvedArr[t][convolvedArrIndex] = 0;                            // init to 0 before accumulation

                    for(k = kernelLength-1, m=0; k >= kOffset; --k, ++m)   // convolve with partial of kernel
                    {
                        convolvedArr[t][convolvedArrIndex] += inputAr[t][inputArrIndex+m] * gKernelX[k];
                    }
                    ++inputArrIndex;                                // next input
                    ++convolvedArrIndex;                               // next output
                    ++kOffset;                              // increase ending index of partial kernel
                }
                inputArrIndex += kCenter;                           // next row
            }
        }
    }
    return convolvedArr;
}

//Vertical convolution.
double* convolutionY(double* inputAr, int inputArWidth, int inputArHeight, double gKernelY[], int kernelLength){

//    double* convolvedArr;
    double sum[inputArWidth];
    int i, j, k, n, len;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* convolvedArr = new double[sizeof(double) * len];

    // clear out array before accumulation
    for(i = 0; i < inputArWidth; ++i)
        sum[i] = 0;

    // find center position of kernel (half of kernel size)
    kCenter = kernelLength >> 1;                          // center index of vertical kernel
    endIndex = inputArHeight - kCenter;                 // index where full kernel convolution should stop

    // init working indices
    int inputArrIndex = 0;
    int currentRowNumber = 0;
//    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

    // ROW FROM index=0 TO index=(kCenter-1)
    kOffset = 0;                                    // starting index of partial kernel varies for each sample
    for(i=0; i < kCenter; ++i)
    {
//        convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
        for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel
        {
            for(j=0; j < inputArWidth; ++j)
            {
                sum[j] += inputAr[inputArrIndex] * gKernelY[k];
                ++inputArrIndex;
            }
        }
        for(n=0; n<inputArWidth; ++n)
        {
            convolvedArr[n+i*inputArWidth] = sum[n];
            sum[n] = 0;
        }

        inputArrIndex = currentRowNumber;                           // reset input array index
        ++kOffset;                                  // increase starting index of kernel
    }
    // ROW FROM index=kCenter TO index=(inputArHeight-kCenter-1)
    for(i = kCenter; i < endIndex; ++i)
    {
        for(k=kernelLength-1; k>=0; --k)             // convolve with full kernel
        {
            for(j = 0; j < inputArWidth; ++j)
            {
                sum[j] += inputAr[inputArrIndex] * gKernelY[k];
                ++inputArrIndex;
            }
        }

        for(n=0; n<inputArWidth; ++n)
        {
            convolvedArr[n+i*inputArWidth] = sum[n];
            sum[n] = 0;
        }

        currentRowNumber += inputArWidth;                // move to next row
        inputArrIndex = currentRowNumber;                           // reset input array index
    }

    // ROW FROM index=(inputArHeight-kCenter) TO index=(inputArHeight-1)
    kOffset = 1;                                    // ending index of partial kernel varies for each sample
    for(i=endIndex; i < inputArHeight; ++i)
    {
        for(k = kernelLength-1; k >= kOffset; --k)        // convolve with partial kernel
        {
            for(j=0; j < inputArWidth; ++j)
            {
                sum[j] += inputAr[inputArrIndex] * gKernelY[k];
                ++inputArrIndex;
            }
        }

        for(n = 0; n < inputArWidth; ++n)              // convert and copy from sum to out
        {
            convolvedArr[n+i*inputArWidth] = sum[n];
            sum[n] = 0;
        }

        currentRowNumber += inputArWidth;                // move to next row
        inputArrIndex = currentRowNumber;              // move to next row
        ++kOffset;                                  // increase ending index of kernel
    }

    return convolvedArr;
}

//Vertical convolution for 3D image.
double* convolution3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength) {
//    double* convolvedArr;
    double sum[inputArWidth];
    int i, j, k, n, len, z;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* convolvedArr = new double[sizeof(double) * len];

    // clear out array before accumulation
    for(i = 0; i < inputArWidth; ++i)
        sum[i] = 0;

    // find center position of kernel (half of kernel size)
    kCenter = kernelLength >> 1;                          // center index of vertical kernel
    endIndex = inputArHeight - kCenter;                 // index where full kernel convolution should stop

//    // init working indices
//    int inputArrIndex = 0;
//    int currentRowNumber = 0;

    for(z=0; z<inputArDepth; z++)
    {
        // init working indices
        int inputArrIndex = 0;
        int currentRowNumber = 0;
    //    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

        // ROW FROM index=0 TO index=(kCenter-1)
        kOffset = 0;                                    // starting index of partial kernel varies for each sample
        for(i=0; i < kCenter; ++i)
        {
    //        convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
            for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel
            {
                for(j=0; j < inputArWidth; ++j)
                {
                    sum[j] += inputAr[inputArrIndex + z*inputArWidth*inputArHeight] * gKernelY[k];
                    ++inputArrIndex;
                }
            }
            for(n=0; n<inputArWidth; ++n)
            {
                convolvedArr[n+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                sum[n] = 0;
            }

            inputArrIndex = currentRowNumber;                           // reset input array index
            ++kOffset;                                                  // increase starting index of kernel
        }
        // ROW FROM index=kCenter TO index=(inputArHeight-kCenter-1)
        for(i = kCenter; i < endIndex; ++i)
        {
            for(k=kernelLength-1; k>=0; --k)             // convolve with full kernel
            {
                for(j = 0; j < inputArWidth; ++j)
                {
                    sum[j] += inputAr[inputArrIndex + z*inputArWidth*inputArHeight] * gKernelY[k];
                    ++inputArrIndex;
                }
            }

            for(n=0; n<inputArWidth; ++n)
            {
                convolvedArr[n+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                sum[n] = 0;
            }

            currentRowNumber += inputArWidth;                // move to next row
            inputArrIndex = currentRowNumber;                           // reset input array index
        }
        // ROW FROM index=(inputArHeight-kCenter) TO index=(inputArHeight-1)
        kOffset = 1;                                    // ending index of partial kernel varies for each sample
        for(i=endIndex; i < inputArHeight; ++i)
        {
            for(k = kernelLength-1; k >= kOffset; --k)        // convolve with partial kernel
            {
                for(j=0; j < inputArWidth; ++j)
                {
                    sum[j] += inputAr[inputArrIndex + z*inputArWidth*inputArHeight] * gKernelY[k];
                    ++inputArrIndex;
                }
            }

            for(n = 0; n < inputArWidth; ++n)              // convert and copy from sum to out
            {
                convolvedArr[n+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                sum[n] = 0;
            }

            currentRowNumber += inputArWidth;                // move to next row
            inputArrIndex = currentRowNumber;              // move to next row
            ++kOffset;                                  // increase ending index of kernel
        }
    }
    return convolvedArr;
}

//Horizontal convolution for 3D image series of 4D data
double** convolution3DY_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength) {
    double sum[inputArWidth];
    int i, j, k, n, len, z;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** convolvedArr = new double*[sizeof(double)*timeLen];
    for(int i = 0; i < timeLen; i++) {
        convolvedArr[i] = new double[sizeof(double)*len];
    }

    for(int t=0; t<timeLen; t++) {
        // clear out array before accumulation
        for(i = 0; i < inputArWidth; ++i)
            sum[i] = 0;

        // find center position of kernel (half of kernel size)
        kCenter = kernelLength >> 1;                          // center index of kernel array
        endIndex = inputArWidth - kCenter;                    // index for full kernel convolution stop

//        // init working indices
//        int inputArrIndex = 0;
//        int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

        for(z=0; z<inputArDepth; z++) {
            // init working indices
            int inputArrIndex = 0;
            int currentRowNumber = 0;
        //    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

            // ROW FROM index=0 TO index=(kCenter-1)
            kOffset = 0;                                    // starting index of partial kernel varies for each sample
            for(i=0; i < kCenter; ++i)
            {
        //        convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
                for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel
                {
                    for(j=0; j < inputArWidth; ++j)
                    {
                        sum[j] += inputAr[t][inputArrIndex + z*inputArWidth*inputArHeight] * gKernelY[k];
                        ++inputArrIndex;
                    }
                }
                for(n=0; n<inputArWidth; ++n)
                {
                    convolvedArr[t][n+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                    sum[n] = 0;
                }

                inputArrIndex = currentRowNumber;                           // reset input array index
                ++kOffset;                                                  // increase starting index of kernel
            }
            // ROW FROM index=kCenter TO index=(inputArHeight-kCenter-1)
            for(i = kCenter; i < endIndex; ++i)
            {
                for(k=kernelLength-1; k>=0; --k)             // convolve with full kernel
                {
                    for(j = 0; j < inputArWidth; ++j)
                    {
                        sum[j] += inputAr[t][inputArrIndex + z*inputArWidth*inputArHeight] * gKernelY[k];
                        ++inputArrIndex;
                    }
                }

                for(n=0; n<inputArWidth; ++n)
                {
                    convolvedArr[t][n+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                    sum[n] = 0;
                }

                currentRowNumber += inputArWidth;                // move to next row
                inputArrIndex = currentRowNumber;                           // reset input array index
            }
            // ROW FROM index=(inputArHeight-kCenter) TO index=(inputArHeight-1)
            kOffset = 1;                                    // ending index of partial kernel varies for each sample
            for(i=endIndex; i < inputArHeight; ++i)
            {
                for(k = kernelLength-1; k >= kOffset; --k)        // convolve with partial kernel
                {
                    for(j=0; j < inputArWidth; ++j)
                    {
                        sum[j] += inputAr[t][inputArrIndex + z*inputArWidth*inputArHeight] * gKernelY[k];
                        ++inputArrIndex;
                    }
                }

                for(n = 0; n < inputArWidth; ++n)              // convert and copy from sum to out
                {
                    convolvedArr[t][n+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                    sum[n] = 0;
                }

                currentRowNumber += inputArWidth;                // move to next row
                inputArrIndex = currentRowNumber;              // move to next row
                ++kOffset;                                  // increase ending index of kernel
            }
        }
    }
    return convolvedArr;
}

//Z direction convolution for 3D image.
double* convolution3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength) {

//    double* convolvedArr;
    double sum[inputArWidth];
    int i, j, k, n, len, h;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* convolvedArr = new double[sizeof(double) * len];

    // clear out array before accumulation
    for(i = 0; i < inputArWidth; ++i)
        sum[i] = 0;

    // find center position of kernel (half of kernel size)
    kCenter = kernelLength >> 1;                          // center index of vertical kernel
    endIndex = inputArDepth - kCenter;                 // index where full kernel convolution should stop

    for(h=0; h<inputArHeight; h++)
    {
        // init working indices
        int inputArrIndex = 0;
        int currentDepthNumber = 0;
    //    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

        // DEPTH FROM 0 TO (kCenter-1)
        kOffset = 0;                                    // starting index of partial kernel varies for each sample
        for(i=0; i < kCenter; ++i)
        {
    //        convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
            for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel
            {
                for(j=0; j < inputArWidth; ++j)
                {
                    sum[j] += inputAr[inputArrIndex + h*inputArWidth] * gKernelY[k];
                    ++inputArrIndex;
                }
                inputArrIndex = inputArrIndex - inputArWidth + inputArWidth*inputArHeight;
            }
            for(n=0; n<inputArWidth; ++n)
            {
                convolvedArr[n+i*inputArWidth*inputArHeight + h*inputArWidth] = sum[n];
                sum[n] = 0;
            }

            inputArrIndex = currentDepthNumber;                           // reset input array index
            ++kOffset;                                                  // increase starting index of kernel
        }
        // ROW FROM index=kCenter TO index=(inputArHeight-kCenter-1)
        for(i = kCenter; i < endIndex; ++i)
        {
            for(k=kernelLength-1; k>=0; --k)             // convolve with full kernel
            {
                for(j = 0; j < inputArWidth; ++j)
                {
                    sum[j] += inputAr[inputArrIndex + h*inputArWidth] * gKernelY[k];
                    ++inputArrIndex;
                }
                inputArrIndex = inputArrIndex - inputArWidth + inputArWidth*inputArHeight;
            }

            for(n=0; n<inputArWidth; ++n)
            {
                convolvedArr[n+i*inputArWidth*inputArHeight + h*inputArWidth] = sum[n];
                sum[n] = 0;
            }

            currentDepthNumber += inputArWidth*inputArHeight;                // move to next row
            inputArrIndex = currentDepthNumber;                           // reset input array index
        }
        // ROW FROM index=(inputArHeight-kCenter) TO index=(inputArHeight-1)
        kOffset = 1;                                    // ending index of partial kernel varies for each sample
        for(i=endIndex; i < inputArHeight; ++i)
        {
            for(k = kernelLength-1; k >= kOffset; --k)        // convolve with partial kernel
            {
                for(j=0; j < inputArWidth; ++j)
                {
                    sum[j] += inputAr[inputArrIndex + h*inputArWidth] * gKernelY[k];
                    ++inputArrIndex;
                }
                inputArrIndex = inputArrIndex - inputArWidth + inputArWidth*inputArHeight;
            }

            for(n = 0; n < inputArWidth; ++n)              // convert and copy from sum to out
            {
                convolvedArr[n+i*inputArWidth*inputArHeight + h*inputArWidth] = sum[n];
                sum[n] = 0;
            }

            currentDepthNumber += inputArWidth*inputArHeight;                // move to next row
            inputArrIndex = currentDepthNumber;              // move to next row
            ++kOffset;                                  // increase ending index of kernel
        }
    }

    return convolvedArr;
}

//Z direction convolution for 3D image.
double** convolution3DZ_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength) {
    double sum[inputArWidth];
    int i, j, k, n, len, h;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** convolvedArr = new double*[sizeof(double)*timeLen];
    for(int i = 0; i < timeLen; i++) {
        convolvedArr[i] = new double[sizeof(double)*len];
    }

    for(int t=0; t<timeLen; t++) {
        // clear out array before accumulation
        for(i = 0; i < inputArWidth; ++i)
            sum[i] = 0;

        // find center position of kernel (half of kernel size)
        kCenter = kernelLength >> 1;                          // center index of vertical kernel
        endIndex = inputArDepth - kCenter;                    // index where full kernel convolution should stop

        for(h=0; h<inputArHeight; h++) {
            // init working indices
            int inputArrIndex = 0;
            int currentDepthNumber = 0;
        //    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

            // DEPTH FROM 0 TO (kCenter-1)
            kOffset = 0;                                    // starting index of partial kernel varies for each sample
            for(i=0; i < kCenter; ++i)
            {
        //        convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
                for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel
                {
                    for(j=0; j < inputArWidth; ++j)
                    {
                        sum[j] += inputAr[t][inputArrIndex + h*inputArWidth] * gKernelY[k];
                        ++inputArrIndex;
                    }
                    inputArrIndex = inputArrIndex - inputArWidth + inputArWidth*inputArHeight;
                }
                for(n=0; n<inputArWidth; ++n)
                {
                    convolvedArr[t][n+i*inputArWidth*inputArHeight + h*inputArWidth] = sum[n];
                    sum[n] = 0;
                }

                inputArrIndex = currentDepthNumber;                           // reset input array index
                ++kOffset;                                                  // increase starting index of kernel
            }
            // ROW FROM index=kCenter TO index=(inputArHeight-kCenter-1)
            for(i = kCenter; i < endIndex; ++i)
            {
                for(k=kernelLength-1; k>=0; --k)             // convolve with full kernel
                {
                    for(j = 0; j < inputArWidth; ++j)
                    {
                        sum[j] += inputAr[t][inputArrIndex + h*inputArWidth] * gKernelY[k];
                        ++inputArrIndex;
                    }
                    inputArrIndex = inputArrIndex - inputArWidth + inputArWidth*inputArHeight;
                }

                for(n=0; n<inputArWidth; ++n)
                {
                    convolvedArr[t][n+i*inputArWidth*inputArHeight + h*inputArWidth] = sum[n];
                    sum[n] = 0;
                }

                currentDepthNumber += inputArWidth*inputArHeight;                // move to next row
                inputArrIndex = currentDepthNumber;                           // reset input array index
            }
            // ROW FROM index=(inputArHeight-kCenter) TO index=(inputArHeight-1)
            kOffset = 1;                                    // ending index of partial kernel varies for each sample
            for(i=endIndex; i < inputArHeight; ++i)
            {
                for(k = kernelLength-1; k >= kOffset; --k)        // convolve with partial kernel
                {
                    for(j=0; j < inputArWidth; ++j)
                    {
                        sum[j] += inputAr[t][inputArrIndex + h*inputArWidth] * gKernelY[k];
                        ++inputArrIndex;
                    }
                    inputArrIndex = inputArrIndex - inputArWidth + inputArWidth*inputArHeight;
                }

                for(n = 0; n < inputArWidth; ++n)              // convert and copy from sum to out
                {
                    convolvedArr[t][n+i*inputArWidth*inputArHeight + h*inputArWidth] = sum[n];
                    sum[n] = 0;
                }

                currentDepthNumber += inputArWidth*inputArHeight;                // move to next row
                inputArrIndex = currentDepthNumber;              // move to next row
                ++kOffset;                                  // increase ending index of kernel
            }
        }
    }
    return convolvedArr;
}

//T direction convolution for 3D image.
double** convolution3DT_of_4D(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth, double gKernelY[], int kernelLength) {
    double sum[inputArWidth];
    int i, j, k, n, len, z;
    int kCenter, kOffset, endIndex;                 // kernel indice

    len  = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** convolvedArr = new double*[sizeof(double)*timeLen];
    for(int i = 0; i < timeLen; i++) {
        convolvedArr[i] = new double[sizeof(double)*len];
    }

    for(int t=0; t<inputArWidth; t++) {
        // clear out array before accumulation
        for(i = 0; i < inputArWidth; ++i)
            sum[i] = 0;

        // find center position of kernel (half of kernel size)
        kCenter = kernelLength >> 1;                          // center index of kernel array
        endIndex = inputArWidth - kCenter;                    // index for full kernel convolution stop

//        // init working indices
//        int inputArrIndex = 0;
//        int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

        for(z=0; z<inputArDepth; z++) {
            // init working indices
            int inputArrIndex = 0;
            int currentVolNumber = 0;
        //    int convolvedArrIndex = 0;                                   // store intermediate results from 1D horizontal convolution

            // ROW FROM index=0 TO index=(kCenter-1)
            kOffset = 0;                                    // starting index of partial kernel varies for each sample
            for(i=0; i < kCenter; ++i)
            {
        //        convolvedArr[convolvedArrIndex] = 0;                              // init to 0 before accumulate
                for(k = kCenter + kOffset; k >= 0; --k)     // convolve with partial kernel
                {
                    for(j=0; j < timeLen; ++j)
                    {
                        sum[j] += inputAr[inputArrIndex][t + z*inputArWidth*inputArHeight] * gKernelY[k];
                        ++inputArrIndex;
                    }
                }
                for(n=0; n<timeLen; ++n)
                {
                    convolvedArr[n][t+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                    sum[n] = 0;
                }

                inputArrIndex = currentVolNumber;                           // reset input array index
                ++kOffset;                                                  // increase starting index of kernel
            }
            // ROW FROM index=kCenter TO index=(inputArHeight-kCenter-1)
            for(i = kCenter; i < endIndex; ++i)
            {
                for(k=kernelLength-1; k>=0; --k)             // convolve with full kernel
                {
                    for(j = 0; j < timeLen; ++j)
                    {
                        sum[j] += inputAr[inputArrIndex][t + z*inputArWidth*inputArHeight] * gKernelY[k];
                        ++inputArrIndex;
                    }
                }

                for(n=0; n<timeLen; ++n)
                {
                    convolvedArr[n][t+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                    sum[n] = 0;
                }

                currentVolNumber += inputArWidth;                // move to next row
                inputArrIndex = currentVolNumber;                           // reset input array index
            }
            // ROW FROM index=(inputArHeight-kCenter) TO index=(inputArHeight-1)
            kOffset = 1;                                    // ending index of partial kernel varies for each sample
            for(i=endIndex; i < inputArHeight; ++i)
            {
                for(k = kernelLength-1; k >= kOffset; --k)        // convolve with partial kernel
                {
                    for(j=0; j < timeLen; ++j)
                    {
                        sum[j] += inputAr[j][t + z*inputArWidth*inputArHeight] * gKernelY[k];
                        ++inputArrIndex;
                    }
                }

                for(n = 0; n < timeLen; ++n)              // convert and copy from sum to out
                {
                    convolvedArr[n][t+i*inputArWidth + z*inputArWidth*inputArHeight] = sum[n];
                    sum[n] = 0;
                }

                currentVolNumber += inputArWidth;                // move to next row
                inputArrIndex = currentVolNumber;              // move to next row
                ++kOffset;                                  // increase ending index of kernel
            }
        }
    }
    return convolvedArr;
}


//Mean Squared Error
double mse(double* imageArray1, double* imageArray2, int imageArrayLen)
{
    double sum_sq = 0;
    double mse;

    for (int i = 0; i < imageArrayLen; ++i)
    {
        double err = imageArray1[i] - imageArray2[i];
        sum_sq += (err * err);
    }
    mse = sum_sq/imageArrayLen;

    return mse;
}

//Mean Squared Error
double mse_4D(double** imageArray1, double** imageArray2, int imageArrayLen, int imgVolume)
{
    double sum_sq = 0;
    double mse;

    for (int v=0;v< imgVolume; v++) {
        for (int i = 0; i < imageArrayLen; ++i)
        {
            double err = imageArray1[i] - imageArray2[i];
            sum_sq += (err * err);
        }
    }
    mse = sum_sq/imageArrayLen*imgVolume;

    return mse;
}

//Mean Squared Error for Dilated Mask Region
double mse_mask_dilatied_region3D(double* imageArray1, double* imageArray2, int imgWidth, int imgHeight, int imgDepth, int* mskPxlLocs)
{
    int dilMskLen=0;
    double sum_sq = 0, mse;
    double *dilatedPxls;

    dilatedPxls = dilatedMask3D(imgWidth, imgHeight, imgDepth, mskPxlLocs, dilMskLen, true);

    for (int i=0; i <dilMskLen; ++i) {
        double err = imageArray1[i] - imageArray2[i];
        sum_sq += (err * err);
    }
    mse = sum_sq/dilMskLen;

    return mse;
}

//Mean Squared Error for Dilated Mask Region
double mse_mask_dilatied_region4D(double** imageArray1, double** imageArray2, int imgWidth, int imgHeight, int imgDepth, int imgVolume, int* mskPxlLocs)
{
    int dilMskLen=0;
    double sum_sq = 0, mse, sum_mse=0;
    double *dilatedPxls;

    for (int v=0;v<imgVolume;v++) {
    dilatedPxls = dilatedMask3D(imgWidth, imgHeight, imgDepth, mskPxlLocs, dilMskLen, true);
        for (int i=0; i <dilMskLen; ++i) {
            double err = imageArray1[i] - imageArray2[i];
            sum_sq += (err * err);
        }
    sum_mse += sum_sq/dilMskLen;
    }
    mse = sum_mse/imgVolume;

    return mse;
}

//Mean Squared Error of Sub Array (Region)
double sub_mse(double* imageArray1, double* imageArray2, double* maskArr, int imageArrayLen)
{
    int maskLen = 0;
    double sum_sq = 0;
    double mse;

    for (int i = 0; i < imageArrayLen; ++i)
    {
        if(maskArr[i] != 0.0) {
            double err = imageArray1[i] - imageArray2[i];
            sum_sq += (err * err);
            maskLen++;
        }
    }
    mse = sum_sq/maskLen;
    return mse;
}

//Mean Squared Error of Sub Array (Region) for 4D (2D(time+data) array) data
double sub_mse_4D(double** imageArray1, double** imageArray2, double** maskArr, int timeLen, int imageArrayLen)
{
    int maskLen = 0;
    double sum_sq = 0;
    double mse;

    for(int t=0; t<timeLen; t++) {
        for (int i = 0; i < imageArrayLen; ++i)
        {
            if(maskArr[t][i] != 0.0) {
                double err = imageArray1[t][i] - imageArray2[t][i];
                sum_sq += (err * err);
                maskLen++;
            }
        }
    }
    mse = sum_sq/maskLen;
    return mse;
}

//Absolute Average Error
double aae(double* imageArray1, double* imageArray2, int imageArrayLen)
{
    double sum_abs = 0;
    double aae;

    for (int i = 0; i < imageArrayLen; ++i)
    {
        double err = imageArray1[i] - imageArray2[i];
        sum_abs += abs(err);
    }
    aae = sum_abs/imageArrayLen;

    return aae;
}

//Absolute Average Error
double aae_4D(double** imageArray1, double** imageArray2, int imageArrayLen, int imgVolume)
{
    double sum_abs = 0;
    double aae;

    for (int v = 0; v<imgVolume; v++) {
        for (int i = 0; i < imageArrayLen; ++i)
        {
            double err = imageArray1[i] - imageArray2[i];
            sum_abs += abs(err);
        }
    }
    aae = sum_abs/imageArrayLen*imgVolume;

    return aae;
}


//Absolute Average Error for Dilated Mask Region
double aae_mask_dilatied_region3D(double* imageArray1, double* imageArray2, int imgWidth, int imgHeight, int imgDepth, int* mskPxlLocs)
{
    int dilMskLen=0;
    double sum_abs = 0, aae;
    double *dilatedPxls;

    dilatedPxls = dilatedMask3D(imgWidth, imgHeight, imgDepth, mskPxlLocs, dilMskLen, true);

    for (int i=0; i <dilMskLen; ++i) {
        double err = imageArray1[i] - imageArray2[i];
        sum_abs += abs(err);
    }
    aae = sum_abs/dilMskLen;

    return aae;
}

double aae_mask_dilatied_region4D(double** imageArray1, double** imageArray2, int imgWidth, int imgHeight, int imgDepth, int imgVolume, int* mskPxlLocs)
{
    int dilMskLen=0;
    double sum_abs = 0, aae, sum_aae=0;
    double *dilatedPxls;

    for (int v = 0; v<imgVolume; v++) {
        dilatedPxls = dilatedMask3D(imgWidth, imgHeight, imgDepth, mskPxlLocs, dilMskLen, true);

        for (int i=0; i <dilMskLen; ++i) {
            double err = imageArray1[i] - imageArray2[i];
            sum_abs += abs(err);
        }
        sum_aae += sum_abs/dilMskLen;
    }
    aae = sum_aae/imgVolume;


    return aae;
}

//Absolute Average Error of Sub Array (Region)
double sub_aae(double* imageArray1, double* imageArray2, double* maskArr, int imageArrayLen)
{
    int maskLen = 0;
    double sum_abs = 0;
    double aae;

    for (int i = 0; i < imageArrayLen; ++i)
    {
        if(maskArr[i] != 0.0) {
            double err = imageArray1[i] - imageArray2[i];
            sum_abs += abs(err);
            maskLen++;
        }
    }
    aae = sum_abs/maskLen;
    return aae;
}

//Absolute Average Error of Sub Array (Region) for 4D (2D(time+data) array) data
double sub_aae_4D(double** imageArray1, double** imageArray2, double** maskArr, int timeLen, int imageArrayLen)
{
    int maskLen = 0;
    double sum_abs = 0;
    double aae;

    for(int t=0; t<timeLen; t++) {
        for (int i = 0; i < imageArrayLen; ++i)
        {
            if(maskArr[t][i] != 0.0) {
                double err = imageArray1[t][i] - imageArray2[t][i];
                sum_abs += abs(err);
                maskLen++;
            }
        }
    }
    aae = sum_abs/maskLen;
    return aae;
}

//Discrete l2 norm
double l2Norm(double* imageArray1, double* imageArray2, int imageArrayLen) {
    double sum_sq = 0;
    double l2norm;

    for (int i = 0; i < imageArrayLen; ++i)
    {
        double err = imageArray1[i] - imageArray2[i];

        sum_sq += (err * err);
    }
    l2norm = sqrt(sum_sq);

    return l2norm;
}

//Discrete l2 norm for 4D (2D(time+data) array) data
double l2Norm_4D(double** imageArray1, double** imageArray2, int timeLen, int imageArrayLen) {
    double sum_sq = 0;
    double l2norm;

    for(int t=0; t<timeLen; t++) {
        for (int i = 0; i < imageArrayLen; ++i) {
            double err = imageArray1[t][i] - imageArray2[t][i];
            sum_sq += (err * err);
        }
    }
    l2norm = sqrt(sum_sq);

    return l2norm;
}

//Charbonnier diffusivity function
double charbonnier_diff(double s_square, double contrastParam) {
    double g;

    g = 1.0/sqrt(1 + s_square/(contrastParam*contrastParam));

    return g;
}

//Perona-Malik diffusivity function
double pm_diff(double s_square, double contrastParam) {
    double g;

    g = (double)1.0/(1.0 + (s_square/(contrastParam*contrastParam)));

    return g;
}

//Aubert
double aubert_diff(double s_square, double contrastParam) {
    double g;

    g = (1.0)/(contrastParam*(contrastParam*contrastParam + s_square));

    return g;
}

//Green
double green_diff(double s, double contPar) {
    double g;

    if(s == 0) {
        s = 0.001;
    }
    g = tanh(s/contPar)/(s*contPar);

    return g;
}

//Geman et Reynolds
double gr_diff(double s_square, double contrastParam) {
    double g;

    g = (2*contrastParam*contrastParam)/((contrastParam*contrastParam + s_square)*(contrastParam*contrastParam + s_square));

    return g;
}

//Perona-Malik diffusivity function 2
double pm_diff2(double s_square, double contrastParam) {
    double g;

    g = exp(-(s_square/(contrastParam*contrastParam)));

    return g;
}

double li1(double s, double contPar) {
    double g;

    g = 1.0/(1+(s/contPar));

    return g;
}

//Forward sum on x axis.
double* sumForwX(double* inputAr, int inputArWidth, int inputArHeight)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len; i++) {
        if((i+1) % inputArWidth == 0) {
            outAr[i] = inputAr[i];
        } else {
            outAr[i] = inputAr[i+1]+inputAr[i];
        }
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;
}

//Backward sum on x axis.
double* sumBackX(double* inputAr, int inputArWidth, int inputArHeight)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len; i++)
    {
        if(i % inputArWidth == 0) {
            outAr[i] = inputAr[i];
        } else {
            outAr[i] = inputAr[i-1]+inputAr[i];
        }
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;
}

//Forward sum on y axis.
double* sumForwY(double* inputAr, int inputArWidth, int inputArHeight)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len-inputArWidth; i++)
        outAr[i] = inputAr[i+inputArWidth]+inputAr[i];
    for(int i=inputArWidth*(inputArHeight-1); i<inputArWidth*inputArHeight; i++)
        outAr[i] = inputAr[i];
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;
}

//Backward derivative on y axis.
double* sumBackY(double* inputAr, int inputArWidth, int inputArHeight)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<inputArWidth; i++)
        outAr[i] = inputAr[i];
    for(int i=inputArWidth; i<len; i++)
        outAr[i] = inputAr[i]+inputAr[i-inputArWidth];
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;
}

// 3D Forward sum on x axis.
double* sumForw3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<(d+1)*slen; i++) {
            if((i+1) % inputArWidth == 0) {
                outAr[i] = inputAr[i];
            } else {
                outAr[i] = inputAr[i+1]+inputAr[i];
            }
        }
    }
    return outAr;
}

// 4D Forward sum on x axis.
double* *sumForw4DX(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }
    for (int v = 0; v < inputArVolume; v++) {
        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<(d+1)*slen; i++) {
                if((i+1) % inputArWidth == 0) {
                    outAr[v][i] = inputAr[v][i];
                } else {
                    outAr[v][i] = inputAr[v][i+1]+inputAr[v][i];
                }
            }
        }
    }
    return outAr;
}

//Backward sum on x axis.
double* sumBack3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<(d+1)*slen; i++) {
            if(i % inputArWidth == 0) {
                outAr[i] = inputAr[i];
            } else {
                outAr[i] = inputAr[i-1]+inputAr[i];
            }
        }
    }
    return outAr;
}

//Backward sum on x axis.
double* *sumBack4DX(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }

    for (int v = 0; v < inputArVolume;v++) {

        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<(d+1)*slen; i++) {
                if(i % inputArWidth == 0) {
                    outAr[v][i] = inputAr[v][i];
                } else {
                    outAr[v][i] = inputAr[v][i-1]+inputAr[v][i];
                }
            }
        }
    }
    return outAr;
}


//3D Forward sum on y axis.
double* sumForw3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<slen-inputArWidth+d*slen; i++)
            outAr[i] = inputAr[i+inputArWidth]+inputAr[i];
        for(int i=inputArWidth*(inputArHeight-1)+d*slen; i<(d+1)*slen; i++)
            outAr[i] = inputAr[i];
    }
    return outAr;
}

//4D Forward sum on y axis.
double* *sumForw4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }

    for (int v = 0; v < inputArVolume; v++) {
        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<slen-inputArWidth+d*slen; i++)
                outAr[v][i] = inputAr[v][i+inputArWidth]+inputAr[v][i];
            for(int i=inputArWidth*(inputArHeight-1)+d*slen; i<(d+1)*slen; i++)
                outAr[v][i] = inputAr[v][i];
        }
    }
    return outAr;
}

//3D Backward derivative on y axis.
double* sumBack3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<d*slen+inputArWidth; i++)
            outAr[i] = inputAr[i];
        for(int i=inputArWidth+d*slen; i<(d+1)*slen; i++)
            outAr[i] = inputAr[i]+inputAr[i-inputArWidth];
    }

    return outAr;
}

//4D Backward derivative on y axis.
double* *sumBack4DY(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }
    for (int v = 0; v < inputArVolume; v++) {

        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<d*slen+inputArWidth; i++)
                outAr[v][i] = inputAr[v][i];
            for(int i=inputArWidth+d*slen; i<(d+1)*slen; i++)
                outAr[v][i] = inputAr[v][i]+inputAr[v][i-inputArWidth];
        }
    }

    return outAr;
}

//3D Forward sum on z axis.
double* sumForw3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len-inputArWidth*inputArHeight; i++) {
        outAr[i] = inputAr[i+inputArWidth*inputArHeight]+inputAr[i];
    }

    for(int i=len-inputArWidth*inputArHeight; i<len; i++) {
        outAr[i] = inputAr[i];
    }

    return outAr;
}

//3D Backward derivative on z axis.
double* sumBack3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=inputArWidth*inputArHeight; i<len; i++) {
        outAr[i] = inputAr[i]+inputAr[i-inputArWidth*inputArHeight];
    }

    for(int i=0; i<inputArWidth*inputArHeight; i++) {
        outAr[i] = inputAr[i];
    }

    return outAr;
}

//4D Forward sum on z axis.
double* *sumForw4DZ(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }
    for (int v = 0; v < inputArVolume; v++) {

        for(int i=0; i<len-inputArWidth*inputArHeight; i++) {
            outAr[v][i] = inputAr[v][i+inputArWidth*inputArHeight]+inputAr[v][i];
        }

        for(int i=len-inputArWidth*inputArHeight; i<len; i++) {
            outAr[v][i] = inputAr[v][i];
        }
    }

    return outAr;
}

//4D Backward derivative on z axis.
double* *sumBack4DZ(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }
    for (int v = 0; v < inputArVolume; v++) {

        for(int i=inputArWidth*inputArHeight; i<len; i++) {
            outAr[v][i] = inputAr[v][i]+inputAr[v][i-inputArWidth*inputArHeight];
        }

        for(int i=0; i<inputArWidth*inputArHeight; i++) {
            outAr[v][i] = inputAr[v][i];
        }
    }
    return outAr;
}

//4D Forward sum on t axis.
double* *sumForw4DT(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }
    int v;
    for (int i = 0; i < len; i++) {

        for(v=0; i<inputArVolume-1; i++) {
            outAr[v][i] = inputAr[v+1][i]+inputAr[v][i];
        }

        outAr[v][i] = inputAr[v][i];

    }

    return outAr;
}

//4D Backward derivative on t axis.
double* *sumBack4DT(double* *inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* *outAr = new double*[sizeof(double) * inputArVolume];
    for (int v = 0; v < inputArVolume; v++) {
        outAr[v] = new double[sizeof(double) * len];
    }
    int v;
    for (int i = 0; i < len; i++) {

        v = 0;
        outAr[v][i] = inputAr[v][i];

        for( v=1; v<inputArVolume; v++) {
            outAr[v][i] = inputAr[v][i]+inputAr[v-1][i];
        }
    }
    return outAr;
}


//Number of primes less than or equal to given intiger calculated by SIEVE OF ERATOSTHENES algorithm
int countPrimes(int n)
{
    if(n == 0 || n == 1) {
        int ans = 0;
        return ans;
    }

    bool *isprime = new bool[n + 1];
    for (int i = 2; i < n + 1; i ++) {
        isprime[i] = true;
    }
    for (int i = 2; i * i < n+1; i ++) {
        if (isprime[i]) {
            for (int j = i * i; j < n+1; j += i) {
                isprime[j] = false;
            }
        }
    }
    int cnt = 0;
    for (int i = 2; i < n+1; i ++) {
        if (isprime[i]) {
            cnt ++;
        }
    }
    return cnt;
}

// Function to generate list of prime numbers less than maxVal using Sieve of Eratosthenes
void SieveOfEratosthenes(vector<short int> &primes, int maxVal)
{
    // Create a boolean array "IsPrime[0..MAX_SIZE]" and
    // initialize all entries it as true. A value in
    // IsPrime[i] will finally be false if i is
    // Not a IsPrime, else true.
    bool IsPrime[maxVal];
    memset(IsPrime, true, sizeof(IsPrime));

    for (int p = 2; p * p < maxVal; p++)
    {
        // If IsPrime[p] is not changed, then it is a prime
        if (IsPrime[p] == true)
        {
            // Update all multiples of p greater than or
            // equal to the square of it
            // numbers which are multiple of p and are
            // less than p^2 are already been marked.
            for (int i = p * p; i <  maxVal; i += p)
                IsPrime[i] = false;
        }
    }

    // Store all prime numbers
    primes.push_back(0);
    for (int p = 2; p < maxVal; p++)
       if (IsPrime[p])
            primes.push_back(p);
}

// Utility function for Sieve of Sundaram
void sieveSundaram(vector<short int> &primes, int maxVal)
{
    // In general Sieve of Sundaram, produces primes smaller
    // than (2*x + 2) for a number given number x. Since
    // we want primes smaller than MAX, we reduce MAX to half
    // This array is used to separate numbers of the form
    // i + j + 2*i*j from others where 1 <= i <= j
    bool marked[maxVal/2 + 100] = {0};

    // Main logic of Sundaram. Mark all numbers which
    // do not generate prime number by doing 2*i+1
    for (int i=1; i<=(sqrt(maxVal)-1)/2; i++)
        for (int j=(i*(i+1))<<1; j<=maxVal/2; j=j+2*i+1)
            marked[j] = true;

    // Since 2 is a prime number
    primes.push_back(2);

    // Print other primes. Remaining primes are of the
    // form 2*i + 1 such that marked[i] is false.
    for (int i=1; i<=maxVal/2; i++)
        if (marked[i] == false)
            primes.push_back(2*i + 1);
}

//Binary algorithm for the GCD
short int gcd(short int a, short int b)
{
    while(b) b ^= a ^= b ^= a %= b;
        return a;
}

// Function to perform Goldbach's conjecture
short int* findPrimes(vector<short int> primes, int n)
{
    short int* prms =  new short int[sizeof(short int)*2];
    // Return if number is not even or less than 3
    if (n<=2 || n%2 != 0)
    {
        qDebug() << "Invalid Input \n";
        return 0;
    }

    // Check only upto half of number
    for (int i=0 ; primes[i] <= n/2; i++)
    {
        // find difference by subtracting current prime from n
        int diff = n - primes[i];

        // Search if the difference is also a prime number
        if (binary_search(primes.begin(), primes.end(), diff))
        {
            // Express as a sum of primes
//            qDebug() << primes[i] << " + " << diff << " = " << n;
            prms[0] = primes[i];
            prms[1] = diff;
            return prms;
        }
    }
}

//Cross product of two vectors
void crossProd(double vec1[3], double vec2[3], double prod[3]) {
//    double norm, square_norm;

    prod[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    prod[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    prod[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];

//    square_norm = prod[0]*prod[0] + prod[1]*prod[1] + prod[2]*prod[2];
//    norm = sqrt(square_norm);

//    prod[0] = prod[0]/norm;
//    prod[1] = prod[1]/norm;
//    prod[2] = prod[2]/norm;

    if(abs(vec1[0]*prod[0]+vec1[1]*prod[1]+vec1[2]*prod[2])>0.0000001 || abs(vec2[0]*prod[0]+vec2[1]*prod[1]+vec2[2]*prod[2])>0.0000001) {
        qDebug() << "From Inside: " << vec1[0]*prod[0]+vec1[1]*prod[1]+vec1[2]*prod[2] << vec2[0]*prod[0]+vec2[1]*prod[1]+vec2[2]*prod[2];
    }
}


// Function to convert 1D array to 2D slicewise array of 3D image
int** linear2twoDim(int* arr, int lenArr, int imgWid, int imgHei, int imgDep)
{
    int* sliceSizeCount = new int[sizeof(int)*imgDep];
    int* slcSzCnt = new int[sizeof(int)*imgDep];

    for (int i=0; i<imgDep; i++) {
        sliceSizeCount[i] = 0;
        slcSzCnt[i] = 0;
    }

    for (int i=0 ; i<lenArr; i++) {
        int i3 = floor(arr[i]/(imgWid*imgHei));
        sliceSizeCount[i3] = sliceSizeCount[i3] + 1;
    }

    // Memory allocation for 2D image slices array.
    int** imageSliceArr = new int*[sizeof(int)*imgDep];
    for(int i = 0; i < imgDep; ++i) {
        imageSliceArr[i] = new int[sizeof(int)*sliceSizeCount[i]];
    }

    for (int i=0 ; i<lenArr; i++) {
        int i3 = floor(arr[i]/(imgWid*imgHei));
        imageSliceArr[i3][slcSzCnt[i3]] = arr[i] - i3*imgWid*imgHei;
        slcSzCnt[i3] = slcSzCnt[i3] + 1;
    }

//    for(int i=0; i<imgDep; i++) {
//        for(int j=0; j<sliceSizeCount[i]; j++) {
//            qDebug() << i << j << imageSliceArr[i][j];
//        }
//    }

    return imageSliceArr;
}


QQueue<int> array2queue(int* arr, int len) {
    QQueue<int> queue;

    for(int i=0; i<len; i++) {
        queue.enqueue(arr[i]);
    }

    return queue;
}

double* nearNeighInit(double* arr, int width, int height, int depth, int* maskInds) {
    int len = width*height*depth;
//    QQueue<int> indQueue;
    //Memory allocation for temporary image arrays.
    double* initMask = new double[sizeof(double) * len];

    int maskSize = 0;
    for(int i=0; i<len; i++) {
        if(i==maskInds[maskSize]) {
            maskSize++;
        }
    }
    qDebug() << maskSize;

    int ind = 0;
    for(int i=0; i<len; i++) {
        if(i == maskInds[ind]) {
            ind++;
            continue;
        }

        int iZ = floor(i/(width*height));
        int iY = i - iZ*width*height;
        iY = floor(iY/width);
        int iX = i - iY*width - iZ*width*height;

        double minDis = 100;
        int minLoc = 0;
        for(int j=0; j<maskSize; j++) {
            int temp = maskInds[j];
            int tempZ = floor(temp/(width*height));
            int tempY = temp - tempZ*width*height;
            tempY = floor(tempY/width);
            int tempX = temp - tempY*width - tempZ*width*height;

            double dis = sqrt((tempX-iX)*(tempX-iX) + (tempY-iY)*(tempY-iY) + (tempZ-iZ)*(tempZ-iZ));
            if(dis < minDis) {
                minDis = dis;
                minLoc = temp;
            }
        }
        initMask[i] = arr[minLoc];

//        qDebug() << "i: " << i << "Nearest mask loc: " << minLoc << "Dis: " << minDis;
    }

//    indQueue = array2queue(maskInds, maskSize);
//    for(int i=0; i<len; i++) {
//        initMask[i] = -1;
//    }
//    for(int i=0; i<maskSize; i++) {
//        int indx = maskInds[i];
//        initMask[indx] = arr[indx];
//    }

//    qDebug() << "Hello1";

//    while (!indQueue.isEmpty()) {
//        int temp = indQueue.dequeue();
//        bool boundary = false;

//        int tempZ = floor(temp/(width*height));
//        int tempY = temp - tempZ*width*height;
//        tempY = floor(tempY/width);
//        int tempX = temp - tempY*width - tempZ*width*height;
// //        qDebug() << temp << tempX << tempY << tempZ;

//        if(tempX==0 || tempX==width || tempY==0 || tempY==height || tempZ==0 || tempZ==depth) {
//            boundary = true;
//        } else {
//            boundary = false;
//            initMask[tempX-1+tempY*width+tempZ*width*height] = arr[temp];
//            initMask[tempX+1+tempY*width+tempZ*width*height] = arr[temp];
//            initMask[tempX+(tempY-1)*width+tempZ*width*height] = arr[temp];
//            initMask[tempX+(tempY+1)*width+tempZ*width*height] = arr[temp];

//            initMask[tempX-1+tempY*width+(tempZ-1)*width*height] = arr[temp];
//            initMask[tempX+1+tempY*width+(tempZ-1)*width*height] = arr[temp];
//            initMask[tempX+(tempY-1)*width+(tempZ-1)*width*height] = arr[temp];
//            initMask[tempX+(tempY+1)*width+(tempZ-1)*width*height] = arr[temp];

//            initMask[tempX-1+tempY*width+(tempZ+1)*width*height] = arr[temp];
//            initMask[tempX+1+tempY*width+(tempZ+1)*width*height] = arr[temp];
//            initMask[tempX+(tempY-1)*width+(tempZ+1)*width*height] = arr[temp];
//            initMask[tempX+(tempY+1)*width+(tempZ+1)*width*height] = arr[temp];

//            indQueue.enqueue(tempX-1+tempY*width+tempZ*width*height);
//            indQueue.enqueue(tempX+1+tempY*width+tempZ*width*height);
//            indQueue.enqueue(tempX+(tempY-1)*width+tempZ*width*height);
//            indQueue.enqueue(tempX+(tempY+1)*width+tempZ*width*height);

//            indQueue.enqueue(tempX-1+tempY*width+(tempZ-1)*width*height);
//            indQueue.enqueue(tempX+1+tempY*width+(tempZ-1)*width*height);
//            indQueue.enqueue(tempX+(tempY-1)*width+(tempZ-1)*width*height);
//            indQueue.enqueue(tempX+(tempY+1)*width+(tempZ-1)*width*height);

//            indQueue.enqueue(tempX-1+tempY*width+(tempZ+1)*width*height);
//            indQueue.enqueue(tempX+1+tempY*width+(tempZ+1)*width*height);
//            indQueue.enqueue(tempX+(tempY-1)*width+(tempZ+1)*width*height);
//            indQueue.enqueue(tempX+(tempY+1)*width+(tempZ+1)*width*height);
//        }
//    }

//    int width1 = 7, height1 = 6, depth1 = 4;
//    int len = width1*height1*depth1;
//    QQueue<int> indQueue;
//    //Memory allocation for temporary image arrays.
//    double* initMask = new double[sizeof(double) * len];
//    double* arr1 = new double[sizeof(double) * len];
//    int* maskInds1 = new int[sizeof(int) * 29];

//    for(int i=0; i<len; i++) {
//        arr1[i] = i;
//    }

//    maskInds1[0] = 0;
//    maskInds1[1] = 2;
//    maskInds1[2] = 11;
//    maskInds1[3] = 19;
//    maskInds1[4] = 24;
//    maskInds1[5] = 32;
//    maskInds1[6] = 36;

//    maskInds1[7] = 43;
//    maskInds1[8] = 48;
//    maskInds1[9] = 51;
//    maskInds1[10] = 53;
//    maskInds1[11] = 60;
//    maskInds1[12] = 64;
//    maskInds1[13] = 69;
//    maskInds1[14] = 71;
//    maskInds1[15] = 80;
//    maskInds1[16] = 82;

//    maskInds1[17] = 84;
//    maskInds1[18] = 95;
//    maskInds1[19] = 97;
//    maskInds1[20] = 99;
//    maskInds1[21] = 103;
//    maskInds1[22] = 107;
//    maskInds1[23] = 108;
//    maskInds1[24] = 119;
//    maskInds1[25] = 120;
//    maskInds1[26] = 121;
//    maskInds1[27] = 125;
//    maskInds1[28] = 126;

//    qDebug() << "Hello";

//    int maskSize = 0;
//    for(int i=0; i<len; i++) {
//        if(i==maskInds1[maskSize]) {
//            maskSize++;
//        }
//    }
//    qDebug() << maskSize;

//    indQueue = array2queue(maskInds1, maskSize);
//    for(int i=0; i<len; i++) {
//        initMask[i] = -1;
//    }
//    for(int i=0; i<maskSize; i++) {
//        int indx = maskInds1[i];
//        initMask[indx] = arr1[indx];
//    }

//    qDebug() << "Hello1";

//    for(int i=0; i<len; i++) {
//        qDebug() << arr1[i];
//    }
//    for(int i=0; i<29; i++) {
//        qDebug() << "Ind Mask1: " << maskInds1[i];
//        int indx = maskInds1[i];
//        qDebug() << "Arr1: " << arr1[indx];
//    }

//    int ind = 0;
//    for(int i=0; i<len; i++) {
//        if(i == maskInds1[ind]) {
//            ind++;
//            continue;
//        }

//        int iZ = floor(i/(width1*height1));
//        int iY = i - iZ*width1*height1;
//        iY = floor(iY/width1);
//        int iX = i - iY*width1 - iZ*width1*height1;

//        double minDis = 100;
//        int minLoc = 0;
//        for(int j=0; j<29; j++) {
//            int temp = maskInds1[j];
//            int tempZ = floor(temp/(width1*height1));
//            int tempY = temp - tempZ*width1*height1;
//            tempY = floor(tempY/width1);
//            int tempX = temp - tempY*width1 - tempZ*width1*height1;

//            double dis = sqrt((tempX-iX)*(tempX-iX) + (tempY-iY)*(tempY-iY) + (tempZ-iZ)*(tempZ-iZ));
// //            qDebug() << dis;
//            if(dis < minDis) {
//                minDis = dis;
//                minLoc = temp;

//                qDebug() << "Disssssssss: " << dis << minDis;
//            }
//        }
//        initMask[i] = arr1[minLoc];

//        qDebug() << "i: " << i << "Nearest mask loc: " << minLoc << "Dis: " << minDis;
//    }


//    while (!indQueue.isEmpty()) {
//        int temp = indQueue.dequeue();
//        qDebug() << "Que: " << temp;
//    }

//    while (!indQueue.isEmpty()) {
//        int temp = indQueue.dequeue();
//        bool boundary = false;

//        qDebug() << "Que: " << temp;

//        int tempZ = floor(temp/(width1*height1));
//        int tempY = temp - tempZ*width1*height1;
//        tempY = floor(tempY/width1);
//        int tempX = temp - tempY*width1 - tempZ*width1*height1;
//         qDebug() << "Que: " << temp << tempX << tempY << tempZ;

//        if(tempX==0 || tempX==(width1-1) || tempY==0 || tempY==(height1-1) || tempZ==0 || tempZ==(depth1-1)) {
//            boundary = true;
//            qDebug() << temp << tempX << tempY << tempZ;
//        } else {
//            boundary = false;
//            initMask[tempX-1+tempY*width1+tempZ*width1*height1] = arr1[temp];
//            initMask[tempX+1+tempY*width1+tempZ*width1*height1] = arr1[temp];
//            initMask[tempX+(tempY-1)*width1+tempZ*width1*height1] = arr1[temp];
//            initMask[tempX+(tempY+1)*width1+tempZ*width1*height1] = arr1[temp];

//            initMask[tempX-1+tempY*width1+(tempZ-1)*width1*height1] = arr1[temp];
//            initMask[tempX+1+tempY*width1+(tempZ-1)*width1*height1] = arr1[temp];
//            initMask[tempX+(tempY-1)*width1+(tempZ-1)*width1*height1] = arr1[temp];
//            initMask[tempX+(tempY+1)*width1+(tempZ-1)*width1*height1] = arr1[temp];

//            initMask[tempX-1+tempY*width1+(tempZ+1)*width1*height1] = arr1[temp];
//            initMask[tempX+1+tempY*width1+(tempZ+1)*width1*height1] = arr1[temp];
//            initMask[tempX+(tempY-1)*width1+(tempZ+1)*width1*height1] = arr1[temp];
//            initMask[tempX+(tempY+1)*width1+(tempZ+1)*width1*height1] = arr1[temp];

//            indQueue.enqueue(tempX-1+tempY*width1+tempZ*width1*height1);
//            indQueue.enqueue(tempX+1+tempY*width1+tempZ*width1*height1);
//            indQueue.enqueue(tempX+(tempY-1)*width1+tempZ*width1*height1);
//            indQueue.enqueue(tempX+(tempY+1)*width1+tempZ*width1*height1);

//            indQueue.enqueue(tempX-1+tempY*width1+(tempZ-1)*width1*height1);
//            indQueue.enqueue(tempX+1+tempY*width1+(tempZ-1)*width1*height1);
//            indQueue.enqueue(tempX+(tempY-1)*width1+(tempZ-1)*width1*height1);
//            indQueue.enqueue(tempX+(tempY+1)*width1+(tempZ-1)*width1*height1);

//            indQueue.enqueue(tempX-1+tempY*width1+(tempZ+1)*width1*height1);
//            indQueue.enqueue(tempX+1+tempY*width1+(tempZ+1)*width1*height1);
//            indQueue.enqueue(tempX+(tempY-1)*width1+(tempZ+1)*width1*height1);
//            indQueue.enqueue(tempX+(tempY+1)*width1+(tempZ+1)*width1*height1);

//            break;

//            qDebug() << tempX-1 << tempY << tempZ;
//            qDebug() << tempX+1 << tempY << tempZ;
//            qDebug() << tempX << tempY-1 << tempZ;
//            qDebug() << tempX << tempY+1 << tempZ;
//        }
//    }

//    double sum = 0;
//    for(int i=0; i<len; i++) {
//        sum = sum + arr[i]-initMask[i];
//    }
//    qDebug() << sum;


//    int ind = 0, maskInd = 0;
//    if(ind == maskPxls[maskInd]) {
//        initMask[ind] = maskPxls[maskInd];
//        maskInd++;
//        continue;
//    }


//    ind = 0;
//    maskInd = 0;
//    while (!imgQueue.isEmpty()) {
//        double temp = imgQueue.dequeue();

//        if(ind == maskPxls[maskInd]) {
//            initMask[ind] = maskPxls[maskInd];
//            maskInd++;
//            continue;
//        }

//        //Find nearest mask pixel to the current pixel

//        initMask[ind] = arr[ind];
//        ind++;
//    }
    return initMask;
}

// Method to compare which one is the more close. It assumes that val2 is greater than val1 and target lies between these two.
int getClosest(int val1, int val2, int target) {
    if (target - val1 >= val2 - target)
        return val2;
    else
        return val1;
}

// Returns element closest to target in array
int findClosest(vector<short int> arr, int n, int target)
{
    // Corner cases
    if (target <= arr[0])
        return arr[0];
    if (target >= arr[n - 1])
        return arr[n - 1];

    // Doing binary search
    int i = 0, j = n, mid = 0;
    while (i < j) {
        mid = (i + j) / 2;

        if (arr[mid] == target)
            return arr[mid];

        /* If target is less than array element,
            then search in left */
        if (target < arr[mid]) {

            // If target is greater than previous
            // to mid, return closest of two
            if (mid > 0 && target > arr[mid - 1])
                return getClosest(arr[mid - 1],
                                  arr[mid], target);

            /* Repeat for left half */
            j = mid;
        }

        // If target is greater than mid
        else {
            if (mid < n - 1 && target < arr[mid + 1])
                return getClosest(arr[mid],
                                  arr[mid + 1], target);
            // update i
            i = mid + 1;
        }
    }

    // Only single element left after search
    return arr[mid];
}

//Convert all voxel values of an image to the closest prime numbers
void primize(double* imgArr, int imgWid, int imgHei, int imgDep, int maxVoxelVal) {
//    int numPrimes = countPrimes(maxVoxelVal);
//    qDebug() << numPrimes << numPrimes+1;

    vector<short int> primeList;
    SieveOfEratosthenes(primeList, maxVoxelVal+50);
//    for(int i=0; i<primeList.size(); i++) {
//        qDebug() << "Primes: " << i << primeList[i];
//    }

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        if(imgArr[i] == 0) {
            continue;
        }
        int tempVal = (int)imgArr[i];
        int nearPrime = findClosest(primeList, primeList.size(), tempVal);
        imgArr[i] = (double)nearPrime;
    }
}

//Convert all voxel values of an image to the number of primes less than or equal to voxel value
void primeCountImage(double* imgArr, int imgWid, int imgHei, int imgDep) {
    int numPrimes = 0;

    for(int i=0; i<imgWid*imgHei*imgDep; i++) {
        int temVal = (int)imgArr[i];
        numPrimes = countPrimes(temVal);
        imgArr[i] = (double)numPrimes;
    }
}

// Calculate Image "qntl" Quantile given in percentage.
double quantlCriterPM(double qntl, double* imgArr, int imgWidth, int imgHeight, int imgDepth, double sigma, int kernelSize) {
    int qntlIndex;
    double gausKernel[kernelSize];
    double *outConvX, *outConvXY, *outConvXYZ, *dervXConv, *dervYConv, *dervZConv;
    std::vector<double> convGradMag;

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    outConvX = convolution3DX(imgArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
    outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
    outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
    //First convolved derivatives
    dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                         // Derivative of convolved image w.r.t. x.
    dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                         // Derivative of convolved image w.r.t. y.
    dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                         // Derivative of convolved image w.r.t. z.

    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        //Calculate the convolved gradients.
        double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
        convGradMag.push_back(normi);
    }
    std::sort(convGradMag.begin(), convGradMag.end());

    qntlIndex = ceil(qntl*imgWidth*imgHeight*imgDepth/100);

//    qDebug() << convGradMag.size();

    delete [] outConvX;
    outConvX = NULL;
    delete [] outConvXY;
    outConvXY = NULL;
    delete [] outConvXYZ;
    outConvXYZ = NULL;
    delete [] dervXConv;
    dervXConv = NULL;
    delete [] dervYConv;
    dervYConv = NULL;
    delete [] dervZConv;
    dervZConv = NULL;

    return convGradMag[qntlIndex];
}

// Calculate Image "qntl" Quantile given in percentage.
double quantlCriterPM4Inpainting4D(double qntl, double* imgArr, double* binMask, int imgWidth, int imgHeight, int imgDepth, int imgVolume,  double sigma, int kernelSize) {
    int qntlIndex;
    double gausKernel[kernelSize], normi, norm_i_square;
    double *outConvX, *outConvXY, *outConvXYZ, *dervXConv, *dervYConv, *dervZConv;
    std::vector<double> convGradMag;

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    outConvX = convolution3DX(imgArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);           // Calculate the convolution on x axis.
    outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
    outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
    //First convolved derivatives
    dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                                // Derivative of convolved image w.r.t. x.
    dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                                // Derivative of convolved image w.r.t. y.
    dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                                // Derivative of convolved image w.r.t. z.

//    int randArrTraceIndex = 0;
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(binMask[i] == 0.0) {
            //Calculate the convolved gradients.
            norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i];
            normi = sqrt(norm_i_square);
            convGradMag.push_back(normi);
        }
    }
    std::sort(convGradMag.begin(), convGradMag.end());

    qntlIndex = ceil(qntl*convGradMag.size()/100);
//    qDebug() << convGradMag.size();

    delete [] outConvX;
    outConvX = NULL;
    delete [] outConvXY;
    outConvXY = NULL;
    delete [] outConvXYZ;
    outConvXYZ = NULL;
    delete [] dervXConv;
    dervXConv = NULL;
    delete [] dervYConv;
    dervYConv = NULL;
    delete [] dervZConv;
    dervZConv = NULL;

    return convGradMag[qntlIndex];
}

// Calculate Image "qntl" Quantile given in percentage.
double quantlCriterPM4Inpainting(double qntl, double* imgArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, double sigma, int kernelSize) {
    int qntlIndex;
    double gausKernel[kernelSize];
    double *outConvX, *outConvXY, *outConvXYZ, *dervXConv, *dervYConv, *dervZConv;
    std::vector<double> convGradMag;

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    outConvX = convolution3DX(imgArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
    outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
    outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
    //First convolved derivatives
    dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                         // Derivative of convolved image w.r.t. x.
    dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                         // Derivative of convolved image w.r.t. y.
    dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,1);                         // Derivative of convolved image w.r.t. z.

    int randArrTraceIndex = 0;
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(i == randPxls[randArrTraceIndex]) {
            randArrTraceIndex++;
            continue;
        }
        //Calculate the convolved gradients.
        double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
        convGradMag.push_back(normi);
    }
    std::sort(convGradMag.begin(), convGradMag.end());

//    //Write grad. magnitudes to file
//    QFile file("./img-outputs/masks/grad_magnitudes");
//    if(!file.open(QFile::WriteOnly | QFile::Text)) {
//        qDebug("Error", "File is NOT open");
//    }
//    QByteArray temp;
//    for(int i=0; i<convGradMag.size(); i++) {
//        char buf[9];
//        ::sprintf(buf, "%d", (int)round(convGradMag[i]));
//        temp.append(buf);
//        temp.append("\n");
//    }
//    file.write(temp);
//    file.flush();
//    file.close();

    qntlIndex = ceil(qntl*convGradMag.size()/100);

//    qDebug() << convGradMag.size();

    delete [] outConvX;
    outConvX = NULL;
    delete [] outConvXY;
    outConvXY = NULL;
    delete [] outConvXYZ;
    outConvXYZ = NULL;
    delete [] dervXConv;
    dervXConv = NULL;
    delete [] dervYConv;
    dervYConv = NULL;
    delete [] dervZConv;
    dervZConv = NULL;

    return convGradMag[qntlIndex];
}

// Calculate Image "qntl" Quantile given in percentage.
double histoCriterPM4Inpainting(double qntl, double* imgArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth) {
    int qntlIndex;
//    double *outConvX, *outConvXY, *outConvXYZ, *dervXConv, *dervYConv, *dervZConv;
    std::vector<double> imgHisto;

    int randArrTraceIndex = 0;
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(i == randPxls[randArrTraceIndex]) {
            randArrTraceIndex++;
            continue;
        }
        imgHisto.push_back(imgArr[i]);
    }
    std::sort(imgHisto.begin(), imgHisto.end());

    qntlIndex = ceil(qntl*imgHisto.size()/100);

//    qDebug() << imgHisto.size();

    return imgHisto[qntlIndex];
}

//Entropy Computation.
double entropy(double* arr, int len)
{
    double entropyValue = 0,temp=0, maxVal=-100000, minVal=100000;

    for(int i=0; i<len; i++) {
        if(arr[i] < minVal) {
            minVal = arr[i];
        } else if(arr[i] > maxVal) {
            maxVal = arr[i];
        }
    }

    qDebug() << minVal << maxVal;

    int histLen = (int)(maxVal-minVal+1);
    // Decleartion and memory allocation for each value
    double* histArr = new double[sizeof(double) * histLen];
    // initialize all intensity values to 0
    for(int i = 0; i<histLen; i++) {
        histArr[i] = 0;
    }

    // Histogram computation
    // Calculate the no of pixels for each intensity values
    double offset = 0;
    if(minVal < 0) {
        offset = abs(minVal);
    }
    for(int i = 0; i<len; i++) {
        arr[i] = offset+arr[i];
        histArr[(int)arr[i]]++;
    }

    for(int i=0; i<histLen; i++) { //the number of times a sybmol has occured
        if(histArr[i] > 0) { //log of zero goes to infinity
            temp = (histArr[i]/len)*(log2(histArr[i]/len));
            entropyValue += temp;
        }
    }

//    qDebug() << "Inside entropy: " << entropyValue*(-1);

    delete [] histArr;
    histArr = NULL;

    entropyValue = entropyValue*(-1);

    return entropyValue;
}

//Sum of squared differences (SSD) array of 4D(3D + q space) data in the x direction.
double* ssd_4DX_central(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=1; i<(slen-1); i++) {
            outAr[i+dlen] = 0;
            for(int t=0; t<timeLen; t++) {
//                qDebug() << "dlen:" << d << "slen" << i << "timeLen: " << t;
                outAr[i+dlen] = outAr[i+dlen] + ((inputAr[t][i+dlen]-inputAr[t][i+1+dlen])*(inputAr[t][i+dlen]-inputAr[t][i+1+dlen]) + (inputAr[t][i+dlen]-inputAr[t][i-1+dlen])*(inputAr[t][i+dlen]-inputAr[t][i-1+dlen]))/2;
            }
            outAr[i+dlen] = outAr[i+dlen]/timeLen;
        }
        for(int k=0; k<inputArHeight; k++) {
            outAr[k*inputArWidth+dlen] = 0;
            outAr[(k+1)*inputArWidth-1+dlen] = 0;
            for(int t=0; t<timeLen; t++) {
                outAr[k*inputArWidth+dlen] = outAr[k*inputArWidth+dlen] + (inputAr[t][k*inputArWidth+dlen]-inputAr[t][k*inputArWidth+1+dlen])*(inputAr[t][k*inputArWidth+dlen]-inputAr[t][k*inputArWidth+1+dlen]);
                outAr[(k+1)*inputArWidth-1+dlen] = outAr[(k+1)*inputArWidth-1+dlen] + (inputAr[t][(k+1)*inputArWidth-1+dlen]-inputAr[t][(k+1)*inputArWidth-2+dlen])*(inputAr[t][(k+1)*inputArWidth-1+dlen]-inputAr[t][(k+1)*inputArWidth-2+dlen]);
            }
            outAr[k*inputArWidth+dlen] = outAr[k*inputArWidth+dlen]/timeLen;
            outAr[(k+1)*inputArWidth-1+dlen] = outAr[(k+1)*inputArWidth-1+dlen]/timeLen;
        }
    }
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128,128) for derivatives
}

//Sum of squared differences (SSD) array of 4D(3D + q space) data in the y direction.
double* ssd_4DY_central(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=inputArWidth; i<slen-inputArWidth; i++) {
            outAr[i+dlen] = 0;
            for(int t=0; t<timeLen; t++) {
                outAr[i+dlen] = outAr[i+dlen] + ((inputAr[t][i+dlen]-inputAr[t][i+inputArWidth+dlen])*(inputAr[t][i+dlen]-inputAr[t][i+inputArWidth+dlen]) + (inputAr[t][i+dlen]-inputAr[t][i-inputArWidth+dlen])*(inputAr[t][i+dlen]-inputAr[t][i-inputArWidth+dlen]))/2;
            }
            outAr[i+dlen] = outAr[i+dlen]/timeLen;
        }
        for(int k=0; k<inputArWidth; k++) {
            outAr[k+dlen] = 0;
            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = 0;
            for(int t=0; t<timeLen; t++) {
                outAr[k+dlen] = outAr[k+dlen] + (inputAr[t][k+dlen]-inputAr[t][k+inputArWidth+dlen])*(inputAr[t][k+dlen]-inputAr[t][k+inputArWidth+dlen]);
                outAr[(inputArHeight-1)*inputArWidth+k+dlen] = outAr[(inputArHeight-1)*inputArWidth+k+dlen] + (inputAr[t][(inputArHeight-1)*inputArWidth+k+dlen]-inputAr[t][(inputArHeight-2)*inputArWidth+k+dlen])*(inputAr[t][(inputArHeight-1)*inputArWidth+k+dlen]-inputAr[t][(inputArHeight-2)*inputArWidth+k+dlen]);
            }
            outAr[k+dlen] = outAr[k+dlen]/timeLen;
            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = outAr[(inputArHeight-1)*inputArWidth+k+dlen]/timeLen;
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//Sum of squared differences (SSD) array of 4D(3D + q space) data in the z direction.
double* ssd_4DZ_central(double** inputAr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    int slen = inputArWidth*inputArHeight;

    for(int i=slen; i<len-slen; i++) {
        outAr[i] = 0;
        for(int t=0; t<timeLen; t++) {
            outAr[i] = outAr[i] + ((inputAr[t][i]-inputAr[t][i-slen])*(inputAr[t][i]-inputAr[t][i-slen])+(inputAr[t][i]-inputAr[t][i+slen])*(inputAr[t][i]-inputAr[t][i+slen]))/2;
        }
        outAr[i] = outAr[i]/timeLen;
    }

    for(int h=0; h<inputArHeight; h++) {
        int hlen = h*inputArWidth;
        for(int k=0; k<inputArWidth; k++) {
            outAr[k+hlen] = 0;
            outAr[len-slen+k+hlen] = 0;
            for(int t=0; t<timeLen; t++) {
                outAr[k+hlen] = outAr[k+hlen] + (inputAr[t][k+slen+hlen]-inputAr[t][k+hlen])*(inputAr[t][k+slen+hlen]-inputAr[t][k+hlen]);
                outAr[len-slen+k+hlen] = outAr[len-slen+k+hlen] + (inputAr[t][len-slen+k+hlen]-inputAr[t][len-2*slen+k+hlen])*(inputAr[t][len-slen+k+hlen]-inputAr[t][len-2*slen+k+hlen]);
            }
            outAr[k+hlen] = outAr[k+hlen]/timeLen;
            outAr[len-slen+k+hlen] = outAr[len-slen+k+hlen]/timeLen;
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//For non-dropouts: Sum of squared differences (SSD) array of 4D(3D + q space) data in the x direction.
double* ssd_4DX_central_nonDropouts(double** inputAr, double** maskArr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=1; i<(slen-1); i++) {
            outAr[i+dlen] = 0;
            int timeCnt = 0;
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][i+dlen] == 1.0 || maskArr[t][i+1+dlen] == 1.0 || maskArr[t][i-1+dlen] == 1.0) {
                    continue;
                } else {
                    outAr[i+dlen] = outAr[i+dlen] + ((inputAr[t][i+dlen]-inputAr[t][i+1+dlen])*(inputAr[t][i+dlen]-inputAr[t][i+1+dlen]) + (inputAr[t][i+dlen]-inputAr[t][i-1+dlen])*(inputAr[t][i+dlen]-inputAr[t][i-1+dlen]))/2;
                    timeCnt++;
                }
            }
            outAr[i+dlen] = outAr[i+dlen]/timeCnt;
        }
        for(int k=0; k<inputArHeight; k++) {
            outAr[k*inputArWidth+dlen] = 0;
            outAr[(k+1)*inputArWidth-1+dlen] = 0;
            int timeCnt1 = 0, timeCnt2 = 0;
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][k*inputArWidth+dlen] == 1.0 || maskArr[t][k*inputArWidth+1+dlen] == 1.0) {
                    continue;
                } else {
                    outAr[k*inputArWidth+dlen] = outAr[k*inputArWidth+dlen] + (inputAr[t][k*inputArWidth+dlen]-inputAr[t][k*inputArWidth+1+dlen])*(inputAr[t][k*inputArWidth+dlen]-inputAr[t][k*inputArWidth+1+dlen]);
                    timeCnt1++;
                }
            }
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][(k+1)*inputArWidth-1+dlen] == 1.0 || maskArr[t][(k+1)*inputArWidth-2+dlen] == 1.0) {
                    continue;
                } else {
                    outAr[(k+1)*inputArWidth-1+dlen] = outAr[(k+1)*inputArWidth-1+dlen] + (inputAr[t][(k+1)*inputArWidth-1+dlen]-inputAr[t][(k+1)*inputArWidth-2+dlen])*(inputAr[t][(k+1)*inputArWidth-1+dlen]-inputAr[t][(k+1)*inputArWidth-2+dlen]);
                    timeCnt2++;
                }
            }
            outAr[k*inputArWidth+dlen] = outAr[k*inputArWidth+dlen]/timeCnt1;
            outAr[(k+1)*inputArWidth-1+dlen] = outAr[(k+1)*inputArWidth-1+dlen]/timeCnt2;
        }
    }
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128,128) for derivatives
}

//For non-dropouts: Sum of squared differences (SSD) array of 4D(3D + q space) data in the y direction.
double* ssd_4DY_central_nonDropouts(double** inputAr, double** maskArr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=inputArWidth; i<slen-inputArWidth; i++) {
            outAr[i+dlen] = 0;
            int timeCnt = 0;
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][i+dlen] == 1.0 || maskArr[t][i+inputArWidth+dlen] == 1.0 || maskArr[t][i-inputArWidth+dlen] == 1.0) {
                    continue;
                } else {
                    outAr[i+dlen] = outAr[i+dlen] + ((inputAr[t][i+dlen]-inputAr[t][i+inputArWidth+dlen])*(inputAr[t][i+dlen]-inputAr[t][i+inputArWidth+dlen]) + (inputAr[t][i+dlen]-inputAr[t][i-inputArWidth+dlen])*(inputAr[t][i+dlen]-inputAr[t][i-inputArWidth+dlen]))/2;
                    timeCnt++;
                }
            }
            outAr[i+dlen] = outAr[i+dlen]/timeCnt;
        }
        for(int k=0; k<inputArWidth; k++) {
            outAr[k+dlen] = 0;
            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = 0;
            int timeCnt1 = 0, timeCnt2 = 0;
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][k+dlen] == 1.0 || maskArr[t][k+inputArWidth+dlen] == 1.0) {
                    continue;
                } else {
                    outAr[k+dlen] = outAr[k+dlen] + (inputAr[t][k+dlen]-inputAr[t][k+inputArWidth+dlen])*(inputAr[t][k+dlen]-inputAr[t][k+inputArWidth+dlen]);
                    timeCnt1++;
                }
            }
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][(inputArHeight-1)*inputArWidth+k+dlen] == 1.0 || maskArr[t][(inputArHeight-2)*inputArWidth+k+dlen] == 1.0) {
                    continue;
                } else {
                    outAr[(inputArHeight-1)*inputArWidth+k+dlen] = outAr[(inputArHeight-1)*inputArWidth+k+dlen] + (inputAr[t][(inputArHeight-1)*inputArWidth+k+dlen]-inputAr[t][(inputArHeight-2)*inputArWidth+k+dlen])*(inputAr[t][(inputArHeight-1)*inputArWidth+k+dlen]-inputAr[t][(inputArHeight-2)*inputArWidth+k+dlen]);
                    timeCnt2++;
                }
            }
            outAr[k+dlen] = outAr[k+dlen]/timeCnt1;
            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = outAr[(inputArHeight-1)*inputArWidth+k+dlen]/timeCnt2;
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//For non-dropouts: Sum of squared differences (SSD) array of 4D(3D + q space) data in the Z direction.
double* ssd_4DZ_central_nonDropouts(double** inputAr, double** maskArr, int timeLen, int inputArWidth, int inputArHeight, int inputArDepth)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    int slen = inputArWidth*inputArHeight;

    for(int i=slen; i<len-slen; i++) {
        outAr[i] = 0;
        int timeCnt = 0;
        for(int t=0; t<timeLen; t++) {
            if(maskArr[t][i] == 1.0 || maskArr[t][i-slen] == 1.0 || maskArr[t][i+slen] == 1.0) {
                continue;
            } else {
                outAr[i] = outAr[i] + ((inputAr[t][i]-inputAr[t][i-slen])*(inputAr[t][i]-inputAr[t][i-slen])+(inputAr[t][i]-inputAr[t][i+slen])*(inputAr[t][i]-inputAr[t][i+slen]))/2;
                timeCnt++;
            }
        }
        outAr[i] = outAr[i]/timeCnt;
    }

    for(int h=0; h<inputArHeight; h++) {
        int hlen = h*inputArWidth;
        for(int k=0; k<inputArWidth; k++) {
            outAr[k+hlen] = 0;
            outAr[len-slen+k+hlen] = 0;
            int timeCnt1 = 0, timeCnt2 = 0;
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][k+hlen] == 1.0 || maskArr[t][k+slen+hlen] == 1.0) {
                    continue;
                } else {
                    outAr[k+hlen] = outAr[k+hlen] + (inputAr[t][k+slen+hlen]-inputAr[t][k+hlen])*(inputAr[t][k+slen+hlen]-inputAr[t][k+hlen]);
                    timeCnt1++;
                }
            }
            for(int t=0; t<timeLen; t++) {
                if(maskArr[t][len-slen+k+hlen] == 1.0 || maskArr[t][len-2*slen+k+hlen] == 1.0) {
                    continue;
                } else {
                    outAr[len-slen+k+hlen] = outAr[len-slen+k+hlen] + (inputAr[t][len-slen+k+hlen]-inputAr[t][len-2*slen+k+hlen])*(inputAr[t][len-slen+k+hlen]-inputAr[t][len-2*slen+k+hlen]);
                    timeCnt2++;
                }
            }
            outAr[k+hlen] = outAr[k+hlen]/timeCnt1;
            outAr[len-slen+k+hlen] = outAr[len-slen+k+hlen]/timeCnt2;
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//~ Create an array of random n different voxel locations (from the same volume) from the mask image for 4D data.
double* randVoxLocs(int locLen, double** maskImg, int volNum, int upperBound)
{
    // Decleartion and memory allocation
    double* randVoxelsArr = new double[sizeof(double) * locLen];
//    srand(time(NULL));    //For generating different random scattered images
    srand(1);

    bool duplicate, nonMask; //variable to check or number is already used
    //for loop to generate a complete set of 10 random numbers
    for (int i = 0; i<locLen; i++) {
        duplicate = false; // set check to false
        nonMask = false;  // set check to false
        // do while loop used to generate random numbers until a distinct random number is generated
        do {
//            randVoxelsArr[i] = rand()%(upperBound-lowerBound+1)+lowerBound; // generates a random number and stores it into randVoxelsArr[i]
            randVoxelsArr[i] = rand()%upperBound + 1; // generates a random number and stores it into randVoxelsArr[i]
            if(maskImg[volNum][int(randVoxelsArr[i])] == 0.0) {
                nonMask = true;
            } else {
                nonMask = false;
            }
            // for loop used to check the other numbers in set for any repeats
            for (int j = i-1; j>-1; j--) // works backwards from the recently generated element to element 0
                if (randVoxelsArr[i] == randVoxelsArr[j]) //checks if number is already used
                    duplicate = true; //sets duplicate to true to indicate there is a repeat
        } while (duplicate || nonMask); //loops until a new, distinct number is generated
    }
    return randVoxelsArr;
}

// Return corrupted volume numbers as a vector of integers
vector<int> corruptedVolumes(double** binMaskArr, int imgWidth, int imgHeight, int imgDepth, int numVols) {
    vector<int> volNumArr;

    for(int i=0; i<numVols; i++) {
//        double *foo = std::find(std::begin(binMaskArr[i]), std::end(binMaskArr[i]), 1);

        // When the element is not found, std::find returns the end of the range
//        if (foo != std::end(binMaskArr[i])) {
//            volNumArr.push_back(i);
//        }
        bool corVol = false;
        for(int j=0; j<imgWidth*imgDepth*imgHeight; j++) {
            if(binMaskArr[i][j] == 1) {
                corVol = true;
                volNumArr.push_back(i);
                break;
            }
        }
    }
    return volNumArr;
}

// Given a volume return corrupted slice numbers as a vector of integers
vector<int> corruptedSlices(double* binMaskArr, int imgWidth, int imgHeight, int imgDepth) {
    vector<int> slcNumArr;

    for(int i=0; i<imgDepth; i++) {
        bool corSlc = false;
        for(int j=imgWidth*imgHeight*i; j<imgWidth*imgHeight*(i+1); j++) {
            if(binMaskArr[j] == 1) {
                corSlc = true;
                slcNumArr.push_back(i);
                break;
            }
        }
    }
    return slcNumArr;
}

// Make a regular grid mask from a slice mask. (Currently, there are 4 modes, TL(0), TR(1), BL(2), BR(3))
void makeRegGrid_on_slice(double* binMaskrr, int imgWidth, int imgHeight, int currSlice, int gridLen, int mode) {
//    // Unmasking
//    for(int i=currSlice*imgHeight*imgWidth; i<(currSlice+1)*imgHeight*imgWidth; i++) {
//        binMaskrr[i] = 0;
//    }

    if(mode == 0) {
        for(int i=0; i<imgHeight; i=i+gridLen) {
            for(int j=0; j<imgWidth; j=j+gridLen) {
                binMaskrr[j + i*imgWidth + currSlice*imgWidth*imgHeight] = 0;
            }
        }
    } else if(mode == 1) {
        for(int i=0; i<imgHeight; i=i+gridLen) {
            for(int j=1; j<imgWidth; j=j+gridLen) {
                binMaskrr[j + i*imgWidth + currSlice*imgWidth*imgHeight] = 0;
            }
        }
    } else if(mode == 2) {
        for(int i=1; i<imgHeight; i=i+gridLen) {
            for(int j=0; j<imgWidth; j=j+gridLen) {
                binMaskrr[j + i*imgWidth + currSlice*imgWidth*imgHeight] = 0;
            }
        }
    } else if(mode == 3) {
        for(int i=1; i<imgHeight; i=i+gridLen) {
            for(int j=1; j<imgWidth; j=j+gridLen) {
                binMaskrr[j + i*imgWidth + currSlice*imgWidth*imgHeight] = 0;
            }
        }
    } else {
        qDebug() << "Invalid mode!";
    }
}


//3D self-created test data
// For Test: Decleration and Memory allocation
//double* testImg = new double[sizeof(double) * 4*4*4];
//double* testRecImg = new double[sizeof(double) * 4*4*4];
//testImg[0] = 10;
//testImg[1] = 8;
//testImg[2] = 6;
//testImg[3] = 3;
//testImg[4] = 7;
//testImg[5] = 9;
//testImg[6] = 6;
//testImg[7] = 7;
//testImg[8] = 1;
//testImg[9] = 0;
//testImg[10] = 0;
//testImg[11] = 0;
//testImg[12] = 0;
//testImg[13] = 0;
//testImg[14] = 2;
//testImg[15] = 1;
//testImg[16] = 2;
//testImg[17] = 3;
//testImg[18] = 4;
//testImg[19] = 5;
//testImg[20] = 2;
//testImg[21] = 3;
//testImg[22] = 4;
//testImg[23] = 5;
//testImg[24] = 5;
//testImg[25] = 1;
//testImg[26] = 2;
//testImg[27] = 1;
//testImg[28] = 2;
//testImg[29] = 6;
//testImg[30] = 7;
//testImg[31] = 8;
//testImg[32] = 9;
//testImg[33] = 10;
//testImg[34] = 11;
//testImg[35] = 14;
//testImg[36] = 18;
//testImg[37] = 16;
//testImg[38] = 13;
//testImg[39] = 17;
//testImg[40] = 19;
//testImg[41] = 16;
//testImg[42] = 17;
//testImg[43] = 11;
//testImg[44] = 10;
//testImg[45] = 10;
//testImg[46] = 10;
//testImg[47] = 10;
//testImg[48] = 20;
//testImg[49] = 22;
//testImg[50] = 21;
//testImg[51] = 12;
//testImg[52] = 13;
//testImg[53] = 14;
//testImg[54] = 15;
//testImg[55] = 12;
//testImg[56] = 13;
//testImg[57] = 14;
//testImg[58] = 15;
//testImg[59] = 15;
//testImg[60] = 11;
//testImg[61] = 12;
//testImg[62] = 11;
//testImg[63] = 12;
//testImg[64] = 10;
//testImg[65] = 8;
//testImg[66] = 6;
//testImg[67] = 3;
//testImg[68] = 7;
//testImg[69] = 9;
//testImg[70] = 6;
//testImg[71] = 7;
//testImg[72] = 1;
//testImg[73] = 0;
//testImg[74] = 0;
//testImg[75] = 0;
//testImg[76] = 0;
//testImg[77] = 0;
//testImg[78] = 2;
//testImg[79] = 1;
//testImg[80] = 2;
//testImg[81] = 3;
//testImg[82] = 4;
//testImg[83] = 5;
//testImg[84] = 2;
//testImg[85] = 3;
//testImg[86] = 4;
//testImg[87] = 5;
//testImg[88] = 5;
//testImg[89] = 1;
//testImg[90] = 2;
//testImg[91] = 1;
//testImg[92] = 2;
//testImg[93] = 6;
//testImg[94] = 7;
//testImg[95] = 8;
//testImg[96] = 9;
//testImg[97] = 10;
//testImg[98] = 11;
//testImg[99] = 14;
//testImg[100] = 18;
//testImg[101] = 16;
//testImg[102] = 13;
//testImg[103] = 17;
//testImg[104] = 19;
//testImg[105] = 16;
//testImg[106] = 17;
//testImg[107] = 11;
//testImg[108] = 10;
//testImg[109] = 10;
//testImg[110] = 10;
//testImg[111] = 10;
//testImg[112] = 20;
//testImg[113] = 22;
//testImg[114] = 21;
//testImg[115] = 12;
//testImg[116] = 13;
//testImg[117] = 14;
//testImg[118] = 15;
//testImg[119] = 12;
//testImg[120] = 13;
//testImg[121] = 14;
//testImg[122] = 15;
//testImg[123] = 15;
//testImg[124] = 11;
//testImg[125] = 12;
//testImg[126] = 11;
//testImg[127] = 12;
//testImg[128] = 18;
//testImg[129] = 16;
//testImg[130] = 13;
//testImg[131] = 17;
//testImg[132] = 19;
//testImg[133] = 16;
//testImg[134] = 17;
//testImg[135] = 11;
//testImg[136] = 10;
//testImg[137] = 10;
//testImg[138] = 10;
//testImg[139] = 10;
//testImg[140] = 20;
//testImg[141] = 22;
//testImg[142] = 21;
//testImg[143] = 12;
//testImg[144] = 13;
//testImg[145] = 14;
//testImg[146] = 15;
//testImg[147] = 12;
//testImg[148] = 13;
//testImg[149] = 14;
//testImg[150] = 15;
//testImg[151] = 15;
//testImg[152] = 11;
//testImg[153] = 12;
//testImg[154] = 11;
//testImg[155] = 12;
//testImg[156] = 14;
//testImg[157] = 15;
//testImg[158] = 15;
//testImg[159] = 11;
//testImg[160] = 12;
//testImg[161] = 11;
//testImg[162] = 12;
//testImg[163] = 18;
//testImg[164] = 16;
//testImg[165] = 13;
//testImg[166] = 17;
//testImg[167] = 19;
//testImg[168] = 16;
//testImg[169] = 17;
//testImg[170] = 11;
//testImg[171] = 10;
//testImg[172] = 10;
//testImg[173] = 10;
//testImg[174] = 10;

//testRecImg[0] = 10;
//testRecImg[1] = 4;
//testRecImg[2] = 6;
//testRecImg[3] = 3;
//testRecImg[4] = 6;
//testRecImg[5] = 9;
//testRecImg[6] = 6;
//testRecImg[7] = 9;
//testRecImg[8] = 2;
//testRecImg[9] = 0;
//testRecImg[10] = 1;
//testRecImg[11] = 2;
//testRecImg[12] = 0;
//testRecImg[13] = 1;
//testRecImg[14] = 1;
//testRecImg[15] = 1;
//testRecImg[16] = 3;
//testRecImg[17] = 4;
//testRecImg[18] = 4;
//testRecImg[19] = 6;
//testRecImg[20] = 1;
//testRecImg[21] = 3;
//testRecImg[22] = 0;
//testRecImg[23] = 0;
//testRecImg[24] = 5;
//testRecImg[25] = 0;
//testRecImg[26] = 0;
//testRecImg[27] = 1;
//testRecImg[28] = 1;
//testRecImg[29] = 8;
//testRecImg[30] = 7;
//testRecImg[31] = 0;
//testRecImg[32] = 10;
//testRecImg[33] = 10;
//testRecImg[34] = 7;
//testRecImg[35] = 7;
//testRecImg[36] = 18;
//testRecImg[37] = 9;
//testRecImg[38] = 17;
//testRecImg[39] = 17;
//testRecImg[40] = 17;
//testRecImg[41] = 16;
//testRecImg[42] = 17;
//testRecImg[43] = 12;
//testRecImg[44] = 11;
//testRecImg[45] = 10;
//testRecImg[46] = 12;
//testRecImg[47] = 13;
//testRecImg[48] = 20;
//testRecImg[49] = 23;
//testRecImg[50] = 22;
//testRecImg[51] = 12;
//testRecImg[52] = 10;
//testRecImg[53] = 13;
//testRecImg[54] = 15;
//testRecImg[55] = 15;
//testRecImg[56] = 13;
//testRecImg[57] = 14;
//testRecImg[58] = 15;
//testRecImg[59] = 15;
//testRecImg[60] = 11;
//testRecImg[61] = 11;
//testRecImg[62] = 10;
//testRecImg[63] = 12;
//testRecImg[64] = 10;
//testRecImg[65] = 8;
//testRecImg[66] = 6;
//testRecImg[67] = 3;
//testRecImg[68] = 7;
//testRecImg[69] = 9;
//testRecImg[70] = 6;
//testRecImg[71] = 7;
//testRecImg[72] = 1;
//testRecImg[73] = 0;
//testRecImg[74] = 0;
//testRecImg[75] = 0;
//testRecImg[76] = 0;
//testRecImg[77] = 0;
//testRecImg[78] = 2;
//testRecImg[79] = 1;
//testRecImg[80] = 2;
//testRecImg[81] = 3;
//testRecImg[82] = 4;
//testRecImg[83] = 5;
//testRecImg[84] = 2;
//testRecImg[85] = 3;
//testRecImg[86] = 4;
//testRecImg[87] = 5;
//testRecImg[88] = 5;
//testRecImg[89] = 1;
//testRecImg[90] = 2;
//testRecImg[91] = 1;
//testRecImg[92] = 2;
//testRecImg[93] = 6;
//testRecImg[94] = 7;
//testRecImg[95] = 8;
//testRecImg[96] = 9;
//testRecImg[97] = 10;
//testRecImg[98] = 11;
//testRecImg[99] = 14;
//testRecImg[100] = 18;
//testRecImg[101] = 16;
//testRecImg[102] = 13;
//testRecImg[103] = 17;
//testRecImg[104] = 19;
//testRecImg[105] = 16;
//testRecImg[106] = 17;
//testRecImg[107] = 11;
//testRecImg[108] = 10;
//testRecImg[109] = 10;
//testRecImg[110] = 10;
//testRecImg[111] = 10;
//testRecImg[112] = 20;
//testRecImg[113] = 22;
//testRecImg[114] = 21;
//testRecImg[115] = 12;
//testRecImg[116] = 13;
//testRecImg[117] = 14;
//testRecImg[118] = 15;
//testRecImg[119] = 12;
//testRecImg[120] = 13;
//testRecImg[121] = 14;
//testRecImg[122] = 15;
//testRecImg[123] = 15;
//testRecImg[124] = 11;
//testRecImg[125] = 12;
//testRecImg[126] = 11;
//testRecImg[127] = 12;
//testRecImg[128] = 18;
//testRecImg[129] = 16;
//testRecImg[130] = 13;
//testRecImg[131] = 17;
//testRecImg[132] = 19;
//testRecImg[133] = 16;
//testRecImg[134] = 17;
//testRecImg[135] = 11;
//testRecImg[136] = 10;
//testRecImg[137] = 10;
//testRecImg[138] = 10;
//testRecImg[139] = 10;
//testRecImg[140] = 20;
//testRecImg[141] = 22;
//testRecImg[142] = 21;
//testRecImg[143] = 12;
//testRecImg[144] = 13;
//testRecImg[145] = 14;
//testRecImg[146] = 15;
//testRecImg[147] = 12;
//testRecImg[148] = 13;
//testRecImg[149] = 14;
//testRecImg[150] = 15;
//testRecImg[151] = 15;
//testRecImg[152] = 11;
//testRecImg[153] = 12;
//testRecImg[154] = 11;
//testRecImg[155] = 12;
//testRecImg[156] = 14;
//testRecImg[157] = 15;
//testRecImg[158] = 15;
//testRecImg[159] = 11;
//testRecImg[160] = 12;
//testRecImg[161] = 11;
//testRecImg[162] = 12;
//testRecImg[163] = 18;
//testRecImg[164] = 16;
//testRecImg[165] = 13;
//testRecImg[166] = 17;
//testRecImg[167] = 19;
//testRecImg[168] = 16;
//testRecImg[169] = 17;
//testRecImg[170] = 11;
//testRecImg[171] = 10;
//testRecImg[172] = 10;
//testRecImg[173] = 10;
//testRecImg[174] = 10;
//int* rpxl1 = new int[sizeof(int) * 24];
//rpxl1[0] = 0;
//rpxl1[1] = 3;
//rpxl1[2] = 6;
//rpxl1[3] = 9;
//rpxl1[4] = 12;
//rpxl1[5] = 15;
//rpxl1[6] = 18;
//rpxl1[7] = 21;
//rpxl1[8] = 24;
//rpxl1[9] = 27;
//rpxl1[10] = 30;
//rpxl1[11] = 33;
//rpxl1[12] = 36;
//rpxl1[13] = 39;
//rpxl1[14] = 42;
//rpxl1[15] = 45;
//rpxl1[16] = 48;
//rpxl1[17] = 51;
//rpxl1[18] = 54;
//rpxl1[19] = 57;
//rpxl1[20] = 60;
//rpxl1[21] = 63;
//rpxl1[22] = 115;
//rpxl1[23] = 168;
//int dilMsk = 0;
//double* rpxl2 = dilatedMask3D(7, 5, 5, rpxl1, dilMsk, false);
//qDebug() << dilMsk;
//for(int i=0; i<140; i++) {
//    qDebug() << "Dilat: " << i << rpxl2[i] << testRecImg[(int)rpxl2[i]];
//}
