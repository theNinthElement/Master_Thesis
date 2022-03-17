#include "derivatives.h"
#include "supplementary_functions.h"
#include "fed.h"

#include <QDebug>

//***************************************************************************************************************************************************************************************************************************//
void eed_3d_smoothing(float timeStepSize, int iterNum, double contPar, double* imageArr, int imgWidth, int imgHeight, int imgDepth) {
    clock_t begin = clock();

    int kernelSize = 3;
    float stepSize = 1;
    double sigma = 1.0, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *OrigImgArray;

    //Memory allocation for diffusion tensor entries dij.
    double* d11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* d33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    //Memory allocation for temporary image array.
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    OrigImgArray = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    for(int i=0; i<imgWidth*imgHeight; i++) {
        OrigImgArray[i] = imageArr[i];
    }

    for(int iter=1; iter<iterNum; iter++) {
        outConvX = convolution3DX(imageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);         // Calculate the convolution on x axis.
        outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
        outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
        //First convolved derivatives
        dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
        dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
        dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.
        //First derivatives
        dervX = derivative3DX(imageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
        dervY = derivative3DY(imageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
        dervZ = derivative3DZ(imageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. z.

        dervForX = dervForw3DX(imageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervForY = dervForw3DY(imageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervForZ = dervForw3DZ(imageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervBacX = dervBack3DX(imageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervBacY = dervBack3DY(imageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervBacZ = dervBack3DZ(imageArr, imgWidth, imgHeight, imgDepth, stepSize);

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            //Define eigenvectors and entries for the diffusion tensor.
            double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3, v1[3], v2[3], v3[3];
            double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);

            if(norm_i_square == 0) {
                v1[0] = 1;
                v1[1] = 0;
                v1[2] = 0;

                v2[0] = 0;
                v2[1] = 1;
                v2[2] = 0;

                v3[0] = 0;
                v3[1] = 0;
                v3[2] = 1;
            } else if(xy_norm_i_square == 0) {
                v1[0] = 0;
                v1[1] = 0;
                v1[2] = dervZConv[i]/normi;

                v2[0] = 0;
                v2[1] = 0;
                v2[2] = 0;

                v3[0] = 0;
                v3[1] = 0;
                v3[2] = 0;
            } else {
                v1[0] = dervXConv[i]/normi;
                v1[1] = dervYConv[i]/normi;
                v1[2] = dervZConv[i]/normi;

                v2[0] = dervYConv[i]/xy_norm_i;
                v2[1] = -dervXConv[i]/xy_norm_i;
                v2[2] = 0;

                v3[0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
                v3[1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
                v3[2] = -xy_norm_i_square/(xy_norm_i*normi);
            }

            d1 = pm_diff(norm_i_square, contPar);
            d2 = 1;
            d3 = 1;

            //Diffusion tensor entries definition with eigenvalues d1, d2, d3 and orthgonal orthonormal eigenvectors v1, v2 and v3.
            d11[i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0] + d3*v3[0]*v3[0];
            d12[i] = (d1*v1[0]*v1[1] + d2*v2[0]*v2[1] + d3*v3[0]*v3[1])*dervY[i];
            d13[i] = (d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2])*dervZ[i];
            d21[i] = (d1*v1[1]*v1[0] + d2*v2[1]*v2[0] + d3*v3[1]*v3[0])*dervX[i];
            d22[i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1] + d3*v3[1]*v3[1];
            d23[i] = (d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2])*dervZ[i];
            d31[i] = (d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2])*dervX[i];
            d32[i] = (d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2])*dervY[i];
            d33[i] = d1*v1[2]*v1[2] + d2*v2[2]*v2[2] + d3*v3[2]*v3[2];
        }
        sumForX = sumForw3DX(d11,imgWidth,imgHeight,imgDepth);
        sumForY = sumForw3DY(d22,imgWidth,imgHeight,imgDepth);
        sumForZ = sumForw3DZ(d33,imgWidth,imgHeight,imgDepth);
        sumBacX = sumBack3DX(d11,imgWidth,imgHeight,imgDepth);
        sumBacY = sumBack3DY(d22,imgWidth,imgHeight,imgDepth);
        sumBacZ = sumBack3DZ(d33,imgWidth,imgHeight,imgDepth);

        dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,1);
        dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,1);
        dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,1);
        dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,1);
        dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,1);
        dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,1);

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            tempImgArrayPrev[i] = imageArr[i];

            imageArr[i] = imageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i]);
        }
        delete [] dervX;
        dervX = NULL;
        delete [] dervY;
        dervY = NULL;
        delete [] dervZ;
        dervZ = NULL;
        delete [] dervXD12;
        dervXD12 = NULL;
        delete [] dervXD13;
        dervXD13 = NULL;
        delete [] dervYD21;
        dervYD21 = NULL;
        delete [] dervYD23;
        dervYD23 = NULL;
        delete [] dervZD31;
        dervZD31 = NULL;
        delete [] dervZD32;
        dervZD32 = NULL;
        delete [] dervForX;
        dervForX = NULL;
        delete [] dervForY;
        dervForY = NULL;
        delete [] dervForZ;
        dervForZ = NULL;
        delete [] dervBacX;
        dervBacX = NULL;
        delete [] dervBacY;
        dervBacY = NULL;
        delete [] dervBacZ;
        dervBacZ = NULL;
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
        delete [] sumForX;
        sumForX = NULL;
        delete [] sumBacX;
        sumBacX = NULL;
        delete [] sumForY;
        sumForY = NULL;
        delete [] sumBacY;
        sumBacY = NULL;
        delete [] sumForZ;
        sumForZ = NULL;
        delete [] sumBacZ;
        sumBacZ = NULL;

        l2normError = l2Norm(tempImgArrayPrev,imageArr,imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";
    }
    mserror = mse(OrigImgArray,imageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(OrigImgArray,imageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] d11;
    d11 = NULL;
    delete [] d12;
    d12 = NULL;
    delete [] d13;
    d13 = NULL;
    delete [] d21;
    d21 = NULL;
    delete [] d22;
    d22 = NULL;
    delete [] d23;
    d23 = NULL;
    delete [] d31;
    d31 = NULL;
    delete [] d32;
    d32 = NULL;
    delete [] d33;
    d33 = NULL;

    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;
    delete [] OrigImgArray;
    OrigImgArray = NULL;


    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
