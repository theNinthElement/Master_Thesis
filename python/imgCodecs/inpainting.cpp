#include <Python.h>
#include "derivatives.h"
#include "supplementary_functions.h"
#include "fed.h"
#include "inpainting.h"

#include <QDebug>
#include <algorithm>
#include <thread>
#include <mutex>
#include <cmath>
#include <QQueue>
#include <iostream>
#include <QFile>
#include <QDir>
#include <QMessageBox>
//#include <ctpl_stl.h>

std::mutex myMutex;

//***************************************************************************************************************************************************************************************************************************//
//~ Edge Enhancing Diffusion-based Inpaiting by Explicit Scheme
void eed_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight)
{
    clock_t begin = clock();

    int kernelSize = 3, stepSize = 1, randArrTraceIndex;
    double sigma = 1.0, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *dervY, *dervX, *dervXConv, *dervYConv;
    double *dervYD21, *dervXD12, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY;
    double *tempImgArrayPrev;

    //Memory allocation for diffusion tensor entries dij.
    double* d11 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d12 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d21 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d22 = new double[sizeof(double) * imgWidth*imgHeight];

    //Memory allocation for temporary image array.
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    while(l2normError > tol) {
        outConvX = convolutionX(scatImageArr, imgWidth, imgHeight, gausKernel, kernelSize);    // Calculate the convolution on x axis.
        outConvXY = convolutionY(outConvX, imgWidth, imgHeight, gausKernel, kernelSize);     // Calculate the convolution on y axis after x axis.
        //First derivatives
        dervXConv = derivativeX(outConvXY,imgWidth,imgHeight,stepSize);                         // Derivative of convolved image w.r.t. x.
        dervYConv = derivativeY(outConvXY,imgWidth,imgHeight,stepSize);                         // Derivative of convolved image w.r.t. y.
        //First derivatives
        dervX = derivativeX(scatImageArr,imgWidth,imgHeight,stepSize);                         // Derivative w.r.t. x.
        dervY = derivativeY(scatImageArr,imgWidth,imgHeight,stepSize);                         // Derivative w.r.t. y.

        dervForX = dervForwX(scatImageArr, imgWidth, imgHeight, stepSize);
        dervForY = dervForwY(scatImageArr, imgWidth, imgHeight, stepSize);
        dervBacX = dervBackX(scatImageArr, imgWidth, imgHeight, stepSize);
        dervBacY = dervBackY(scatImageArr, imgWidth, imgHeight, stepSize);

        for(int i=0; i<imgWidth*imgHeight; i++) {
            tempImgArrayPrev[i] = scatImageArr[i];

            //Define eigenvectors and entries for the diffusion tensor.
            double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i];
            double normi = sqrt(norm_i_square);
            double d1, d2;
            double v1[2], v2[2];

            if(normi == 0) {
                v1[0] = 1;
                v1[1] = 0;

                v2[0] = 0;
                v2[1] = 1;
            } else {
                v1[0] = dervXConv[i]/normi;
                v1[1] = dervYConv[i]/normi;

                v2[0] = dervYConv[i]/normi;
                v2[1] = -dervXConv[i]/normi;
            }

            d1 = charbonnier_diff(norm_i_square, 0.1);
            d2 = 1;

            //Diffusion tensor entries definition with eigenvalues d1, d2 and orthgonal orthonormal eigenvectors v1 and v2.
            d11[i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0];
            d12[i] = (d1*v1[0]*v1[1] + d2*v2[0]*v2[1])*dervY[i];
            d21[i] = (d1*v1[1]*v1[0] + d2*v2[1]*v2[0])*dervX[i];
            d22[i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1];
        }
        sumForX = sumForwX(d11,imgWidth,imgHeight);
        sumForY = sumForwY(d22,imgWidth,imgHeight);
        sumBacX = sumBackX(d11,imgWidth,imgHeight);
        sumBacY = sumBackY(d22,imgWidth,imgHeight);

        dervXD12 = derivativeX(d12,imgWidth,imgHeight,1);
        dervYD21 = derivativeY(d21,imgWidth,imgHeight,1);

        randArrTraceIndex = 0;
        for(int i=0; i<imgWidth*imgHeight; i++) {
            if(i == randPxls[randArrTraceIndex]) {
                randArrTraceIndex++;
                continue;
            }
            scatImageArr[i] = scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i])/2 + dervXD12[i] + dervYD21[i]);
        }

        delete [] dervX;
        dervX = NULL;
        delete [] dervY;
        dervY = NULL;
        delete [] dervXD12;
        dervXD12 = NULL;
        delete [] dervYD21;
        dervYD21 = NULL;
        delete [] dervForX;
        dervForX = NULL;
        delete [] dervForY;
        dervForY = NULL;
        delete [] dervBacX;
        dervBacX = NULL;
        delete [] dervBacY;
        dervBacY = NULL;
        delete [] outConvX;
        outConvX = NULL;
        delete [] outConvXY;
        outConvXY = NULL;
        delete [] dervXConv;
        dervXConv = NULL;
        delete [] dervYConv;
        dervYConv = NULL;
        delete [] sumForX;
        sumForX = NULL;
        delete [] sumBacX;
        sumBacX = NULL;
        delete [] sumForY;
        sumForY = NULL;
        delete [] sumBacY;
        sumBacY = NULL;

        l2normError = l2Norm(tempImgArrayPrev,scatImageArr,imgWidth*imgHeight);
        qDebug() << l2normError << "\n";
    }
    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] d11;
    d11 = NULL;
    delete [] d12;
    d12 = NULL;
    delete [] d21;
    d21 = NULL;
    delete [] d22;
    d22 = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//~ Edge Enhancing Diffusion-based Inpaiting by Fast Explicit Scheme (Steps)
void eed_inpaintingFED_steps(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight)
{
    clock_t begin = clock();

    int kernelSize = 3, stepSize = 1;  // for sigma=1.0 it is better to choose kernelsize=7,  we have 0.05% of the curve’s area outside the discrete kernel. (For different sigma, refer to: http://dev.theomader.com/gaussian-kernel-calculator/)
    int   N;        /* Number of steps                               */
    int randArrTraceIndex;
    float *tau;                                                       /* Vector of FED time step sizes                 */
    double sigma = 1.0;
    double* outConvX;
    double* outConvXY;
    double* dervX;
    double* dervY;
    double* dervXConv;
    double* dervYConv;
    double* dervXD12;
    double* dervYD21;
    double gausKernel[kernelSize];
    double mserror;
//    double aaerror;

    double* dervForX;
    double* dervBacX;
    double* dervForY;
    double* dervBacY;
    double* sumForX;
    double* sumBacX;
    double* sumForY;
    double* sumBacY;

    //Memory allocation for diffusion tensor entries dij.
    double* d11 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d12 = new double[sizeof(double) * imgWidth*imgHeight];
//    double* d21 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d22 = new double[sizeof(double) * imgWidth*imgHeight];

    // Decleartion and memory allocation
    double* tempImgArray = new double[sizeof(double) * imgWidth*imgHeight];

    double* a = new double[sizeof(double) * imgWidth*imgHeight];
    double* b = new double[sizeof(double) * imgWidth*imgHeight];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    /* Initialise step sizes for process with                        */
    /* - number of steps per cycle n:              10                */
    /* - stability limit for 2-D diffusion:        0.25              */
    N = fed_tau_by_steps(numSteps, timeStepSize, 1, &tau);

    float s=0;
    printf("%d\n", N);
    for(int i=0; i<N; i++) {
        s = s + tau[i];
        printf("%f\n", tau[i]);
    }
    printf("Total time step: %f\n", s);

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "Error:" << mserror;
//    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
//    printf("Error: %lf\n", aaerror);

    while(mserror > tol) {
        outConvX = convolutionX(scatImageArr, imgWidth, imgHeight, gausKernel, kernelSize);    // Calculate the convolution on x axis.
        outConvXY = convolutionY(outConvX, imgWidth, imgHeight, gausKernel, kernelSize);     // Calculate the convolution on y axis after x axis.

        //First derivatives
        dervXConv = derivativeX(outConvXY,imgWidth,imgHeight,stepSize);                         // Derivative of convolved image w.r.t. x.
        dervYConv = derivativeY(outConvXY,imgWidth,imgHeight,stepSize);                         // Derivative of convolved image w.r.t. y.

        for(int i=0; i<imgWidth*imgHeight; i++) {
            //Define eigenvectors and entries for the diffusion tensor.
            double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i];
            double normi = sqrt(norm_i_square);
            double d1, d2;
            double v1[2], v2[2];

            if(normi == 0) {
                v1[0] = 1;
                v1[1] = 0;

                v2[0] = 0;
                v2[1] = 1;
            } else {
                v1[0] = dervXConv[i]/normi;
                v1[1] = dervYConv[i]/normi;

                v2[0] = dervYConv[i]/normi;
                v2[1] = -dervXConv[i]/normi;
            }

            d1 = charbonnier_diff(norm_i_square, 0.1);
//                d1 = diffusivity_function(norm_i_square, 3.5, 4, 3.31488);
            d2 = 1;

            //Diffusion tensor entries definition with eigenvalues d1, d2 and orthgonal orthonormal eigenvectors v1 and v2.
            d11[i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0];
            d12[i] = d1*v1[0]*v1[1] + d2*v2[0]*v2[1];
//            d21[i] = d1*v1[1]*v1[0] + d2*v2[1]*v2[0];
            d22[i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1];
        }
        sumForX = sumForwX(d11,imgWidth,imgHeight);
        sumForY = sumForwY(d22,imgWidth,imgHeight);
        sumBacX = sumBackX(d11,imgWidth,imgHeight);
        sumBacY = sumBackY(d22,imgWidth,imgHeight);

        for(int n=0; n<N; n++) {
            //First derivatives
            dervX = derivativeX(scatImageArr,imgWidth,imgHeight,stepSize);                         // Derivative w.r.t. x.
            dervY = derivativeY(scatImageArr,imgWidth,imgHeight,stepSize);                         // Derivative w.r.t. y.

            for(int i=0; i<imgWidth*imgHeight; i++) {
                a[i] = d12[i]*dervY[i];
                b[i] = d12[i]*dervX[i];
            }

            dervXD12 = derivativeX(a,imgWidth,imgHeight,1);
            dervYD21 = derivativeY(b,imgWidth,imgHeight,1);

            dervForX = dervForwX(scatImageArr, imgWidth, imgHeight, stepSize);
            dervForY = dervForwY(scatImageArr, imgWidth, imgHeight, stepSize);
            dervBacX = dervBackX(scatImageArr, imgWidth, imgHeight, stepSize);
            dervBacY = dervBackY(scatImageArr, imgWidth, imgHeight, stepSize);

            randArrTraceIndex = 0;
            for(int i=0; i<imgWidth*imgHeight; i++) {
                if(i == randPxls[randArrTraceIndex]) {
                    randArrTraceIndex++;
                    continue;
                }
                scatImageArr[i] = scatImageArr[i] + tau[n]*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i])/2 + dervXD12[i] + dervYD21[i]);
            }
            delete [] dervX;
            dervX = NULL;
            delete [] dervY;
            dervY = NULL;
            delete [] dervXD12;
            dervXD12 = NULL;
            delete [] dervYD21;
            dervYD21 = NULL;
            delete [] dervForX;
            dervForX = NULL;
            delete [] dervForY;
            dervForY = NULL;
            delete [] dervBacX;
            dervBacX = NULL;
            delete [] dervBacY;
            dervBacY = NULL;
        }
        delete [] outConvX;
        outConvX = NULL;
        delete [] outConvXY;
        outConvXY = NULL;
        delete [] dervXConv;
        dervXConv = NULL;
        delete [] dervYConv;
        dervYConv = NULL;
        delete [] sumForX;
        sumForX = NULL;
        delete [] sumBacX;
        sumBacX = NULL;
        delete [] sumForY;
        sumForY = NULL;
        delete [] sumBacY;
        sumBacY = NULL;

        mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
        qDebug() << "Error:" << mserror;
//        aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
//        printf("Error: %lf\n", aaerror);
    }

    delete [] d11;
    d11 = NULL;
    delete [] d12;
    d12 = NULL;
//    delete [] d21;
//    d21 = NULL;
    delete [] d22;
    d22 = NULL;
    delete [] a;
    a = NULL;
    delete [] b;
    b = NULL;
    delete [] tempImgArray;
    tempImgArray = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//~ Edge Enhancing Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void eed_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight)
{
    clock_t begin = clock();

    int kernelSize = 5, stepSize = 1, randArrTraceIndex, N;
    double sigma = 1.0, gausKernel[kernelSize], mserror, aaerror, l2normError, alpha;
    double *outConvX, *outConvXY, *dervY, *dervX, *dervXConv, *dervYConv;
    double *dervYD21, *dervXD12, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY;
    double *tempImgArrayPrev, *tempImgArrayCurr;

    //Memory allocation for diffusion tensor entries dij.
    double* d11 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d12 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d21 = new double[sizeof(double) * imgWidth*imgHeight];
    double* d22 = new double[sizeof(double) * imgWidth*imgHeight];

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    N = numSteps;
    while(l2normError > tol) {
        for(int n=0; n<N; n++) {
            outConvX = convolutionX(scatImageArr, imgWidth, imgHeight, gausKernel, kernelSize);    // Calculate the convolution on x axis.
            outConvXY = convolutionY(outConvX, imgWidth, imgHeight, gausKernel, kernelSize);     // Calculate the convolution on y axis after x axis.
            //First derivatives
            dervXConv = derivativeX(outConvXY,imgWidth,imgHeight,stepSize);                         // Derivative of convolved image w.r.t. x.
            dervYConv = derivativeY(outConvXY,imgWidth,imgHeight,stepSize);                         // Derivative of convolved image w.r.t. y.
            //First derivatives
            dervX = derivativeX(scatImageArr,imgWidth,imgHeight,stepSize);                         // Derivative w.r.t. x.
            dervY = derivativeY(scatImageArr,imgWidth,imgHeight,stepSize);                         // Derivative w.r.t. y.

            dervForX = dervForwX(scatImageArr, imgWidth, imgHeight, stepSize);
            dervForY = dervForwY(scatImageArr, imgWidth, imgHeight, stepSize);
            dervBacX = dervBackX(scatImageArr, imgWidth, imgHeight, stepSize);
            dervBacY = dervBackY(scatImageArr, imgWidth, imgHeight, stepSize);

            for(int i=0; i<imgWidth*imgHeight; i++) {
                if(n == 0) {
                    tempImgArrayPrev[i] = scatImageArr[i];
                }
                tempImgArrayCurr[i] = scatImageArr[i];

                //Define eigenvectors and entries for the diffusion tensor.
                double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i];
                double normi = sqrt(norm_i_square);
                double d1, d2;
                double v1[2], v2[2];

                if(normi == 0) {
                    v1[0] = 1;
                    v1[1] = 0;

                    v2[0] = 0;
                    v2[1] = 1;
                } else {
                    v1[0] = dervXConv[i]/normi;
                    v1[1] = dervYConv[i]/normi;

                    v2[0] = dervYConv[i]/normi;
                    v2[1] = -dervXConv[i]/normi;
                }

                d1 = charbonnier_diff(norm_i_square, 0.1);
                d2 = 1;

                //Diffusion tensor entries definition with eigenvalues d1, d2 and orthgonal orthonormal eigenvectors v1 and v2.
                d11[i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0];
                d12[i] = (d1*v1[0]*v1[1] + d2*v2[0]*v2[1])*dervY[i];
                d21[i] = (d1*v1[1]*v1[0] + d2*v2[1]*v2[0])*dervX[i];
                d22[i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1];
            }
            sumForX = sumForwX(d11,imgWidth,imgHeight);
            sumForY = sumForwY(d22,imgWidth,imgHeight);
            sumBacX = sumBackX(d11,imgWidth,imgHeight);
            sumBacY = sumBackY(d22,imgWidth,imgHeight);

            dervXD12 = derivativeX(d12,imgWidth,imgHeight,1);
            dervYD21 = derivativeY(d21,imgWidth,imgHeight,1);

            alpha = (double)(4*n+2)/(2*n+3);

            randArrTraceIndex = 0;
            for(int i=0; i<imgWidth*imgHeight; i++) {
                if(i == randPxls[randArrTraceIndex]) {
                    randArrTraceIndex++;
                    continue;
                }
                scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i])/2 + dervXD12[i] + dervYD21[i])) + (1-alpha)*tempImgArrayPrev[i];
                tempImgArrayPrev[i] = tempImgArrayCurr[i];
            }
            delete [] dervX;
            dervX = NULL;
            delete [] dervY;
            dervY = NULL;
            delete [] dervXD12;
            dervXD12 = NULL;
            delete [] dervYD21;
            dervYD21 = NULL;
            delete [] dervForX;
            dervForX = NULL;
            delete [] dervForY;
            dervForY = NULL;
            delete [] dervBacX;
            dervBacX = NULL;
            delete [] dervBacY;
            dervBacY = NULL;
            delete [] outConvX;
            outConvX = NULL;
            delete [] outConvXY;
            outConvXY = NULL;
            delete [] dervXConv;
            dervXConv = NULL;
            delete [] dervYConv;
            dervYConv = NULL;
            delete [] sumForX;
            sumForX = NULL;
            delete [] sumBacX;
            sumBacX = NULL;
            delete [] sumForY;
            sumForY = NULL;
            delete [] sumBacY;
            sumBacY = NULL;
        }
        l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight);
        qDebug() << l2normError << "\n";
    }
    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] d11;
    d11 = NULL;
    delete [] d12;
    d12 = NULL;
    delete [] d21;
    d21 = NULL;
    delete [] d22;
    d22 = NULL;

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
void eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask) {
    clock_t begin = clock();

    int kernelSize = 3, randArrTraceIndex, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror, aaerror, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3], ones[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;

    qDebug() << gridSpcX;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0, locsVectLenn = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }
//    int locsVectLenn = locsVect.size();
//    for(int i=0; i<locsVectLenn; i++) {
//        qDebug() << "Updt: " << updatedPxlLocs[i];
//    }



    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

//    //Count number of zeros in the original image
//    int cnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(imageArr[i] == 0) {
//            cnt++;
//        }
//    }
//    qDebug() << "Number of zeros: " << cnt;

//    //Initialization with average mask pixels
//    double vxlSum = 0;
//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == randPxls[randArrTraceIndex]) {
//            vxlSum = vxlSum + imageArr[i];
//            randArrTraceIndex++;
//        }
//        continue;
//    }
//    double maskSize = randArrTraceIndex;
//    qDebug() << "Unknown voxel initialization: " << vxlSum/maskSize;
//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == randPxls[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        scatImageArr[i] = vxlSum/maskSize;
//    }

//    double vxlSum = 0, count = 1;
//    randArrTraceIndex = 0;
//    for(int d=0; d<imgDepth; d++) {
//        for(int i=imgWidth*imgHeight*d; i<imgWidth*imgHeight*(d+1); i++) {
//            if(i == randPxls[randArrTraceIndex]) {
//                vxlSum = vxlSum + scatImageArr[i];
//                count++;
//                randArrTraceIndex++;
//                continue;
//            }
//        }
//        for(int i=imgWidth*imgHeight*d; i<imgWidth*imgHeight*(d+1); i++) {
//            scatImageArr[i] = vxlSum/(count-1);
//        }
//        qDebug() << "Unknown voxel initialization: " << vxlSum/(count-1);
//        vxlSum = 0;
//        count = 1;
//    }
//    //End of initialization process

    //Initialization mask with known neighbors pixels
//    double* arr = new double[sizeof(double) * 10];
//    double* arrCp = new double[sizeof(double) * 10];
//    for(int i=0; i<10; i++) {
//        arr[i] = (i+sqrt(i))/2;
//        qDebug() << arr[i];
//    }

//    qDebug() << "          ";

//    QQueue<double> arrQueue;
//    arrQueue = array2queue(arr, 10);
//    int i = 0;
//    while (!arrQueue.isEmpty()) {
//        arrCp[i] = arrQueue.dequeue();
//        qDebug() << "Q: " << arrCp[i];
//        i++;
//    }


//    int maskSize = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i==randPxls[maskSize]) {
//            maskSize++;
//        }
//    }
//    qDebug() << maskSize;



//    double* initMask;
//    int count=0;
//    initMask = nearNeighInit(scatImageArr, imgWidth,imgHeight,imgDepth, randPxls);
//    for(int i=0; i<7*6*4; i++) {
//        qDebug() << "After Init Mask: " << initMask[i];
//    }
//    return;

//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(scatImageArr[i]-initMask != 0) {
//            count++;
//        }
//    }
//    qDebug() << count;

//    QQueue<double> imgQueue;
//    imgQueue = array2queue(scatImageArr, imgWidth*imgHeight*imgDepth);


//    int i=0;
//    int count=0;
//    while (!imgQueue.isEmpty()) {
//        double temp = imgQueue.dequeue();
// //        qDebug() << temp << scatImageArr[i]-temp;
//        if(scatImageArr[i]-temp != 0) {
//            count++;
//        }
//        i++;
//    }
//    qDebug() << count;
    //End of initialization process

//    //TV calculation
//    int indx = 0;
//    double tvImg = 0;
//    double* outConvX_tv = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
//    double* outConvXY_tv = convolution3DY(outConvX_tv, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
//    double* outConvXYZ_tv = convolution3DZ(outConvXY_tv, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
//    //First convolved derivatives
//    double* dervXConv_tv = derivative3DX(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//    double* dervYConv_tv = derivative3DY(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//    double* dervZConv_tv = derivative3DZ(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[indx]) {                        //i == randPxls[indx]
//            indx++;
//            continue;
//        }
//        tvImg = tvImg + sqrt(dervXConv_tv[i]*dervXConv_tv[i] + dervYConv_tv[i]*dervYConv_tv[i] + dervZConv_tv[i]*dervZConv_tv[i]);
//    }
//    qDebug() << "Total Variation: " << tvImg/(indx) << indx;
//    delete [] outConvX_tv;
//    outConvX_tv = NULL;
//    delete [] outConvXY_tv;
//    outConvXY_tv = NULL;
//    delete [] outConvXYZ_tv;
//    outConvXYZ_tv = NULL;
//    delete [] dervXConv_tv;
//    dervXConv_tv = NULL;
//    delete [] dervYConv_tv;
//    dervYConv_tv = NULL;
//    delete [] dervZConv_tv;
//    dervZConv_tv = NULL;

//    double tvImg = 0;
//    //First convolved derivatives
//    double* dervX_tv = derivative3DX(imageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//    double* dervY_tv = derivative3DY(imageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//    double* dervZ_tv = derivative3DZ(imageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        tvImg = tvImg + sqrt(dervX_tv[i]*dervX_tv[i] + dervY_tv[i]*dervY_tv[i] + dervZ_tv[i]*dervZ_tv[i]);
//    }
//    qDebug() << "TV: " << tvImg/(imgWidth*imgHeight*imgDepth);



    //Contrast Paremeter estimation
    double quantile = 90.0;
    double estContPar;
//    estContPar = quantlCriterPM(quantile, imageArr,imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    if(zeros2mask) {
        estContPar = quantlCriterPM4Inpainting(quantile, imageArr, updatedPxlLocs, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    } else {
        estContPar = quantlCriterPM4Inpainting(quantile, imageArr, randPxls, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    }
//    estContPar = histoCriterPM4Inpainting(quantile, imageArr, randPxls, imgWidth,imgHeight,imgDepth);
    estContPar = estContPar/25;
    qDebug() << "Estimated Contrast Parameter: " << estContPar;


    double firstSteps = 3, itrCount = 0;
    if(zeros2mask) {
//        N = numSteps;
        while(l2normError > tol) {
            if(itrCount > 3) {
                N = numSteps;
            } else {
                N = firstSteps;
            }
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                //First convolved derivatives
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

                //First derivatives
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
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

//                    double contPar = 0.65364;
//                    double contPar = 0.1;
                    d1 = charbonnier_diff(norm_i_square, estContPar);
                    d2 = 1;
                    d3 = 1;

//                    qDebug() << "Zeros Included - EED";

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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == updatedPxlLocs[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
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
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";

            itrCount++;
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                //First convolved derivatives
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

    //            //Convolved First derivatives
    //            dervForXConv = dervForw3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
    //            dervForYConv = dervForw3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
    //            dervForZConv = dervForw3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);
    //            dervBacXConv = dervBack3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
    //            dervBacYConv = dervBack3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
    //            dervBacZConv = dervBack3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);

                //First derivatives
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);

    // //             Adaptive contrast parameter selection******************************
    //            double tv_norm = 0, contPar;
    //            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
    //                double norm_i_square, normi;

    //                norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i];
    //                normi = sqrt(norm_i_square);
    //                tv_norm = tv_norm + normi;
    //            }
    //            tv_norm = tv_norm/(imgWidth*imgHeight*imgDepth);
    //            contPar = 0.002*tv_norm;
    // //            qDebug() << tv_norm << contPar;
    // //            ********************************************************************

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);

    //                                    qDebug() << "v1" << v1[0] << v1[1] << v1[2];
    //                                    qDebug() << "v2" << v2[0] << v2[1] << v2[2];
    //                                    qDebug() << "v3" << v3[0] << v3[1] << v3[2];
    //                                    qDebug() << "V12" << v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
    //                                    qDebug() << "V13" << v1[0]*v3[0]+v1[1]*v3[1]+v1[2]*v3[2];
    //                                    qDebug() << "V32" << v3[0]*v2[0]+v3[1]*v2[1]+v3[2]*v2[2];

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


    //                double sm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
    //                if(sm == 1) {
    //                    qDebug() << "One";
    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    //                } else {
    //                    qDebug() << "Not One";
    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    //                }
    //                qDebug() << v2[0]*v2[0] << v2[1]*v2[1] << v2[2]*v2[2] << "Sum: " << v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
    //                qDebug() << v3[0]*v3[0] << v3[1]*v3[1] << v3[2]*v3[2] << "Sum: " << v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];

    //                double tempVal = std::max(dervForXConv[i]*dervBacXConv[i],0.0) + std::max(dervForYConv[i]*dervBacYConv[i],0.0) + std::max(dervForZConv[i]*dervBacZConv[i],0.0);
    //                d1 = charbonnier_diff(tempVal, 0.1);

// //                    double contPar = 0.561268;
                    d1 = charbonnier_diff(norm_i_square, estContPar);
                    d2 = 1;
                    d3 = 1;

//                    //Take inverse of (v1,v2,v3) matrix
//                    double mat[3][3], invmat[3][3], determinant=0, g[3];
//                    mat[0][0] = v1[0];
//                    mat[1][0] = v1[1];
//                    mat[2][0] = v1[2];
//                    mat[0][1] = v2[0];
//                    mat[1][1] = v2[1];
//                    mat[2][1] = v2[2];
//                    mat[0][2] = v3[0];
//                    mat[1][2] = v3[1];
//                    mat[2][2] = v3[2];
//                    for(int i = 0; i < 3; i++) {
//                        determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
//                    }
//                    for(int i = 0; i < 3; i++){
//                        for(int j = 0; j < 3; j++)
//                          invmat[i][j] = ((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/determinant;
//                    }
//                    if(i==79468) {
//                        qDebug() << "Smoothed norm: " << norm_i_square << "Soothed XY norm: " << xy_norm_i_square << "Smoothed Z norm: " << dervZConv[i];

//                        for(int k = 0; k < 3; k++){
//                            qDebug() << "v1: " << v1[k];
//                            qDebug() << "v2: " << v2[k];
//                            qDebug() << "v3: " << v3[k];
//                        }
//                        for(int k = 0; k < 3; k++){
//                            qDebug() << "\n";
//                            for(int j = 0; j < 3; j++)
//                                qDebug() << mat[k][j];
//                        }
//                        qDebug() << "Inverse: " << i;
//                        for(int k = 0; k < 3; k++){
//                            for(int j = 0; j < 3; j++)
//                                qDebug() << invmat[k][j] << "\t";
//                                qDebug() << "\n";
//                        }
//                    }
//                    g[0] = invmat[0][0]*dervX[i] + invmat[0][1]*dervY[i] + invmat[0][2]*dervZ[i];
//                    g[1] = invmat[1][0]*dervX[i] + invmat[1][1]*dervY[i] + invmat[1][2]*dervZ[i];
//                    g[2] = invmat[2][0]*dervX[i] + invmat[2][1]*dervY[i] + invmat[2][2]*dervZ[i];
//                    if(i==79468) {
//                        qDebug() << "Coefficients: " << g[0] << g[1] << g[2];
//                    }
//                    d1 = charbonnier_diff(norm_i_square, estContPar);
//                    if(g[1] == 0) {
//                        d2 = 1;
//                    } else {
//                        d2 = 1/g[1];
//                    }
//                    if(g[2] == 0) {
//                        d3 = 1;
//                    } else {
//                        d3 = 1/g[2];
//                    }
//                    if(norm_i_square == 0 || xy_norm_i_square == 0) {
//                        d2 = 1;
//                        d3 = 1;
//                    } else {
//                        if(g[1] == 0) {
//                            d2 = 1;
//                        } else {
//                            d2 = 1/g[1];
//                        }
//                        if(g[2] == 0) {
//                            d3 = 1;
//                        } else {
//                            d3 = 1/g[2];
//                        }
// //                        d2 = 1/g[1];
// //                        d3 = 1/g[2];
//                    }
//                    if(i==79468) {
//                        qDebug() << "Tensor eigenvalues: " << d1 << d2 << d3;
//                    }

//                    qDebug() << "Zeros Not Included - EED";

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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == randPxls[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
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

    //            delete [] dervForXConv;
    //            dervForXConv = NULL;
    //            delete [] dervForYConv;
    //            dervForYConv = NULL;
    //            delete [] dervForZConv;
    //            dervForZConv = NULL;
    //            delete [] dervBacXConv;
    //            dervBacXConv = NULL;
    //            delete [] dervBacYConv;
    //            dervBacYConv = NULL;
    //            delete [] dervBacZConv;
    //            dervBacZConv = NULL;
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    //Make voxel values nonnegative
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(scatImageArr[i] < 0) {
            scatImageArr[i] = 0;
        }
    }

//    //Residual and its entropy calculated here*************************************************************************
//    //uint16 max residual
//    double* resiImg = new double[sizeof(double) * (imgWidth*imgHeight*imgDepth-locsVectLenn)];
//    double maxVal = *max_element(imageArr , imageArr + imgWidth*imgHeight*imgDepth);
//    randArrTraceIndex = 0;
//    int indCnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        resiImg[indCnt] = round(scatImageArr[i] - imageArr[i]);
//        if(resiImg[indCnt] < 0) {
//            resiImg[indCnt] = resiImg[indCnt] + maxVal;
//        }
//        indCnt++;
//    }
//    //Entropy of residual
//    double entrpResi = entropy(resiImg, (imgWidth*imgHeight*imgDepth-locsVectLenn));
//    qDebug() << "Residual Entropy: " << entrpResi << "Max: " << maxVal;
//    delete [] resiImg;
//    resiImg = NULL;
//    //End***************************************************************************************************************

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";


    if(zeros2mask) {
        mserror = mse_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,updatedPxlLocs);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,updatedPxlLocs);
        qDebug() << "AAE error: " << aaerror << "\n";
    } else {
        mserror = mse_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,randPxls);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,randPxls);
        qDebug() << "AAE error: " << aaerror << "\n";
    }


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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
void spatial_dti_eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask, double* spatialDti) {
    clock_t begin = clock();

    int kernelSize = 3, randArrTraceIndex, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror, aaerror, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;

    qDebug() << gridSpcX;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0, locsVectLenn = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }
//    int locsVectLenn = locsVect.size();
//    for(int i=0; i<locsVectLenn; i++) {
//        qDebug() << "Updt: " << updatedPxlLocs[i];
//    }



    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

//    //Count number of zeros in the original image
//    int cnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(imageArr[i] == 0) {
//            cnt++;
//        }
//    }
//    qDebug() << "Number of zeros: " << cnt;

//    //Initialization with average mask pixels
//    double vxlSum = 0;
//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == randPxls[randArrTraceIndex]) {
//            vxlSum = vxlSum + imageArr[i];
//            randArrTraceIndex++;
//        }
//        continue;
//    }
//    double maskSize = randArrTraceIndex;
//    qDebug() << "Unknown voxel initialization: " << vxlSum/maskSize;
//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == randPxls[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        scatImageArr[i] = vxlSum/maskSize;
//    }

//    double vxlSum = 0, count = 1;
//    randArrTraceIndex = 0;
//    for(int d=0; d<imgDepth; d++) {
//        for(int i=imgWidth*imgHeight*d; i<imgWidth*imgHeight*(d+1); i++) {
//            if(i == randPxls[randArrTraceIndex]) {
//                vxlSum = vxlSum + scatImageArr[i];
//                count++;
//                randArrTraceIndex++;
//                continue;
//            }
//        }
//        for(int i=imgWidth*imgHeight*d; i<imgWidth*imgHeight*(d+1); i++) {
//            scatImageArr[i] = vxlSum/(count-1);
//        }
//        qDebug() << "Unknown voxel initialization: " << vxlSum/(count-1);
//        vxlSum = 0;
//        count = 1;
//    }
//    //End of initialization process

    //Initialization mask with known neighbors pixels
//    double* arr = new double[sizeof(double) * 10];
//    double* arrCp = new double[sizeof(double) * 10];
//    for(int i=0; i<10; i++) {
//        arr[i] = (i+sqrt(i))/2;
//        qDebug() << arr[i];
//    }

//    qDebug() << "          ";

//    QQueue<double> arrQueue;
//    arrQueue = array2queue(arr, 10);
//    int i = 0;
//    while (!arrQueue.isEmpty()) {
//        arrCp[i] = arrQueue.dequeue();
//        qDebug() << "Q: " << arrCp[i];
//        i++;
//    }


//    int maskSize = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i==randPxls[maskSize]) {
//            maskSize++;
//        }
//    }
//    qDebug() << maskSize;



//    double* initMask;
//    int count=0;
//    initMask = nearNeighInit(scatImageArr, imgWidth,imgHeight,imgDepth, randPxls);
//    for(int i=0; i<7*6*4; i++) {
//        qDebug() << "After Init Mask: " << initMask[i];
//    }
//    return;

//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(scatImageArr[i]-initMask != 0) {
//            count++;
//        }
//    }
//    qDebug() << count;

//    QQueue<double> imgQueue;
//    imgQueue = array2queue(scatImageArr, imgWidth*imgHeight*imgDepth);


//    int i=0;
//    int count=0;
//    while (!imgQueue.isEmpty()) {
//        double temp = imgQueue.dequeue();
// //        qDebug() << temp << scatImageArr[i]-temp;
//        if(scatImageArr[i]-temp != 0) {
//            count++;
//        }
//        i++;
//    }
//    qDebug() << count;
    //End of initialization process

//    //TV calculation
//    int indx = 0;
//    double tvImg = 0;
//    double* outConvX_tv = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
//    double* outConvXY_tv = convolution3DY(outConvX_tv, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
//    double* outConvXYZ_tv = convolution3DZ(outConvXY_tv, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
//    //First convolved derivatives
//    double* dervXConv_tv = derivative3DX(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//    double* dervYConv_tv = derivative3DY(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//    double* dervZConv_tv = derivative3DZ(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[indx]) {                        //i == randPxls[indx]
//            indx++;
//            continue;
//        }
//        tvImg = tvImg + sqrt(dervXConv_tv[i]*dervXConv_tv[i] + dervYConv_tv[i]*dervYConv_tv[i] + dervZConv_tv[i]*dervZConv_tv[i]);
//    }
//    qDebug() << "Total Variation: " << tvImg/(indx) << indx;
//    delete [] outConvX_tv;
//    outConvX_tv = NULL;
//    delete [] outConvXY_tv;
//    outConvXY_tv = NULL;
//    delete [] outConvXYZ_tv;
//    outConvXYZ_tv = NULL;
//    delete [] dervXConv_tv;
//    dervXConv_tv = NULL;
//    delete [] dervYConv_tv;
//    dervYConv_tv = NULL;
//    delete [] dervZConv_tv;
//    dervZConv_tv = NULL;

//    double tvImg = 0;
//    //First convolved derivatives
//    double* dervX_tv = derivative3DX(imageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//    double* dervY_tv = derivative3DY(imageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//    double* dervZ_tv = derivative3DZ(imageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        tvImg = tvImg + sqrt(dervX_tv[i]*dervX_tv[i] + dervY_tv[i]*dervY_tv[i] + dervZ_tv[i]*dervZ_tv[i]);
//    }
//    qDebug() << "TV: " << tvImg/(imgWidth*imgHeight*imgDepth);



    //Contrast Paremeter estimation
    double quantile = 90.0;
    double estContPar;
//    estContPar = quantlCriterPM(quantile, imageArr,imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    if(zeros2mask) {
        estContPar = quantlCriterPM4Inpainting(quantile, imageArr, updatedPxlLocs, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    } else {
        estContPar = quantlCriterPM4Inpainting(quantile, imageArr, randPxls, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    }
//    estContPar = histoCriterPM4Inpainting(quantile, imageArr, randPxls, imgWidth,imgHeight,imgDepth);
    estContPar = estContPar/25;
    qDebug() << "Estimated Contrast Parameter: " << estContPar;


    double firstSteps = 3, itrCount = 0;
    if(zeros2mask) {
//        N = numSteps;
        while(l2normError > tol) {
            if(itrCount > 3) {
                N = numSteps;
            } else {
                N = firstSteps;
            }
            for(int n=0; n<N; n++) {
//                outConvX = convolution3DX(spatialDti, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
//                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
//                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
//                //First convolved derivatives
//                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
                //First convolved derivatives
                dervXConv = derivative3DX(spatialDti,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(spatialDti,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(spatialDti,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

                //First derivatives
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
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

//                    double contPar = 0.65364;
//                    double contPar = 0.1;
                    d1 = charbonnier_diff(norm_i_square, estContPar);
                    d2 = 1;
                    d3 = 1;

//                    qDebug() << "Zeros Included - EED";

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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == updatedPxlLocs[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
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
//                delete [] outConvX;
//                outConvX = NULL;
//                delete [] outConvXY;
//                outConvXY = NULL;
//                delete [] outConvXYZ;
//                outConvXYZ = NULL;
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
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";

            itrCount++;
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                //                outConvX = convolution3DX(spatialDti, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                //                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                //                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                //                //First convolved derivatives
                //                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                //                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                //                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
                                //First convolved derivatives
                dervXConv = derivative3DX(spatialDti,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(spatialDti,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(spatialDti,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
    //            //Convolved First derivatives
    //            dervForXConv = dervForw3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
    //            dervForYConv = dervForw3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
    //            dervForZConv = dervForw3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);
    //            dervBacXConv = dervBack3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
    //            dervBacYConv = dervBack3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
    //            dervBacZConv = dervBack3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);

                //First derivatives
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);

    // //             Adaptive contrast parameter selection******************************
    //            double tv_norm = 0, contPar;
    //            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
    //                double norm_i_square, normi;

    //                norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i];
    //                normi = sqrt(norm_i_square);
    //                tv_norm = tv_norm + normi;
    //            }
    //            tv_norm = tv_norm/(imgWidth*imgHeight*imgDepth);
    //            contPar = 0.002*tv_norm;
    // //            qDebug() << tv_norm << contPar;
    // //            ********************************************************************

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);

    //                                    qDebug() << "v1" << v1[0] << v1[1] << v1[2];
    //                                    qDebug() << "v2" << v2[0] << v2[1] << v2[2];
    //                                    qDebug() << "v3" << v3[0] << v3[1] << v3[2];
    //                                    qDebug() << "V12" << v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
    //                                    qDebug() << "V13" << v1[0]*v3[0]+v1[1]*v3[1]+v1[2]*v3[2];
    //                                    qDebug() << "V32" << v3[0]*v2[0]+v3[1]*v2[1]+v3[2]*v2[2];

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


    //                double sm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
    //                if(sm == 1) {
    //                    qDebug() << "One";
    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    //                } else {
    //                    qDebug() << "Not One";
    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    //                }
    //                qDebug() << v2[0]*v2[0] << v2[1]*v2[1] << v2[2]*v2[2] << "Sum: " << v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
    //                qDebug() << v3[0]*v3[0] << v3[1]*v3[1] << v3[2]*v3[2] << "Sum: " << v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];

    //                double tempVal = std::max(dervForXConv[i]*dervBacXConv[i],0.0) + std::max(dervForYConv[i]*dervBacYConv[i],0.0) + std::max(dervForZConv[i]*dervBacZConv[i],0.0);
    //                d1 = charbonnier_diff(tempVal, 0.1);

//                    double contPar = 0.561268;
                    d1 = charbonnier_diff(norm_i_square, estContPar);
                    d2 = 1;
                    d3 = 1;

//                    qDebug() << "Zeros Not Included - EED";

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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == randPxls[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
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
//                delete [] outConvX;
//                outConvX = NULL;
//                delete [] outConvXY;
//                outConvXY = NULL;
//                delete [] outConvXYZ;
//                outConvXYZ = NULL;
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

    //            delete [] dervForXConv;
    //            dervForXConv = NULL;
    //            delete [] dervForYConv;
    //            dervForYConv = NULL;
    //            delete [] dervForZConv;
    //            dervForZConv = NULL;
    //            delete [] dervBacXConv;
    //            dervBacXConv = NULL;
    //            delete [] dervBacYConv;
    //            dervBacYConv = NULL;
    //            delete [] dervBacZConv;
    //            dervBacZConv = NULL;
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    //Make voxel values nonnegative
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(scatImageArr[i] < 0) {
            scatImageArr[i] = 0;
        }
    }

//    //Residual and its entropy calculated here*************************************************************************
//    //uint16 max residual
//    double* resiImg = new double[sizeof(double) * (imgWidth*imgHeight*imgDepth-locsVectLenn)];
//    double maxVal = *max_element(imageArr , imageArr + imgWidth*imgHeight*imgDepth);
//    randArrTraceIndex = 0;
//    int indCnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        resiImg[indCnt] = round(scatImageArr[i] - imageArr[i]);
//        if(resiImg[indCnt] < 0) {
//            resiImg[indCnt] = resiImg[indCnt] + maxVal;
//        }
//        indCnt++;
//    }
//    //Entropy of residual
//    double entrpResi = entropy(resiImg, (imgWidth*imgHeight*imgDepth-locsVectLenn));
//    qDebug() << "Residual Entropy: " << entrpResi << "Max: " << maxVal;
//    delete [] resiImg;
//    resiImg = NULL;
//    //End***************************************************************************************************************

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    if(zeros2mask) {
        mserror = mse_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,updatedPxlLocs);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,updatedPxlLocs);
        qDebug() << "AAE error: " << aaerror << "\n";
    } else {
        mserror = mse_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,randPxls);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,randPxls);
        qDebug() << "AAE error: " << aaerror << "\n";
    }

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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
void dti_eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask, double **eigVals, double ***eigVecs) {
    clock_t begin = clock();

    int kernelSize = 3, randArrTraceIndex, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror, aaerror, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;

    qDebug() << gridSpcX;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0, locsVectLenn = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }
//    int locsVectLenn = locsVect.size();
//    for(int i=0; i<locsVectLenn; i++) {
//        qDebug() << "Updt: " << updatedPxlLocs[i];
//    }



    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
//    l2normError = 1000000000;
    qDebug() << "l2 norm error: " << l2normError << "\n";

//    //Count number of zeros in the original image
//    int cnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(imageArr[i] == 0) {
//            cnt++;
//        }
//    }
//    qDebug() << "Number of zeros: " << cnt;

//    //Initialization with average mask pixels
//    double vxlSum = 0;
//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == randPxls[randArrTraceIndex]) {
//            vxlSum = vxlSum + imageArr[i];
//            randArrTraceIndex++;
//        }
//        continue;
//    }
//    double maskSize = randArrTraceIndex;
//    qDebug() << "Unknown voxel initialization: " << vxlSum/maskSize;
//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == randPxls[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        scatImageArr[i] = vxlSum/maskSize;
//    }

//    double vxlSum = 0, count = 1;
//    randArrTraceIndex = 0;
//    for(int d=0; d<imgDepth; d++) {
//        for(int i=imgWidth*imgHeight*d; i<imgWidth*imgHeight*(d+1); i++) {
//            if(i == randPxls[randArrTraceIndex]) {
//                vxlSum = vxlSum + scatImageArr[i];
//                count++;
//                randArrTraceIndex++;
//                continue;
//            }
//        }
//        for(int i=imgWidth*imgHeight*d; i<imgWidth*imgHeight*(d+1); i++) {
//            scatImageArr[i] = vxlSum/(count-1);
//        }
//        qDebug() << "Unknown voxel initialization: " << vxlSum/(count-1);
//        vxlSum = 0;
//        count = 1;
//    }
//    //End of initialization process

    //Initialization mask with known neighbors pixels
//    double* arr = new double[sizeof(double) * 10];
//    double* arrCp = new double[sizeof(double) * 10];
//    for(int i=0; i<10; i++) {
//        arr[i] = (i+sqrt(i))/2;
//        qDebug() << arr[i];
//    }

//    qDebug() << "          ";

//    QQueue<double> arrQueue;
//    arrQueue = array2queue(arr, 10);
//    int i = 0;
//    while (!arrQueue.isEmpty()) {
//        arrCp[i] = arrQueue.dequeue();
//        qDebug() << "Q: " << arrCp[i];
//        i++;
//    }


//    int maskSize = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i==randPxls[maskSize]) {
//            maskSize++;
//        }
//    }
//    qDebug() << maskSize;



//    double* initMask;
//    int count=0;
//    initMask = nearNeighInit(scatImageArr, imgWidth,imgHeight,imgDepth, randPxls);
//    for(int i=0; i<7*6*4; i++) {
//        qDebug() << "After Init Mask: " << initMask[i];
//    }
//    return;

//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(scatImageArr[i]-initMask != 0) {
//            count++;
//        }
//    }
//    qDebug() << count;

//    QQueue<double> imgQueue;
//    imgQueue = array2queue(scatImageArr, imgWidth*imgHeight*imgDepth);


//    int i=0;
//    int count=0;
//    while (!imgQueue.isEmpty()) {
//        double temp = imgQueue.dequeue();
// //        qDebug() << temp << scatImageArr[i]-temp;
//        if(scatImageArr[i]-temp != 0) {
//            count++;
//        }
//        i++;
//    }
//    qDebug() << count;
    //End of initialization process

//    //TV calculation
//    int indx = 0;
//    double tvImg = 0;
//    double* outConvX_tv = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
//    double* outConvXY_tv = convolution3DY(outConvX_tv, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
//    double* outConvXYZ_tv = convolution3DZ(outConvXY_tv, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
//    //First convolved derivatives
//    double* dervXConv_tv = derivative3DX(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//    double* dervYConv_tv = derivative3DY(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//    double* dervZConv_tv = derivative3DZ(outConvXYZ_tv,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[indx]) {                        //i == randPxls[indx]
//            indx++;
//            continue;
//        }
//        tvImg = tvImg + sqrt(dervXConv_tv[i]*dervXConv_tv[i] + dervYConv_tv[i]*dervYConv_tv[i] + dervZConv_tv[i]*dervZConv_tv[i]);
//    }
//    qDebug() << "Total Variation: " << tvImg/(indx) << indx;
//    delete [] outConvX_tv;
//    outConvX_tv = NULL;
//    delete [] outConvXY_tv;
//    outConvXY_tv = NULL;
//    delete [] outConvXYZ_tv;
//    outConvXYZ_tv = NULL;
//    delete [] dervXConv_tv;
//    dervXConv_tv = NULL;
//    delete [] dervYConv_tv;
//    dervYConv_tv = NULL;
//    delete [] dervZConv_tv;
//    dervZConv_tv = NULL;

//    double tvImg = 0;
//    //First convolved derivatives
//    double* dervX_tv = derivative3DX(imageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//    double* dervY_tv = derivative3DY(imageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//    double* dervZ_tv = derivative3DZ(imageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        tvImg = tvImg + sqrt(dervX_tv[i]*dervX_tv[i] + dervY_tv[i]*dervY_tv[i] + dervZ_tv[i]*dervZ_tv[i]);
//    }
//    qDebug() << "TV: " << tvImg/(imgWidth*imgHeight*imgDepth);



    //Contrast Paremeter estimation
    double quantile = 90.0;
    double estContPar;
//    estContPar = quantlCriterPM(quantile, imageArr,imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    if(zeros2mask) {
        estContPar = quantlCriterPM4Inpainting(quantile, imageArr, updatedPxlLocs, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    } else {
        estContPar = quantlCriterPM4Inpainting(quantile, imageArr, randPxls, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
    }
//    estContPar = histoCriterPM4Inpainting(quantile, imageArr, randPxls, imgWidth,imgHeight,imgDepth);
    estContPar = estContPar/25;
    qDebug() << "Estimated Contrast Parameter: " << estContPar;


//    for(int i=0; i<21; i++) {
//        qDebug() << eigVals[i][0] << eigVals[i][1] << eigVals[i][2] << imageArr[i];
//    }

//    for(int i=0; i<54501; i++) {
//        for(int j=0; j<3; j++) {
//            for(int k=0; k<3; k++) {
//                qDebug() << i << j << k << eigVecs[i][j][k];
//            }
//            qDebug() << "-----------";
//        }
//        qDebug() << "********************************";
//    }


    double firstSteps = 3, itrCount = 0;
    if(zeros2mask) {
//        N = numSteps;
        while(l2normError > tol) {
            if(itrCount > 3) {
                N = numSteps;
            } else {
                N = firstSteps;
            }
            for(int n=0; n<N; n++) {
//                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
//                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
//                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
//                //First convolved derivatives
//                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

                //First derivatives
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double d1, d2, d3;
//                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
//                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
//                    if(norm_i_square == 0) {
//                        v1[0] = 1;
//                        v1[1] = 0;
//                        v1[2] = 0;

//                        v2[0] = 0;
//                        v2[1] = 1;
//                        v2[2] = 0;

//                        v3[0] = 0;
//                        v3[1] = 0;
//                        v3[2] = 1;
//                    } else if(xy_norm_i_square == 0) {
//                        v1[0] = 0;
//                        v1[1] = 0;
//                        v1[2] = dervZConv[i]/normi;

//                        v2[0] = 0;
//                        v2[1] = 0;
//                        v2[2] = 0;

//                        v3[0] = 0;
//                        v3[1] = 0;
//                        v3[2] = 0;
//                    } else {
//                        v1[0] = dervXConv[i]/normi;
//                        v1[1] = dervYConv[i]/normi;
//                        v1[2] = dervZConv[i]/normi;

//                        v2[0] = dervYConv[i]/xy_norm_i;
//                        v2[1] = -dervXConv[i]/xy_norm_i;
//                        v2[2] = 0;

//                        v3[0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
//                        v3[1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
//                        v3[2] = -xy_norm_i_square/(xy_norm_i*normi);
//                    }
                    v1[0] = eigVecs[i][0][0];
                    v1[1] = eigVecs[i][0][1];
                    v1[2] = eigVecs[i][0][2];
                    v2[0] = eigVecs[i][1][0];
                    v2[1] = eigVecs[i][1][1];
                    v2[2] = eigVecs[i][1][2];
                    v3[0] = eigVecs[i][2][0];
                    v3[1] = eigVecs[i][2][1];
                    v3[2] = eigVecs[i][2][2];

// //                    double contPar = 0.65364;
// //                    double contPar = 0.1;
//                    d1 = charbonnier_diff(norm_i_square, estContPar);
//                    d2 = 1;
//                    d3 = 1;
                    d1 = eigVals[i][0];
                    d2 = eigVals[i][1];
                    d3 = eigVals[i][2];

//                    qDebug() << "Zeros Included - EED";

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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == updatedPxlLocs[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
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
//                delete [] outConvX;
//                outConvX = NULL;
//                delete [] outConvXY;
//                outConvXY = NULL;
//                delete [] outConvXYZ;
//                outConvXYZ = NULL;
//                delete [] dervXConv;
//                dervXConv = NULL;
//                delete [] dervYConv;
//                dervYConv = NULL;
//                delete [] dervZConv;
//                dervZConv = NULL;
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
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";

            itrCount++;
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
//                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
//                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
//                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
//                //First convolved derivatives
//                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
//                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
//                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

//    //            //Convolved First derivatives
//    //            dervForXConv = dervForw3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
//    //            dervForYConv = dervForw3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
//    //            dervForZConv = dervForw3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);
//    //            dervBacXConv = dervBack3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
//    //            dervBacYConv = dervBack3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
//    //            dervBacZConv = dervBack3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);

                //First derivatives
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, gridSpcZ);

    // //             Adaptive contrast parameter selection******************************
    //            double tv_norm = 0, contPar;
    //            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
    //                double norm_i_square, normi;

    //                norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i];
    //                normi = sqrt(norm_i_square);
    //                tv_norm = tv_norm + normi;
    //            }
    //            tv_norm = tv_norm/(imgWidth*imgHeight*imgDepth);
    //            contPar = 0.002*tv_norm;
    // //            qDebug() << tv_norm << contPar;
    // //            ********************************************************************

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double d1, d2, d3;
//                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
//                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);

    //                                    qDebug() << "v1" << v1[0] << v1[1] << v1[2];
    //                                    qDebug() << "v2" << v2[0] << v2[1] << v2[2];
    //                                    qDebug() << "v3" << v3[0] << v3[1] << v3[2];
    //                                    qDebug() << "V12" << v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
    //                                    qDebug() << "V13" << v1[0]*v3[0]+v1[1]*v3[1]+v1[2]*v3[2];
    //                                    qDebug() << "V32" << v3[0]*v2[0]+v3[1]*v2[1]+v3[2]*v2[2];

//                    if(norm_i_square == 0) {
//                        v1[0] = 1;
//                        v1[1] = 0;
//                        v1[2] = 0;

//                        v2[0] = 0;
//                        v2[1] = 1;
//                        v2[2] = 0;

//                        v3[0] = 0;
//                        v3[1] = 0;
//                        v3[2] = 1;
//                    } else if(xy_norm_i_square == 0) {
//                        v1[0] = 0;
//                        v1[1] = 0;
//                        v1[2] = dervZConv[i]/normi;

//                        v2[0] = 0;
//                        v2[1] = 0;
//                        v2[2] = 0;

//                        v3[0] = 0;
//                        v3[1] = 0;
//                        v3[2] = 0;
//                    } else {
//                        v1[0] = dervXConv[i]/normi;
//                        v1[1] = dervYConv[i]/normi;
//                        v1[2] = dervZConv[i]/normi;

//                        v2[0] = dervYConv[i]/xy_norm_i;
//                        v2[1] = -dervXConv[i]/xy_norm_i;
//                        v2[2] = 0;

//                        v3[0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
//                        v3[1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
//                        v3[2] = -xy_norm_i_square/(xy_norm_i*normi);
//                    }
                    v1[0] = eigVecs[i][0][0];
                    v1[1] = eigVecs[i][1][0];
                    v1[2] = eigVecs[i][2][0];
                    v2[0] = eigVecs[i][0][1];
                    v2[1] = eigVecs[i][1][1];
                    v2[2] = eigVecs[i][2][1];
                    v3[0] = eigVecs[i][0][2];
                    v3[1] = eigVecs[i][1][2];
                    v3[2] = eigVecs[i][2][2];

//                    if(abs(sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])-1)>0.001 || abs(sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]))>0.001 || abs(sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]))>0.001) {
//                        qDebug() << sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]) << sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]) << sqrt(v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2]);
//                    }


//    //                double sm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
//    //                if(sm == 1) {
//    //                    qDebug() << "One";
//    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
//    //                } else {
//    //                    qDebug() << "Not One";
//    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
//    //                }
//    //                qDebug() << v2[0]*v2[0] << v2[1]*v2[1] << v2[2]*v2[2] << "Sum: " << v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
//    //                qDebug() << v3[0]*v3[0] << v3[1]*v3[1] << v3[2]*v3[2] << "Sum: " << v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];

//    //                double tempVal = std::max(dervForXConv[i]*dervBacXConv[i],0.0) + std::max(dervForYConv[i]*dervBacYConv[i],0.0) + std::max(dervForZConv[i]*dervBacZConv[i],0.0);
//    //                d1 = charbonnier_diff(tempVal, 0.1);

// //                    double contPar = 0.561268;
//                    d1 = charbonnier_diff(norm_i_square, estContPar);
//                    d2 = 1;
//                    d3 = 1;
//                    d1 = 1 - charbonnier_diff(eigVals[i][0], estContPar);
//                    d2 = 1 - charbonnier_diff(eigVals[i][1], estContPar);
//                    d3 = 1 - charbonnier_diff(eigVals[i][2], estContPar);

                    //Eigenvalues resolution change
                    double eigVal1_tmp = eigVals[i][0]*100000000000000;
                    double eigVal2_tmp = eigVals[i][1]*100000000000000;
                    double eigVal3_tmp = eigVals[i][2]*100000000000000;
                    d1 = eigVal1_tmp/sqrt(1+eigVal1_tmp*eigVal1_tmp);
                    d2 = d1*eigVal2_tmp/eigVal1_tmp;
                    d3 = d1*eigVal3_tmp/eigVal1_tmp;
//                    d2 = eigVal2_tmp/sqrt(1+eigVal2_tmp*eigVal2_tmp);
//                    d3 = eigVal3_tmp/sqrt(1+eigVal3_tmp*eigVal3_tmp);
//                    d1=1, d2=1, d3=1;
//                    //Eigenvalues transformation
//                    if(eigVal1_tmp<0) {
//                        d1 = 0;
//                    } else if(eigVal1_tmp>1) {
//                        d1 = 1;
//                    } else {
//                        d1 = eigVal1_tmp;
//                    }
//                    if(eigVal2_tmp<0) {
//                        d2 = 0;
//                    } else if(eigVal2_tmp>1) {
//                        d2 = 1;
//                    } else {
//                        d2 = eigVal2_tmp;
//                    }
//                    if(eigVal3_tmp<0) {
//                        d3 = 0;
//                    } else if(eigVal3_tmp>1) {
//                        d3 = 1;
//                    } else {
//                        d3 = eigVal3_tmp;
//                    }

                    if(i == 121546) {
                        qDebug() << i << v1[0] << v1[1] << v1[2] << d1 << eigVal1_tmp;
                        qDebug() << i << v2[0] << v2[1] << v2[2] << d2 << eigVal2_tmp;
                        qDebug() << i << v3[0] << v3[1] << v3[2] << d3 << eigVal3_tmp;
                        qDebug() << """""""";
                    }


//                    qDebug() << "Zeros Not Included - EED";

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
//                    d11[i] = d1*v1[0]*v1[0] + d1*v2[0]*v2[0] + d1*v3[0]*v3[0];
//                    d12[i] = (d1*v1[0]*v1[1] + d1*v2[0]*v2[1] + d1*v3[0]*v3[1])*dervY[i];
//                    d13[i] = (d1*v1[0]*v1[2] + d1*v2[0]*v2[2] + d1*v3[0]*v3[2])*dervZ[i];
//                    d21[i] = (d2*v1[1]*v1[0] + d2*v2[1]*v2[0] + d2*v3[1]*v3[0])*dervX[i];
//                    d22[i] = d2*v1[1]*v1[1] + d2*v2[1]*v2[1] + d2*v3[1]*v3[1];
//                    d23[i] = (d2*v1[1]*v1[2] + d2*v2[1]*v2[2] + d2*v3[1]*v3[2])*dervZ[i];
//                    d31[i] = (d3*v1[0]*v1[2] + d3*v2[0]*v2[2] + d3*v3[0]*v3[2])*dervX[i];
//                    d32[i] = (d3*v1[1]*v1[2] + d3*v2[1]*v2[2] + d3*v3[1]*v3[2])*dervY[i];
//                    d33[i] = d3*v1[2]*v1[2] + d3*v2[2]*v2[2] + d3*v3[2]*v3[2];
//                    d11[i] = d1*v1[0]*v1[0] + d2*v1[1]*v1[1] + d3*v1[2]*v1[2];
//                    d12[i] = (d1*v1[0]*v2[1] + d2*v1[1]*v2[1] + d3*v1[2]*v2[2])*dervY[i];
//                    d13[i] = (d1*v1[0]*v3[0] + d2*v1[1]*v3[1] + d3*v1[2]*v3[2])*dervZ[i];
//                    d21[i] = (d1*v1[0]*v1[0] + d2*v1[1]*v1[1] + d3*v1[2]*v1[2])*dervX[i];
//                    d22[i] = d1*v1[0]*v2[0] + d2*v1[1]*v2[1] + d3*v1[2]*v2[2];
//                    d23[i] = (d1*v1[0]*v3[0] + d2*v1[1]*v3[1] + d3*v1[2]*v3[2])*dervZ[i];
//                    d31[i] = (d1*v1[0]*v3[0] + d2*v1[1]*v3[1] + d3*v1[2]*v3[2])*dervX[i];
//                    d32[i] = (d1*v2[0]*v3[0] + d2*v2[1]*v3[1] + d3*v2[2]*v3[1])*dervY[i];
//                    d33[i] = d1*v3[0]*v3[0] + d2*v3[1]*v3[1] + d3*v3[2]*v3[2];
                }
                sumForX = sumForw3DX(d11,imgWidth,imgHeight,imgDepth);
                sumForY = sumForw3DY(d22,imgWidth,imgHeight,imgDepth);
                sumForZ = sumForw3DZ(d33,imgWidth,imgHeight,imgDepth);
                sumBacX = sumBack3DX(d11,imgWidth,imgHeight,imgDepth);
                sumBacY = sumBack3DY(d22,imgWidth,imgHeight,imgDepth);
                sumBacZ = sumBack3DZ(d33,imgWidth,imgHeight,imgDepth);

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == randPxls[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
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
//                delete [] outConvX;
//                outConvX = NULL;
//                delete [] outConvXY;
//                outConvXY = NULL;
//                delete [] outConvXYZ;
//                outConvXYZ = NULL;
//                delete [] dervXConv;
//                dervXConv = NULL;
//                delete [] dervYConv;
//                dervYConv = NULL;
//                delete [] dervZConv;
//                dervZConv = NULL;
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

    //            delete [] dervForXConv;
    //            dervForXConv = NULL;
    //            delete [] dervForYConv;
    //            dervForYConv = NULL;
    //            delete [] dervForZConv;
    //            dervForZConv = NULL;
    //            delete [] dervBacXConv;
    //            dervBacXConv = NULL;
    //            delete [] dervBacYConv;
    //            dervBacYConv = NULL;
    //            delete [] dervBacZConv;
    //            dervBacZConv = NULL;
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    //Make voxel values nonnegative
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(scatImageArr[i] < 0) {
            scatImageArr[i] = 0;
        }
    }

//    //Residual and its entropy calculated here*************************************************************************
//    //uint16 max residual
//    double* resiImg = new double[sizeof(double) * (imgWidth*imgHeight*imgDepth-locsVectLenn)];
//    double maxVal = *max_element(imageArr , imageArr + imgWidth*imgHeight*imgDepth);
//    randArrTraceIndex = 0;
//    int indCnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        resiImg[indCnt] = round(scatImageArr[i] - imageArr[i]);
//        if(resiImg[indCnt] < 0) {
//            resiImg[indCnt] = resiImg[indCnt] + maxVal;
//        }
//        indCnt++;
//    }
//    //Entropy of residual
//    double entrpResi = entropy(resiImg, (imgWidth*imgHeight*imgDepth-locsVectLenn));
//    qDebug() << "Residual Entropy: " << entrpResi << "Max: " << maxVal;
//    delete [] resiImg;
//    resiImg = NULL;
//    //End***************************************************************************************************************

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    if(zeros2mask) {
        mserror = mse_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,updatedPxlLocs);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,updatedPxlLocs);
        qDebug() << "AAE error: " << aaerror << "\n";
    } else {
        mserror = mse_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,randPxls);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region3D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth,randPxls);
        qDebug() << "AAE error: " << aaerror << "\n";
    }

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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
void linear_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask) {
    clock_t begin = clock();

    int stepSize = 1, randArrTraceIndex, N;
    double mserror, aaerror, l2normError, alpha;
    double *dervYY, *dervXX, *dervZZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";



//    randArrTraceIndex = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        scatImageArr[i] = 0;
//        if(i == randPxls[randArrTraceIndex]) {
//            scatImageArr[i] = imageArr[i];
//            randArrTraceIndex++;
//        }
//    }
//    return;

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0, locsVectLenn = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }

    if(zeros2mask) {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                //First derivatives
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    if(i == updatedPxlLocs[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*(dervXX[i]+dervYY[i]+dervZZ[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
                }
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                //First derivatives
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    if(i == randPxls[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] + timeStepSize*(dervXX[i]+dervYY[i]+dervZZ[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
                }
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
            }
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
void eed_3d_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth) {
    clock_t begin = clock();

    int kernelSize = 5, randArrTraceIndex;
    float stepSize = 1;
    double sigma = 1.0, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev;

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

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    while(l2normError > tol) {
        outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);         // Calculate the convolution on x axis.
        outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
        outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
        //First convolved derivatives
        dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
        dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
        dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.
        //First derivatives
        dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
        dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
        dervZ = derivative3DZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. z.

        dervForX = dervForw3DX(scatImageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervForY = dervForw3DY(scatImageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervForZ = dervForw3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervBacX = dervBack3DX(scatImageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervBacY = dervBack3DY(scatImageArr, imgWidth, imgHeight, imgDepth, stepSize);
        dervBacZ = dervBack3DZ(scatImageArr, imgWidth, imgHeight, imgDepth, stepSize);

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

            d1 = charbonnier_diff(norm_i_square, 3.2);
            d2 = 1;
            d3 = 1;


//            //Test
//            if(i == 564478) {
//                qDebug() << "Eig.value 1: " << d1 << v1[0] << v1[1] << v1[2];
//                qDebug() << "Eig.value 2: " << d2 << v2[0] << v2[1] << v2[2];
//                qDebug() << "Eig.value 3: " << d3 << v3[0] << v3[1] << v3[2];
//            }


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

        randArrTraceIndex = 0;
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            tempImgArrayPrev[i] = scatImageArr[i];
            if(i == randPxls[randArrTraceIndex]) {
                randArrTraceIndex++;
                continue;
            }
            scatImageArr[i] = scatImageArr[i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i]);
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

        l2normError = l2Norm(tempImgArrayPrev,scatImageArr,imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";

//        //Test
//        qDebug() << "Original: " << imageArr[564478];
//        qDebug() << "Rec.: " << scatImageArr[564478];
//        break;
    }
    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
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

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
//~ Linear Homogenous Fourth Order Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth)
{
    clock_t begin = clock();

    int randArrTraceIndex;
    float stepSize = 1;
    double *dervXX, *dervYY, *dervZZ, *dervXXXX, *dervYYYY, *dervZZZZ, *dervXXYY, *dervZZXX, *dervYYZZ;

    double mserror, aaerror, l2normError, alpha;
    double* tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << l2normError << "\n";

    while(l2normError > tol) {
        for(int n=0; n<numSteps; n++) {
            //Second and Fourth derivatives
            dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. xx.
            dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. yy.
            dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. zz.
            dervXXXX = derivative3DXX(dervXX,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. xxxx.
            dervYYYY = derivative3DYY(dervYY,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. yyyy.
            dervZZZZ = derivative3DZZ(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. zzzz.
            dervXXYY = derivative3DYY(dervXX,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. xxyy.
            dervZZXX = derivative3DXX(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. zzxx.
            dervYYZZ = derivative3DZZ(dervYY,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. yyzz.

            alpha = (double)(4*n+2)/(2*n+3);

            randArrTraceIndex = 0;
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                if(n == 0) {
                    tempImgArrayPrev[i] = scatImageArr[i];
                }
                tempImgArrayCurr[i] = scatImageArr[i];

                if(i == randPxls[randArrTraceIndex]) {
                    randArrTraceIndex++;
                    continue;
                }
                scatImageArr[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
                tempImgArrayPrev[i] = tempImgArrayCurr[i];
            }

            delete [] dervXX;
            dervXX = NULL;
            delete [] dervYY;
            dervYY = NULL;
            delete [] dervZZ;
            dervZZ = NULL;
            delete [] dervXXYY;
            dervXXYY = NULL;
            delete dervYYZZ;
            dervYYZZ = NULL;
            delete [] dervZZXX;
            dervZZXX = NULL;
            delete [] dervXXXX;
            dervXXXX = NULL;
            delete [] dervYYYY;
            dervYYYY = NULL;
            delete [] dervZZZZ;
            dervZZZZ = NULL;
        }

        l2normError = l2Norm(scatImageArr,tempImgArrayCurr,imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE: " << aaerror << "\n";

    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;
    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
//~ Fourth Order Edge Enhancing Anisotropic Diffusion based Inpaiting by Fast Semi-Iterative Scheme with fourth order diffusion tensor
void foeed_3d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask)
{
    clock_t begin = clock();

    int kernelSize = 3, N, randArrTraceIndex;
    float stepSize = 1;
    double sigma = 1.0, alpha, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervXX, *dervYY, *dervZZ, *dervXY, *dervYZ, *dervX, *dervY, *dervXZ, *dervXConv, *dervYConv, *dervZConv, *t11, *t12, *t13, *t21, *t22, *t23, *t31, *t32, *t33;
    double *d1111, *d1112, *d1113, *d1121, *d1122, *d1123, *d1131, *d1132, *d1133, *d1211, *d1212, *d1213, *d1221, *d1222, *d1223, *d1231, *d1232, *d1233, *d1311, *d1312, *d1313, *d1321, *d1322, *d1323, *d1331, *d1332, *d1333;
    double *d2111, *d2112, *d2113, *d2121, *d2122, *d2123, *d2131, *d2132, *d2133, *d2211, *d2212, *d2213, *d2221, *d2222, *d2223, *d2231, *d2232, *d2233, *d2311, *d2312, *d2313, *d2321, *d2322, *d2323, *d2331, *d2332, *d2333;
    double *d3111, *d3112, *d3113, *d3121, *d3122, *d3123, *d3131, *d3132, *d3133, *d3211, *d3212, *d3213, *d3221, *d3222, *d3223, *d3231, *d3232, *d3233, *d3311, *d3312, *d3313, *d3321, *d3322, *d3323, *d3331, *d3332, *d3333;
    double *dervXXD11, *dervYYD22, *dervZZD33, *dervYD21, *dervXYD21, *dervXD12, *dervYXD12, *dervZD31, *dervXZD31, *dervXD13, *dervZXD13, *dervYD23, *dervZYD23, *dervZD32, *dervYZD32;
    double *tempImgArrayCurr, *tempImgArrayPrev;

    qDebug() << "3D FOEED: " << "Kernel Size: " << kernelSize;

//    double* dt[3][3][3][3];
//    for(int i=0; i<3; i++) {
//        for(int j=0; j<3; j++) {
//            for(int k=0; k<3; k++) {
//                for(int l=0; l<3; l++) {
//                    dt[i][j][k][l] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
//                }
//            }
//        }
//    }

    //************************************************************
    //Memory allocation for 4th order diffusion tensor entries.
    d1111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    //************************************************************
//    double* te[3][3];
//    for(int i=0; i<3; i++) {
//        for(int j=0; j<3; j++) {
//            te[i][j] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
//        }
//    }

    t11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    //Memory allocation for temporary image array.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        int locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }

//    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
//    qDebug() << "MSE error: " << mserror << "\n";
//    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
//    printf("Error: %lf\n", aaerror);
//    qDebug() << "Here";
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    double firstSteps = 3, itrCount = 0;
    if(zeros2mask) {
        N = numSteps;
        while(l2normError > tol) {
            if(itrCount > 4) {
                N = numSteps;
            } else {
                N = firstSteps;
            }
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);    // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on z axis after x and y axis.

                //First derivatives of convolved mask
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.

                //Second derivatives
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
                dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xy.
                dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xz.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
                dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. yz.


                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                    double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                    double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];
    //                double en[6][3][3], vn[3][3];

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

    //                    vn[0][0] = 1;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = 0;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 1;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 1;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = 0;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = dervZConv[i]/normi;;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 0;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 0;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = dervXConv[i]/normi;
    //                    vn[0][1] = dervYConv[i]/normi;
    //                    vn[0][2] = dervZConv[i]/normi;

    //                    vn[1][0] = dervYConv[i]/xy_norm_i;
    //                    vn[1][1] = -dervXConv[i]/xy_norm_i;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][2] = -xy_norm_i_square/(xy_norm_i*normi);

    //                    if(v1[2]-vn[0][2] < 0) {
    //                        qDebug() << v1[2]-vn[0][2];
    //                    }

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
                    }

    //                for(int i=0; i<3; i++) {
    //                    if(v1[i]-vn[0][i] < 0) {
    //                        qDebug() << v1[i]-vn[0][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v2[i]-vn[1][i] < 0) {
    //                        qDebug() << v2[i]-vn[1][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v3[i]-vn[2][i] < 0) {
    //                        qDebug() << v3[i]-vn[2][i];
    //                    }
    //                }


    //                m1 = green_diff(normi, 1);
    //                m1 = aubert_diff(norm_i_square, 1);
    //                m1 = li1(normi, 0.2);
    //                m1 = pm_diff(norm_i_square, 0.2);
    //                m1 = gr_diff(norm_i_square, 1);
    //                m1 = pm_diff2(norm_i_square, 0.2);

                    m1 = charbonnier_diff(norm_i_square, 0.65364);
                    m2 = 1;
                    m3 = 1;
                    m4 = sqrt(m1*m2);
                    //m4 = (m1+m2)/2;
                    //m4 = m1;
                    //m4 = 1;

                    m5 = sqrt(m1*m3);
    //                m5 = 0;

                    m6 = 1;
    //                m6 = 0;

//                   qDebug() << "Zeros Included - FOEED!";

    //                for(int k=0; k<3; k++) {
    //                    for(int i=0; i<3; i++) {
    //                        for(int j=0; j<3; j++) {
    //                            en[k][i][j] = v;
    //                        }
    //                    }
    //                }



    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e1[i][j] = v1[i]*v1[j];
    //                    }
    //                }
                    e1[0][0] = v1[0]*v1[0];
                    e1[0][1] = v1[0]*v1[1];
                    e1[0][2] = v1[0]*v1[2];
                    e1[1][0] = v1[1]*v1[0];
                    e1[1][1] = v1[1]*v1[1];
                    e1[1][2] = v1[1]*v1[2];
                    e1[2][0] = v1[2]*v1[0];
                    e1[2][1] = v1[2]*v1[1];
                    e1[2][2] = v1[2]*v1[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e2[i][j] = v2[i]*v2[j];
    //                    }
    //                }
                    e2[0][0] = v2[0]*v2[0];
                    e2[0][1] = v2[0]*v2[1];
                    e2[0][2] = v2[0]*v2[2];
                    e2[1][0] = v2[1]*v2[0];
                    e2[1][1] = v2[1]*v2[1];
                    e2[1][2] = v2[1]*v2[2];
                    e2[2][0] = v2[2]*v2[0];
                    e2[2][1] = v2[2]*v2[1];
                    e2[2][2] = v2[2]*v2[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e3[i][j] = v3[i]*v3[j];
    //                    }
    //                }
                    e3[0][0] = v3[0]*v3[0];
                    e3[0][1] = v3[0]*v3[1];
                    e3[0][2] = v3[0]*v3[2];
                    e3[1][0] = v3[1]*v3[0];
                    e3[1][1] = v3[1]*v3[1];
                    e3[1][2] = v3[1]*v3[2];
                    e3[2][0] = v3[2]*v3[0];
                    e3[2][1] = v3[2]*v3[1];
                    e3[2][2] = v3[2]*v3[2];

                    e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                    e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                    e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                    e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                    e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                    e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                    e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                    e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                    e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                    e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                    e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                    e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                    e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                    e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                    e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                    e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                    e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                    e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                    e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                    e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                    e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                    e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                    e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                    e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                    e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                    e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                    e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                    d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                    d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                    d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                    d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                    d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                    d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                    d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                    d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                    d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                    d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                    d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                    d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                    d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                    d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                    d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                    d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                    d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                    d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                    d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                    d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                    d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                    d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                    d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                    d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                    d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                    d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                    d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                    d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                    d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                    d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                    d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                    d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                    d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                    d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                    d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                    d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                    d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                    d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                    d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                    d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                    d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                    d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                    d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                    d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                    d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                    d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                    d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                    d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                    d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                    d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                    d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                    d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                    d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                    d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                    d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                    d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                    d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                    d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                    d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                    d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                    d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                    d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                    d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                    d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                    d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                    d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                    d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                    d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                    d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                    d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                    d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                    d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                    d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                    d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                    d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                    d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                    d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                    d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                    d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                    d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                    d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                    t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                    t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                    t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                    t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                    t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                    t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                    t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                    t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                    t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
                }

                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD21 = derivative3DX(t21,imgWidth,imgHeight,imgDepth,stepSize);
                dervXYD21 = derivative3DY(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD31 = derivative3DX(t31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXZD31 = derivative3DZ(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD12 = derivative3DY(t12,imgWidth,imgHeight,imgDepth,stepSize);
                dervYXD12 = derivative3DX(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD32 = derivative3DY(t32,imgWidth,imgHeight,imgDepth,stepSize);
                dervYZD32 = derivative3DZ(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD13 = derivative3DZ(t13,imgWidth,imgHeight,imgDepth,stepSize);
                dervZXD13 = derivative3DX(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD23 = derivative3DZ(t23,imgWidth,imgHeight,imgDepth,stepSize);
                dervZYD23 = derivative3DY(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == updatedPxlLocs[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStep*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
                }
                delete [] dervX;
                dervX = NULL;
                delete [] dervY;
                dervY = NULL;
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
                delete [] dervXY;
                dervXY = NULL;
                delete [] dervXZ;
                dervXZ = NULL;
                delete [] dervYZ;
                dervYZ = NULL;
                delete [] dervZZD33;
                dervZZD33 = NULL;
                delete [] dervXXD11;
                dervXXD11 = NULL;
                delete [] dervXYD21;
                dervXYD21 = NULL;
                delete [] dervYXD12;
                dervYXD12 = NULL;
                delete [] dervYD21;
                dervYD21 = NULL;
                delete [] dervZD31;
                dervZD31 = NULL;
                delete [] dervYYD22;
                dervYYD22 = NULL;
                delete [] dervXZD31;
                dervXZD31 = NULL;
                delete [] dervYZD32;
                dervYZD32 = NULL;
                delete [] dervZXD13;
                dervZXD13 = NULL;
                delete [] dervZYD23;
                dervZYD23 = NULL;
                delete [] dervXD12;
                dervXD12 = NULL;
                delete [] dervZD32;
                dervZD32 = NULL;
                delete [] dervXD13;
                dervXD13 = NULL;
                delete [] dervYD23;
                dervYD23 = NULL;
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
            }
    //        double l2normErrorOri = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
    //        printf("%lf\n", l2normErrorOri);
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);    // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on z axis after x and y axis.

                //First derivatives of convolved mask
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.

                //Second derivatives
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
                dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xy.
                dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xz.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
                dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. yz.


                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                    double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                    double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];
    //                double en[6][3][3], vn[3][3];

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

    //                    vn[0][0] = 1;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = 0;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 1;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 1;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = 0;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = dervZConv[i]/normi;;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 0;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 0;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = dervXConv[i]/normi;
    //                    vn[0][1] = dervYConv[i]/normi;
    //                    vn[0][2] = dervZConv[i]/normi;

    //                    vn[1][0] = dervYConv[i]/xy_norm_i;
    //                    vn[1][1] = -dervXConv[i]/xy_norm_i;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][2] = -xy_norm_i_square/(xy_norm_i*normi);

    //                    if(v1[2]-vn[0][2] < 0) {
    //                        qDebug() << v1[2]-vn[0][2];
    //                    }

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
                    }

    //                for(int i=0; i<3; i++) {
    //                    if(v1[i]-vn[0][i] < 0) {
    //                        qDebug() << v1[i]-vn[0][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v2[i]-vn[1][i] < 0) {
    //                        qDebug() << v2[i]-vn[1][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v3[i]-vn[2][i] < 0) {
    //                        qDebug() << v3[i]-vn[2][i];
    //                    }
    //                }


    //                m1 = green_diff(normi, 1);
    //                m1 = aubert_diff(norm_i_square, 1);
    //                m1 = li1(normi, 0.2);
    //                m1 = pm_diff(norm_i_square, 0.2);
    //                m1 = gr_diff(norm_i_square, 1);
    //                m1 = pm_diff2(norm_i_square, 0.2);

                    m1 = charbonnier_diff(norm_i_square, 0.561268);
                    m2 = 1;
                    m3 = 1;
                    m4 = sqrt(m1*m2);
                    //m4 = (m1+m2)/2;
                    //m4 = m1;
                    //m4 = 1;

                    m5 = sqrt(m1*m3);
    //                m5 = 0;

                    m6 = 1;
    //                m6 = 0;

//            qDebug() << "Ana";

    //                for(int k=0; k<3; k++) {
    //                    for(int i=0; i<3; i++) {
    //                        for(int j=0; j<3; j++) {
    //                            en[k][i][j] = v;
    //                        }
    //                    }
    //                }



    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e1[i][j] = v1[i]*v1[j];
    //                    }
    //                }
                    e1[0][0] = v1[0]*v1[0];
                    e1[0][1] = v1[0]*v1[1];
                    e1[0][2] = v1[0]*v1[2];
                    e1[1][0] = v1[1]*v1[0];
                    e1[1][1] = v1[1]*v1[1];
                    e1[1][2] = v1[1]*v1[2];
                    e1[2][0] = v1[2]*v1[0];
                    e1[2][1] = v1[2]*v1[1];
                    e1[2][2] = v1[2]*v1[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e2[i][j] = v2[i]*v2[j];
    //                    }
    //                }
                    e2[0][0] = v2[0]*v2[0];
                    e2[0][1] = v2[0]*v2[1];
                    e2[0][2] = v2[0]*v2[2];
                    e2[1][0] = v2[1]*v2[0];
                    e2[1][1] = v2[1]*v2[1];
                    e2[1][2] = v2[1]*v2[2];
                    e2[2][0] = v2[2]*v2[0];
                    e2[2][1] = v2[2]*v2[1];
                    e2[2][2] = v2[2]*v2[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e3[i][j] = v3[i]*v3[j];
    //                    }
    //                }
                    e3[0][0] = v3[0]*v3[0];
                    e3[0][1] = v3[0]*v3[1];
                    e3[0][2] = v3[0]*v3[2];
                    e3[1][0] = v3[1]*v3[0];
                    e3[1][1] = v3[1]*v3[1];
                    e3[1][2] = v3[1]*v3[2];
                    e3[2][0] = v3[2]*v3[0];
                    e3[2][1] = v3[2]*v3[1];
                    e3[2][2] = v3[2]*v3[2];

                    e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                    e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                    e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                    e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                    e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                    e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                    e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                    e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                    e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                    e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                    e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                    e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                    e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                    e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                    e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                    e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                    e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                    e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                    e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                    e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                    e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                    e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                    e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                    e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                    e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                    e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                    e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                    d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                    d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                    d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                    d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                    d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                    d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                    d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                    d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                    d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                    d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                    d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                    d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                    d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                    d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                    d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                    d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                    d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                    d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                    d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                    d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                    d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                    d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                    d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                    d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                    d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                    d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                    d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                    d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                    d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                    d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                    d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                    d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                    d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                    d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                    d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                    d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                    d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                    d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                    d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                    d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                    d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                    d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                    d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                    d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                    d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                    d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                    d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                    d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                    d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                    d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                    d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                    d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                    d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                    d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                    d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                    d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                    d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                    d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                    d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                    d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                    d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                    d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                    d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                    d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                    d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                    d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                    d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                    d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                    d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                    d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                    d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                    d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                    d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                    d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                    d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                    d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                    d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                    d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                    d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                    d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                    d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                    t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                    t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                    t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                    t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                    t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                    t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                    t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                    t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                    t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
                }

                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD21 = derivative3DX(t21,imgWidth,imgHeight,imgDepth,stepSize);
                dervXYD21 = derivative3DY(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD31 = derivative3DX(t31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXZD31 = derivative3DZ(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD12 = derivative3DY(t12,imgWidth,imgHeight,imgDepth,stepSize);
                dervYXD12 = derivative3DX(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD32 = derivative3DY(t32,imgWidth,imgHeight,imgDepth,stepSize);
                dervYZD32 = derivative3DZ(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD13 = derivative3DZ(t13,imgWidth,imgHeight,imgDepth,stepSize);
                dervZXD13 = derivative3DX(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD23 = derivative3DZ(t23,imgWidth,imgHeight,imgDepth,stepSize);
                dervZYD23 = derivative3DY(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == randPxls[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStep*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
                }
                delete [] dervX;
                dervX = NULL;
                delete [] dervY;
                dervY = NULL;
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
                delete [] dervXY;
                dervXY = NULL;
                delete [] dervXZ;
                dervXZ = NULL;
                delete [] dervYZ;
                dervYZ = NULL;
                delete [] dervZZD33;
                dervZZD33 = NULL;
                delete [] dervXXD11;
                dervXXD11 = NULL;
                delete [] dervXYD21;
                dervXYD21 = NULL;
                delete [] dervYXD12;
                dervYXD12 = NULL;
                delete [] dervYD21;
                dervYD21 = NULL;
                delete [] dervZD31;
                dervZD31 = NULL;
                delete [] dervYYD22;
                dervYYD22 = NULL;
                delete [] dervXZD31;
                dervXZD31 = NULL;
                delete [] dervYZD32;
                dervYZD32 = NULL;
                delete [] dervZXD13;
                dervZXD13 = NULL;
                delete [] dervZYD23;
                dervZYD23 = NULL;
                delete [] dervXD12;
                dervXD12 = NULL;
                delete [] dervZD32;
                dervZD32 = NULL;
                delete [] dervXD13;
                dervXD13 = NULL;
                delete [] dervYD23;
                dervYD23 = NULL;
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
            }
    //        double l2normErrorOri = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
    //        printf("%lf\n", l2normErrorOri);
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    //Make nonnegative voxel values
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(scatImageArr[i] < 0) {
            scatImageArr[i] = 0;
        }
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] d1111;
    d1111 = NULL;
    delete [] d1112;
    d1112 = NULL;
    delete [] d1113;
    d1113 = NULL;
    delete [] d1121;
    d1121 = NULL;
    delete [] d1122;
    d1122 = NULL;
    delete [] d1123;
    d1123 = NULL;
    delete [] d1131;
    d1131 = NULL;
    delete [] d1132;
    d1132 = NULL;
    delete [] d1133;
    d1133 = NULL;

    delete [] d1211;
    d1211 = NULL;
    delete [] d1212;
    d1212 = NULL;
    delete [] d1213;
    d1213 = NULL;
    delete [] d1221;
    d1221 = NULL;
    delete [] d1222;
    d1222 = NULL;
    delete [] d1223;
    d1223 = NULL;
    delete [] d1231;
    d1231 = NULL;
    delete [] d1232;
    d1232 = NULL;
    delete [] d1233;
    d1233 = NULL;

    delete [] d1311;
    d1311 = NULL;
    delete [] d1312;
    d1312 = NULL;
    delete [] d1313;
    d1313 = NULL;
    delete [] d1321;
    d1321 = NULL;
    delete [] d1322;
    d1322 = NULL;
    delete [] d1323;
    d1323 = NULL;
    delete [] d1331;
    d1331 = NULL;
    delete [] d1332;
    d1332 = NULL;
    delete [] d1333;
    d1333 = NULL;


    delete [] d2111;
    d2111 = NULL;
    delete [] d2112;
    d2112 = NULL;
    delete [] d2113;
    d2113 = NULL;
    delete [] d2121;
    d2121 = NULL;
    delete [] d2122;
    d2122 = NULL;
    delete [] d2123;
    d2123 = NULL;
    delete [] d2131;
    d2131 = NULL;
    delete [] d2132;
    d2132 = NULL;
    delete [] d2133;
    d2133 = NULL;

    delete [] d2211;
    d2211 = NULL;
    delete [] d2212;
    d2212 = NULL;
    delete [] d2213;
    d2213 = NULL;
    delete [] d2221;
    d2221 = NULL;
    delete [] d2222;
    d2222 = NULL;
    delete [] d2223;
    d2223 = NULL;
    delete [] d2231;
    d2231 = NULL;
    delete [] d2232;
    d2232 = NULL;
    delete [] d2233;
    d2233 = NULL;

    delete [] d2311;
    d2311 = NULL;
    delete [] d2312;
    d2312 = NULL;
    delete [] d2313;
    d2313 = NULL;
    delete [] d2321;
    d2321 = NULL;
    delete [] d2322;
    d2322 = NULL;
    delete [] d2323;
    d2323 = NULL;
    delete [] d2331;
    d2331 = NULL;
    delete [] d2332;
    d2332 = NULL;
    delete [] d2333;
    d2333 = NULL;


    delete [] d3111;
    d3111 = NULL;
    delete [] d3112;
    d3112 = NULL;
    delete [] d3113;
    d3113 = NULL;
    delete [] d3121;
    d3121 = NULL;
    delete [] d3122;
    d3122 = NULL;
    delete [] d3123;
    d3123 = NULL;
    delete [] d3131;
    d3131 = NULL;
    delete [] d3132;
    d3132 = NULL;
    delete [] d3133;
    d3133 = NULL;

    delete [] d3211;
    d3211 = NULL;
    delete [] d3212;
    d3212 = NULL;
    delete [] d3213;
    d3213 = NULL;
    delete [] d3221;
    d3221 = NULL;
    delete [] d3222;
    d3222 = NULL;
    delete [] d3223;
    d3223 = NULL;
    delete [] d3231;
    d3231 = NULL;
    delete [] d3232;
    d3232 = NULL;
    delete [] d3233;
    d3233 = NULL;

    delete [] d3311;
    d3311 = NULL;
    delete [] d3312;
    d3312 = NULL;
    delete [] d3313;
    d3313 = NULL;
    delete [] d3321;
    d3321 = NULL;
    delete [] d3322;
    d3322 = NULL;
    delete [] d3323;
    d3323 = NULL;
    delete [] d3331;
    d3331 = NULL;
    delete [] d3332;
    d3332 = NULL;
    delete [] d3333;
    d3333 = NULL;

    delete [] t11;
    t11 = NULL;
    delete [] t12;
    t12 = NULL;
    delete [] t13;
    t13 = NULL;
    delete [] t21;
    t21 = NULL;
    delete [] t22;
    t22 = NULL;
    delete [] t23;
    t23 = NULL;
    delete [] t31;
    t31 = NULL;
    delete [] t32;
    t32 = NULL;
    delete [] t33;
    t33 = NULL;

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}


void st_eed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, int imgVolume, float gridSpcX, float gridSpcY, float gridSpcZ, float gridSpcT, bool zeros2mask)
{
    clock_t begin = clock();

    int kernelSize = 3, randArrTraceIndex, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror, aaerror, l2normError, alpha;
    double* *outConvX;
    double* *outConvXY;
    double* *outConvXYZ;
    double* *outConvXYZT;
    double* *dervY;
    double* *dervX;
    double* *dervZ;
    double* *dervT;
    double* *dervXConv;
    double* *dervYConv;
    double* *dervZConv;
    double* *dervTConv;
    double* *dervYD21;
    double* *dervXD12;
    double* *dervXD13;
    double* *dervXD14;
    double* *dervYD23;
    double* *dervYD24;
    double* *dervZD31;
    double* *dervZD32;
    double* *dervZD34;
    double* *dervTD41;
    double* *dervTD42;
    double* *dervTD43;
    double* *dervForX;
    double* *dervBacX;
    double* *dervForY;
    double* *dervBacY;
    double* *sumForX;
    double* *sumBacX;
    double* *sumForY;
    double* *sumBacY;
    double* *sumForZ;
    double* *sumBacZ;
    double* *sumForT;
    double* *sumBacT;
    double* *dervForZ;
    double* *dervForT;
    double* *dervBacZ;
    double* *dervBacT;
    double* *tempImgArrayPrev;
    double* *tempImgArrayCurr;
    double v1[4], v2[4], v3[4], v4[4], ones[4];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

    //Memory allocation for diffusion tensor entries dij.
    double** d11 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d11[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d12 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d12[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d13 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d13[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d14 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d14[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    double** d21 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d21[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d22 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d22[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d23 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d23[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d24 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d24[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    double** d31 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d31[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d32 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d32[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d33 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d33[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d34 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d34[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    double** d41 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d41[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d42 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d42[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d43 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d43[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double** d44 = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        d44[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];



    gridSpcX = gridSpcY = gridSpcZ = gridSpcT =  1;

    qDebug() << gridSpcX;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        tempImgArrayCurr[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double*[sizeof(double) * imgVolume];
    for(int i = 0; i<imgVolume; i++)
        tempImgArrayPrev[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0, locsVectLenn = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";

        for(int i=0; i<imgWidth*imgHeight*imgDepth*imgVolume; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth*imgVolume; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }

    l2normError = l2Norm_4D(imageArr,scatImageArr,imgVolume, imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    //Contrast Paremeter estimation
    double quantile = 90.0;
    double* estContPar;
//    estContPar = quantlCriterPM(quantile, imageArr,imgWidth,imgHeight,imgDepth,sigma,kernelSize);

    if(zeros2mask) {
            for (int volume=0; volume<imgVolume;volume++)
                estContPar[volume] = quantlCriterPM4Inpainting(quantile, imageArr[volume], updatedPxlLocs, imgWidth,imgHeight,imgDepth,sigma,kernelSize);
        } else {
                for (int volume=0; volume<imgVolume;volume++)
                    estContPar[volume] = quantlCriterPM4Inpainting(quantile, imageArr[volume], randPxls, imgWidth, imgHeight, imgDepth, sigma, kernelSize);
        }
    for (int volume=0; volume<imgVolume;volume++) {
        estContPar[volume] = histoCriterPM4Inpainting(quantile, imageArr[volume], randPxls, imgWidth, imgHeight, imgDepth);
        estContPar[volume] = estContPar[volume]/25;
        qDebug() << "Estimated Contrast Parameter: for volume " << volume << "    "<< estContPar[volume];
    }


    double firstSteps = 3, itrCount = 0;
    if(zeros2mask) {
//        N = numSteps;
        while(l2normError > tol) {
            if(itrCount > 3) {
                N = numSteps;
            } else {
                N = firstSteps;
            }
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX_of_4D(scatImageArr, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                outConvXY = convolution3DY_of_4D(outConvX, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ_of_4D(outConvXY, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                outConvXYZT = convolution3DT_of_4D(outConvXY, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on t axis after x and y and z axis.
                //First convolved derivatives
                dervXConv = derivative4DX(outConvXYZ,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative4DY(outConvXYZ,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative4DZ(outConvXYZ,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);                         // Derivative of convolved image w.r.t. z.
                dervTConv = derivative4DT(outConvXYZ,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);                         // Derivative of convolved image w.r.t. t.

                //First derivatives
                dervX = derivative4DX(scatImageArr,imgWidth,imgHeight,imgDepth,imgVolume, gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative4DY(scatImageArr,imgWidth,imgHeight,imgDepth,imgVolume, gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative4DZ(scatImageArr,imgWidth,imgHeight,imgDepth,imgVolume, gridSpcZ);                           // Derivative w.r.t. z.
                dervT = derivative4DT(scatImageArr,imgWidth,imgHeight,imgDepth,imgVolume, gridSpcZ);                           // Derivative w.r.t. t.

                dervForX = dervForw4DX(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcX);
                dervForY = dervForw4DY(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcY);
                dervForZ = dervForw4DZ(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);
                dervForT = dervForw4DT(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);
                dervBacX = dervBack4DX(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcX);
                dervBacY = dervBack4DY(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcY);
                dervBacZ = dervBack4DZ(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);
                dervBacT = dervBack4DT(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);

                for (int volume=0; volume< imgVolume; volume++){
                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(n == 0) {
                            tempImgArrayPrev[volume][i] = scatImageArr[volume][i];
                        }
                        tempImgArrayCurr[volume][i] = scatImageArr[volume][i];

                        //Define eigenvectors and entries for the diffusion tensor.
                        double norm_i_square = dervXConv[volume][i]*dervXConv[volume][i] + dervYConv[volume][i]*dervYConv[volume][i] + dervZConv[volume][i]*dervZConv[volume][i]
                                + dervTConv[volume][i]*dervTConv[volume][i],
                                normi = sqrt(norm_i_square), d1, d2, d3, d4;
                        double xy_norm_i_square = dervXConv[volume][i]*dervXConv[volume][i] + dervYConv[volume][i]*dervYConv[volume][i], xy_norm_i = sqrt(xy_norm_i_square);
                        double yz_norm_i_square = dervYConv[volume][i]*dervYConv[volume][i] + dervZConv[volume][i]*dervZConv[volume][i], yz_norm_i = sqrt(yz_norm_i_square);
                        double xyz_norm_i_square = dervXConv[volume][i]*dervXConv[volume][i] + dervYConv[volume][i]*dervYConv[volume][i] + dervZConv[volume][i]*dervZConv[volume][i]
                                , xyz_norm_i = sqrt(xyz_norm_i_square);

                        if(norm_i_square == 0) {
                            v1[0] = 1;
                            v1[1] = 0;
                            v1[2] = 0;
                            v1[3] = 0;

                            v2[0] = 0;
                            v2[1] = 1;
                            v2[2] = 0;
                            v2[3] = 0;

                            v3[0] = 0;
                            v3[1] = 0;
                            v3[2] = 1;
                            v3[3] = 0;

                            v4[0] = 0;
                            v4[1] = 0;
                            v4[2] = 0;
                            v4[3] = 1;

                        } else if(xyz_norm_i_square == 0) {
                            v1[0] = 0;
                            v1[1] = 0;
                            v1[2] = 0;
                            v1[3] = dervTConv[volume][i]/normi;

                            v2[0] = 0;
                            v2[1] = 0;
                            v2[2] = 0;
                            v2[3] = 0;

                            v3[0] = 0;
                            v3[1] = 0;
                            v3[2] = 0;
                            v3[3] = 0;

                            v4[0] = 0;
                            v4[1] = 0;
                            v4[2] = 0;
                            v4[3] = 0;
                        } else {
                            v1[0] = dervXConv[volume][i]/normi;
                            v1[1] = dervYConv[volume][i]/normi;
                            v1[2] = dervZConv[volume][i]/normi;
                            v1[3] = dervTConv[volume][i]/normi;

                            v2[0] = dervYConv[volume][i]/xy_norm_i;
                            v2[1] = -dervXConv[volume][i]/xy_norm_i;
                            v2[2] = 0;
                            v2[3] = 0;

                            v3[0] = 0;
                            v3[1] = - dervZConv[volume][i]/yz_norm_i;
                            v3[2] = dervYConv[volume][i]/yz_norm_i;
                            v3[3] = 0;

                            v4[0] = (dervXConv[volume][i]*dervTConv[volume][i])/(xyz_norm_i*normi);
                            v4[1] = (dervYConv[volume][i]*dervTConv[volume][i])/(xyz_norm_i*normi);
                            v4[2] = (dervZConv[volume][i]*dervTConv[volume][i])/(xyz_norm_i*normi);
                            v4[3] = -xyz_norm_i_square/(xyz_norm_i*normi);
                        }

//                        double contPar = 0.65364;
//                        double contPar = 0.1;
                        d1 = charbonnier_diff(norm_i_square, estContPar[volume]);
                        d2 = 1;
                        d3 = 1;
                        d4 = 1;

    //                    qDebug() << "Zeros Included - EED";

                        //Diffusion tensor entries definition with eigenvalues d1, d2, d3, d4 and orthgonal orthonormal eigenvectors v1, v2, v3 and v4.
                        d11[volume][i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0] + d3*v3[0]*v3[0] + d4*v4[0]*v4[0];
                        d12[volume][i] = (d1*v1[0]*v1[1] + d2*v2[0]*v2[1] + d3*v3[0]*v3[1] + d4*v4[0]*v4[1])*dervY[volume][i];
                        d13[volume][i] = (d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2] + d4*v4[0]*v4[2])*dervZ[volume][i];
                        d14[volume][i] = (d1*v1[0]*v1[3] + d2*v2[0]*v2[3] + d3*v3[0]*v3[3] + d4*v4[0]*v4[3])*dervT[volume][i];

                        d21[volume][i] = (d1*v1[1]*v1[0] + d2*v2[1]*v2[0] + d3*v3[1]*v3[0] + d4*v4[1]*v4[0])*dervX[volume][i];
                        d22[volume][i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1] + d3*v3[1]*v3[1] + d4*v4[1]*v4[1];
                        d23[volume][i] = (d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2] + d4*v4[1]*v4[2])*dervZ[volume][i];
                        d24[volume][i] = (d1*v1[1]*v1[3] + d2*v2[1]*v2[3] + d3*v3[1]*v3[3] + d4*v4[1]*v4[3])*dervT[volume][i];

                        d31[volume][i] = (d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2] + d4*v4[0]*v4[2] )*dervX[volume][i];
                        d32[volume][i] = (d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2] + d4*v4[1]*v4[2] )*dervY[volume][i];
                        d33[volume][i] = d1*v1[2]*v1[2] + d2*v2[2]*v2[2] + d3*v3[2]*v3[2] + d4*v4[2]*v4[2];
                        d34[volume][i] = (d1*v1[3]*v1[2] + d2*v2[3]*v2[2] + d3*v3[3]*v3[2] + d4*v4[3]*v4[2] )*dervT[volume][i];

                        d41[volume][i] = (d1*v1[0]*v1[3] + d2*v2[0]*v2[3] + d3*v3[0]*v3[3] + d4*v4[0]*v4[3] )*dervX[volume][i];
                        d42[volume][i] = (d1*v1[1]*v1[3] + d2*v2[1]*v2[3] + d3*v3[1]*v3[3] + d4*v4[1]*v4[3] )*dervY[volume][i];
                        d43[volume][i] = (d1*v1[2]*v1[3] + d2*v2[2]*v2[3] + d3*v3[2]*v3[3] + d4*v4[2]*v4[3] )*dervZ[volume][i];
                        d44[volume][i] = (d1*v1[3]*v1[3] + d2*v2[3]*v2[3] + d3*v3[3]*v3[3] + d4*v4[3]*v4[3] );
                    }
                }
                sumForX = sumForw4DX(d11,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForY = sumForw4DY(d22,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForZ = sumForw4DZ(d33,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForT = sumForw4DT(d33,imgWidth,imgHeight,imgDepth, imgVolume);
                sumBacX = sumBack4DX(d11,imgWidth,imgHeight,imgDepth, imgVolume);
                sumBacY = sumBack4DY(d22,imgWidth,imgHeight,imgDepth, imgVolume);
                sumBacZ = sumBack4DZ(d33,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForT = sumBack4DT(d33,imgWidth,imgHeight,imgDepth, imgVolume);

                dervXD12 = derivative4DX(d12,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);
                dervXD13 = derivative4DX(d13,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);
                dervXD14 = derivative4DX(d14,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);
                dervYD21 = derivative4DY(d21,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);
                dervYD23 = derivative4DY(d23,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);
                dervYD24 = derivative4DY(d24,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);
                dervZD31 = derivative4DZ(d31,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervZD32 = derivative4DZ(d32,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervZD34 = derivative4DZ(d34,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervTD41 = derivative4DZ(d41,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervTD42 = derivative4DZ(d42,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervTD43 = derivative4DZ(d43,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);


                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for (int volume = 0; volume < imgVolume; volume++){
                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(i == updatedPxlLocs[randArrTraceIndex]) {
                            randArrTraceIndex++;
                            continue;
                        }
                        scatImageArr[volume][i] = alpha*(scatImageArr[volume][i] + timeStepSize*((sumForX[volume][i]*dervForX[volume][i] - sumBacX[volume][i]*dervBacX[volume][i] +
                                                      sumForY[volume][i]*dervForY[volume][i] - sumBacY[volume][i]*dervBacY[volume][i] +
                                                      sumForZ[volume][i]*dervForZ[volume][i] - sumBacZ[volume][i]*dervBacZ[volume][i] +
                                                      sumForT[volume][i]*dervForZ[volume][i] - sumBacT[volume][i]*dervBacZ[volume][i])/2 +
                                                      dervXD12[volume][i] + dervXD13[volume][i] + dervXD14[volume][i] +
                                                      dervYD21[volume][i] + dervYD23[volume][i] + dervYD24[volume][i] +
                                                      dervZD31[volume][i] + dervZD32[volume][i] + dervZD34[volume][i] +
                                                      dervTD41[volume][i] + dervTD42[volume][i] + dervTD43[volume][i])) +
                                                      (1-alpha)*tempImgArrayPrev[volume][i];
                        tempImgArrayPrev[volume][i] = tempImgArrayCurr[volume][i];
                    }
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
            }
            l2normError = l2Norm_4D(tempImgArrayCurr,scatImageArr, imgVolume, imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";

            itrCount++;
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX_of_4D(scatImageArr, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                outConvXY = convolution3DY_of_4D(outConvX, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ_of_4D(outConvXY, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                outConvXYZT = convolution3DT_of_4D(outConvXY, imgVolume, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on t axis after x and y and z axis.

                //First convolved derivatives
                dervXConv = derivative4DX(outConvXYZT,imgWidth,imgHeight,imgDepth,imgVolume,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative4DY(outConvXYZT,imgWidth,imgHeight,imgDepth,imgVolume,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative4DZ(outConvXYZT,imgWidth,imgHeight,imgDepth,imgVolume,gridSpcZ);                         // Derivative of convolved image w.r.t. z.
                dervTConv = derivative4DT(outConvXYZT,imgWidth,imgHeight,imgDepth,imgVolume,gridSpcZ);                         // Derivative of convolved image w.r.t. t.

    //            //Convolved First derivatives
    //            dervForXConv = dervForw3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
    //            dervForYConv = dervForw3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
    //            dervForZConv = dervForw3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);
    //            dervBacXConv = dervBack3DX(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcX);
    //            dervBacYConv = dervBack3DY(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcY);
    //            dervBacZConv = dervBack3DZ(outConvXYZ, imgWidth, imgHeight, imgDepth, gridSpcZ);

                //First derivatives
                dervX = derivative4DX(scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative4DY(scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative4DZ(scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);                           // Derivative w.r.t. z.
                dervT = derivative4DT(scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);                           // Derivative w.r.t. t.

                dervForX = dervForw4DX(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcX);
                dervForY = dervForw4DY(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcY);
                dervForZ = dervForw4DZ(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);
                dervForT = dervForw4DT(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);
                dervBacX = dervBack4DX(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcX);
                dervBacY = dervBack4DY(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcY);
                dervBacZ = dervBack4DZ(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);
                dervBacT = dervBack4DT(scatImageArr, imgWidth, imgHeight, imgDepth, imgVolume, gridSpcZ);

    // //             Adaptive contrast parameter selection******************************
    //            double tv_norm = 0, contPar;
    //            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
    //                double norm_i_square, normi;

    //                norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i];
    //                normi = sqrt(norm_i_square);
    //                tv_norm = tv_norm + normi;
    //            }
    //            tv_norm = tv_norm/(imgWidth*imgHeight*imgDepth);
    //            contPar = 0.002*tv_norm;
    // //            qDebug() << tv_norm << contPar;
    // //            ********************************************************************
                for (int volume = 0 ; volume < imgVolume; volume++) {
                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(n == 0) {
                            tempImgArrayPrev[i] = scatImageArr[i];
                        }
                        tempImgArrayCurr[i] = scatImageArr[i];

                        //Define eigenvectors and entries for the diffusion tensor.
                        double norm_i_square = dervXConv[volume][i]*dervXConv[volume][i] + dervYConv[volume][i]*dervYConv[volume][i] +
                                dervZConv[volume][i]*dervZConv[volume][i] + dervTConv[volume][i]*dervTConv[volume][i],
                                normi = sqrt(norm_i_square), d1, d2, d3, d4;
                        double xy_norm_i_square = dervXConv[volume][i]*dervXConv[volume][i] + dervYConv[volume][i]*dervYConv[volume][i], xy_norm_i = sqrt(xy_norm_i_square);
                        double yz_norm_i_square = dervYConv[volume][i]*dervYConv[volume][i] + dervZConv[volume][i]*dervZConv[volume][i], yz_norm_i = sqrt(yz_norm_i_square);
                        double xyz_norm_i_square = dervXConv[volume][i]*dervXConv[volume][i] + dervYConv[volume][i]*dervYConv[volume][i] + dervZConv[volume][i]*dervZConv[volume][i]
                                , xyz_norm_i = sqrt(xyz_norm_i_square);

                        if(norm_i_square == 0) {
                            v1[0] = 1;
                            v1[1] = 0;
                            v1[2] = 0;
                            v1[3] = 0;

                            v2[0] = 0;
                            v2[1] = 1;
                            v2[2] = 0;
                            v2[3] = 0;

                            v3[0] = 0;
                            v3[1] = 0;
                            v3[2] = 1;
                            v3[3] = 0;

                            v4[0] = 0;
                            v4[1] = 0;
                            v4[2] = 0;
                            v4[3] = 1;

                        } else if(xyz_norm_i_square == 0) {
                            v1[0] = 0;
                            v1[1] = 0;
                            v1[2] = 0;
                            v1[3] = dervTConv[volume][i]/normi;

                            v2[0] = 0;
                            v2[1] = 0;
                            v2[2] = 0;
                            v2[3] = 0;

                            v3[0] = 0;
                            v3[1] = 0;
                            v3[2] = 0;
                            v3[3] = 0;

                            v4[0] = 0;
                            v4[1] = 0;
                            v4[2] = 0;
                            v4[3] = 0;
                        } else {
                            v1[0] = dervXConv[volume][i]/normi;
                            v1[1] = dervYConv[volume][i]/normi;
                            v1[2] = dervZConv[volume][i]/normi;
                            v1[3] = dervTConv[volume][i]/normi;

                            v2[0] = dervYConv[volume][i]/xy_norm_i;
                            v2[1] = -dervXConv[volume][i]/xy_norm_i;
                            v2[2] = 0;
                            v2[3] = 0;

                            v3[0] = 0;
                            v3[1] = - dervZConv[volume][i]/yz_norm_i;
                            v3[2] = dervYConv[volume][i]/yz_norm_i;
                            v3[3] = 0;

                            v4[0] = (dervXConv[volume][i]*dervTConv[volume][i])/(xyz_norm_i*normi);
                            v4[1] = (dervYConv[volume][i]*dervTConv[volume][i])/(xyz_norm_i*normi);
                            v4[2] = (dervZConv[volume][i]*dervTConv[volume][i])/(xyz_norm_i*normi);
                            v4[3] = -xyz_norm_i_square/(xyz_norm_i*normi);
                        }

                        d1 = charbonnier_diff(norm_i_square, estContPar[volume]);
                        d2 = 1;
                        d3 = 1;
                        d4 = 1;

                        //Diffusion tensor entries definition with eigenvalues d1, d2, d3 and orthgonal orthonormal eigenvectors v1, v2 and v3.
                        d11[volume][i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0] + d3*v3[0]*v3[0] + d4*v4[0]*v4[0];
                        d12[volume][i] = (d1*v1[0]*v1[1] + d2*v2[0]*v2[1] + d3*v3[0]*v3[1] + d4*v4[0]*v4[1])*dervY[volume][i];
                        d13[volume][i] = (d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2] + d4*v4[0]*v4[2])*dervZ[volume][i];
                        d14[volume][i] = (d1*v1[0]*v1[3] + d2*v2[0]*v2[3] + d3*v3[0]*v3[3] + d4*v4[0]*v4[3])*dervT[volume][i];

                        d21[volume][i] = (d1*v1[1]*v1[0] + d2*v2[1]*v2[0] + d3*v3[1]*v3[0] + d4*v4[1]*v4[0])*dervX[volume][i];
                        d22[volume][i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1] + d3*v3[1]*v3[1] + d4*v4[1]*v4[1];
                        d23[volume][i] = (d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2] + d4*v4[1]*v4[2])*dervZ[volume][i];
                        d24[volume][i] = (d1*v1[1]*v1[3] + d2*v2[1]*v2[3] + d3*v3[1]*v3[3] + d4*v4[1]*v4[3])*dervT[volume][i];

                        d31[volume][i] = (d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2] + d4*v4[0]*v4[2] )*dervX[volume][i];
                        d32[volume][i] = (d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2] + d4*v4[1]*v4[2] )*dervY[volume][i];
                        d33[volume][i] = d1*v1[2]*v1[2] + d2*v2[2]*v2[2] + d3*v3[2]*v3[2] + d4*v4[2]*v4[2];
                        d34[volume][i] = (d1*v1[3]*v1[2] + d2*v2[3]*v2[2] + d3*v3[3]*v3[2] + d4*v4[3]*v4[2] )*dervT[volume][i];

                        d41[volume][i] = (d1*v1[0]*v1[3] + d2*v2[0]*v2[3] + d3*v3[0]*v3[3] + d4*v4[0]*v4[3] )*dervX[volume][i];
                        d42[volume][i] = (d1*v1[1]*v1[3] + d2*v2[1]*v2[3] + d3*v3[1]*v3[3] + d4*v4[1]*v4[3] )*dervY[volume][i];
                        d43[volume][i] = (d1*v1[2]*v1[3] + d2*v2[2]*v2[3] + d3*v3[2]*v3[3] + d4*v4[2]*v4[3] )*dervZ[volume][i];
                        d44[volume][i] = (d1*v1[3]*v1[3] + d2*v2[3]*v2[3] + d3*v3[3]*v3[3] + d4*v4[3]*v4[3] );
                    }
                }
                sumForX = sumForw4DX(d11,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForY = sumForw4DY(d22,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForZ = sumForw4DZ(d33,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForT = sumForw4DT(d33,imgWidth,imgHeight,imgDepth, imgVolume);
                sumBacX = sumBack4DX(d11,imgWidth,imgHeight,imgDepth, imgVolume);
                sumBacY = sumBack4DY(d22,imgWidth,imgHeight,imgDepth, imgVolume);
                sumBacZ = sumBack4DZ(d33,imgWidth,imgHeight,imgDepth, imgVolume);
                sumForT = sumBack4DT(d33,imgWidth,imgHeight,imgDepth, imgVolume);

                dervXD12 = derivative4DX(d12,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);
                dervXD13 = derivative4DX(d13,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);
                dervXD14 = derivative4DX(d14,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcX);
                dervYD21 = derivative4DY(d21,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);
                dervYD23 = derivative4DY(d23,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);
                dervYD24 = derivative4DY(d24,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcY);
                dervZD31 = derivative4DZ(d31,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervZD32 = derivative4DZ(d32,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervZD34 = derivative4DZ(d34,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervTD41 = derivative4DZ(d41,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervTD42 = derivative4DZ(d42,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);
                dervTD43 = derivative4DZ(d43,imgWidth,imgHeight,imgDepth, imgVolume, gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);
                for (int volume = 0; volume < imgVolume; volume++){
                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(i == updatedPxlLocs[randArrTraceIndex]) {
                            randArrTraceIndex++;
                            continue;
                        }
                        scatImageArr[volume][i] = alpha*(scatImageArr[volume][i] + timeStepSize*((sumForX[volume][i]*dervForX[volume][i] - sumBacX[volume][i]*dervBacX[volume][i] +
                                                      sumForY[volume][i]*dervForY[volume][i] - sumBacY[volume][i]*dervBacY[volume][i] +
                                                      sumForZ[volume][i]*dervForZ[volume][i] - sumBacZ[volume][i]*dervBacZ[volume][i] +
                                                      sumForT[volume][i]*dervForZ[volume][i] - sumBacT[volume][i]*dervBacZ[volume][i])/2 +
                                                      dervXD12[volume][i] + dervXD13[volume][i] + dervXD14[volume][i] +
                                                      dervYD21[volume][i] + dervYD23[volume][i] + dervYD24[volume][i] +
                                                      dervZD31[volume][i] + dervZD32[volume][i] + dervZD34[volume][i] +
                                                      dervTD41[volume][i] + dervTD42[volume][i] + dervTD43[volume][i])) +
                                                      (1-alpha)*tempImgArrayPrev[volume][i];
                        tempImgArrayPrev[volume][i] = tempImgArrayCurr[volume][i];
                    }
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

    //            delete [] dervForXConv;
    //            dervForXConv = NULL;
    //            delete [] dervForYConv;
    //            dervForYConv = NULL;
    //            delete [] dervForZConv;
    //            dervForZConv = NULL;
    //            delete [] dervBacXConv;
    //            dervBacXConv = NULL;
    //            delete [] dervBacYConv;
    //            dervBacYConv = NULL;
    //            delete [] dervBacZConv;
    //            dervBacZConv = NULL;
            }
            l2normError = l2Norm_4D(tempImgArrayCurr,scatImageArr, imgVolume, imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    //Make voxel values nonnegative
    for (int v = 0; v < imgVolume; v++) {
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(scatImageArr[v][i] < 0) {
                scatImageArr[v][i] = 0;
            }
        }
    }

//    //Residual and its entropy calculated here*************************************************************************
//    //uint16 max residual
//    double* resiImg = new double[sizeof(double) * (imgWidth*imgHeight*imgDepth-locsVectLenn)];
//    double maxVal = *max_element(imageArr , imageArr + imgWidth*imgHeight*imgDepth);
//    randArrTraceIndex = 0;
//    int indCnt = 0;
//    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//        if(i == updatedPxlLocs[randArrTraceIndex]) {
//            randArrTraceIndex++;
//            continue;
//        }
//        resiImg[indCnt] = round(scatImageArr[i] - imageArr[i]);
//        if(resiImg[indCnt] < 0) {
//            resiImg[indCnt] = resiImg[indCnt] + maxVal;
//        }
//        indCnt++;
//    }
//    //Entropy of residual
//    double entrpResi = entropy(resiImg, (imgWidth*imgHeight*imgDepth-locsVectLenn));
//    qDebug() << "Residual Entropy: " << entrpResi << "Max: " << maxVal;
//    delete [] resiImg;
//    resiImg = NULL;
//    //End***************************************************************************************************************

    mserror = mse_4D(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth, imgVolume);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae_4D(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth, imgVolume);
    qDebug() << "AAE error: " << aaerror << "\n";


    if(zeros2mask) {
        mserror = mse_mask_dilatied_region4D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume,updatedPxlLocs);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region4D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume,updatedPxlLocs);
        qDebug() << "AAE error: " << aaerror << "\n";
    } else {
        mserror = mse_mask_dilatied_region4D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume,randPxls);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae_mask_dilatied_region4D(imageArr,scatImageArr,imgWidth,imgHeight,imgDepth, imgVolume,randPxls);
        qDebug() << "AAE error: " << aaerror << "\n";
    }


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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);

}


//***************************************************************************************************************************************************************************************************************************//
//~ Fourth Order Edge Enhancing Anisotropic Diffusion based Inpaiting by Fast Semi-Iterative Scheme with fourth order diffusion tensor for 4D volumes
void st_foeed_4d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, int imgVolume, bool zeros2mask)
{
    clock_t begin = clock();

    int kernelSize = 3, N, randArrTraceIndex;
    float stepSize = 1;
    double sigma = 1.0, alpha, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervXX, *dervYY, *dervZZ, *dervXY, *dervYZ, *dervX, *dervY, *dervXZ, *dervXConv, *dervYConv, *dervZConv, *t11, *t12, *t13, *t21, *t22, *t23, *t31, *t32, *t33;
    double *d1111, *d1112, *d1113, *d1121, *d1122, *d1123, *d1131, *d1132, *d1133, *d1211, *d1212, *d1213, *d1221, *d1222, *d1223, *d1231, *d1232, *d1233, *d1311, *d1312, *d1313, *d1321, *d1322, *d1323, *d1331, *d1332, *d1333;
    double *d2111, *d2112, *d2113, *d2121, *d2122, *d2123, *d2131, *d2132, *d2133, *d2211, *d2212, *d2213, *d2221, *d2222, *d2223, *d2231, *d2232, *d2233, *d2311, *d2312, *d2313, *d2321, *d2322, *d2323, *d2331, *d2332, *d2333;
    double *d3111, *d3112, *d3113, *d3121, *d3122, *d3123, *d3131, *d3132, *d3133, *d3211, *d3212, *d3213, *d3221, *d3222, *d3223, *d3231, *d3232, *d3233, *d3311, *d3312, *d3313, *d3321, *d3322, *d3323, *d3331, *d3332, *d3333;
    double *dervXXD11, *dervYYD22, *dervZZD33, *dervYD21, *dervXYD21, *dervXD12, *dervYXD12, *dervZD31, *dervXZD31, *dervXD13, *dervZXD13, *dervYD23, *dervZYD23, *dervZD32, *dervYZD32;
    double *tempImgArrayCurr, *tempImgArrayPrev;

    qDebug() << "4D FOEED: " << "Kernel Size: " << kernelSize;

    //************************************************************
    //Memory allocation for 4th order diffusion tensor entries.
    d1111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    //************************************************************
//    double* te[3][3];
//    for(int i=0; i<3; i++) {
//        for(int j=0; j<3; j++) {
//            te[i][j] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
//        }
//    }

    t11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    //Memory allocation for temporary image array.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    int* updatedPxlLocs;
    vector<int> locsVect;
    int numZeros = 0, arrTrcIndx = 0;
    if(zeros2mask) {
        //Add zero voxels to Mask
        qDebug() << "Zeros considered as known pixels!";
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(imageArr[i] == 0) {
//                scatImageArr[i] = imageArr[i];
                locsVect.push_back(i);
                numZeros++;
            }
        }
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[arrTrcIndx]) {
                locsVect.push_back(i);
                arrTrcIndx++;
            }
        }
        sort(locsVect.begin(), locsVect.end());
        locsVect.erase(unique(locsVect.begin(), locsVect.end()), locsVect.end() );

        int locsVectLenn = locsVect.size();
        updatedPxlLocs = new int[sizeof(int) * locsVectLenn];
        for(int i=0; i<locsVectLenn; i++) {
            updatedPxlLocs[i] = locsVect[i];
        }
    }

//    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
//    qDebug() << "MSE error: " << mserror << "\n";
//    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
//    printf("Error: %lf\n", aaerror);
//    qDebug() << "Here";
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    double firstSteps = 3, itrCount = 0;
    if(zeros2mask) {
        N = numSteps;
        while(l2normError > tol) {
            if(itrCount > 4) {
                N = numSteps;
            } else {
                N = firstSteps;
            }
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);    // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on z axis after x and y axis.

                //First derivatives of convolved mask
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.

                //Second derivatives
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
                dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xy.
                dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xz.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
                dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. yz.


                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                    double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                    double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];
    //                double en[6][3][3], vn[3][3];

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

    //                    vn[0][0] = 1;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = 0;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 1;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 1;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = 0;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = dervZConv[i]/normi;;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 0;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 0;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = dervXConv[i]/normi;
    //                    vn[0][1] = dervYConv[i]/normi;
    //                    vn[0][2] = dervZConv[i]/normi;

    //                    vn[1][0] = dervYConv[i]/xy_norm_i;
    //                    vn[1][1] = -dervXConv[i]/xy_norm_i;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][2] = -xy_norm_i_square/(xy_norm_i*normi);

    //                    if(v1[2]-vn[0][2] < 0) {
    //                        qDebug() << v1[2]-vn[0][2];
    //                    }

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
                    }

    //                for(int i=0; i<3; i++) {
    //                    if(v1[i]-vn[0][i] < 0) {
    //                        qDebug() << v1[i]-vn[0][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v2[i]-vn[1][i] < 0) {
    //                        qDebug() << v2[i]-vn[1][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v3[i]-vn[2][i] < 0) {
    //                        qDebug() << v3[i]-vn[2][i];
    //                    }
    //                }


    //                m1 = green_diff(normi, 1);
    //                m1 = aubert_diff(norm_i_square, 1);
    //                m1 = li1(normi, 0.2);
    //                m1 = pm_diff(norm_i_square, 0.2);
    //                m1 = gr_diff(norm_i_square, 1);
    //                m1 = pm_diff2(norm_i_square, 0.2);

                    m1 = charbonnier_diff(norm_i_square, 0.65364);
                    m2 = 1;
                    m3 = 1;
                    m4 = sqrt(m1*m2);
                    //m4 = (m1+m2)/2;
                    //m4 = m1;
                    //m4 = 1;

                    m5 = sqrt(m1*m3);
    //                m5 = 0;

                    m6 = 1;
    //                m6 = 0;

//                   qDebug() << "Zeros Included - FOEED!";

    //                for(int k=0; k<3; k++) {
    //                    for(int i=0; i<3; i++) {
    //                        for(int j=0; j<3; j++) {
    //                            en[k][i][j] = v;
    //                        }
    //                    }
    //                }



    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e1[i][j] = v1[i]*v1[j];
    //                    }
    //                }
                    e1[0][0] = v1[0]*v1[0];
                    e1[0][1] = v1[0]*v1[1];
                    e1[0][2] = v1[0]*v1[2];
                    e1[1][0] = v1[1]*v1[0];
                    e1[1][1] = v1[1]*v1[1];
                    e1[1][2] = v1[1]*v1[2];
                    e1[2][0] = v1[2]*v1[0];
                    e1[2][1] = v1[2]*v1[1];
                    e1[2][2] = v1[2]*v1[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e2[i][j] = v2[i]*v2[j];
    //                    }
    //                }
                    e2[0][0] = v2[0]*v2[0];
                    e2[0][1] = v2[0]*v2[1];
                    e2[0][2] = v2[0]*v2[2];
                    e2[1][0] = v2[1]*v2[0];
                    e2[1][1] = v2[1]*v2[1];
                    e2[1][2] = v2[1]*v2[2];
                    e2[2][0] = v2[2]*v2[0];
                    e2[2][1] = v2[2]*v2[1];
                    e2[2][2] = v2[2]*v2[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e3[i][j] = v3[i]*v3[j];
    //                    }
    //                }
                    e3[0][0] = v3[0]*v3[0];
                    e3[0][1] = v3[0]*v3[1];
                    e3[0][2] = v3[0]*v3[2];
                    e3[1][0] = v3[1]*v3[0];
                    e3[1][1] = v3[1]*v3[1];
                    e3[1][2] = v3[1]*v3[2];
                    e3[2][0] = v3[2]*v3[0];
                    e3[2][1] = v3[2]*v3[1];
                    e3[2][2] = v3[2]*v3[2];

                    e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                    e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                    e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                    e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                    e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                    e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                    e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                    e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                    e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                    e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                    e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                    e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                    e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                    e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                    e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                    e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                    e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                    e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                    e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                    e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                    e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                    e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                    e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                    e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                    e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                    e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                    e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                    d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                    d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                    d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                    d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                    d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                    d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                    d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                    d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                    d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                    d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                    d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                    d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                    d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                    d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                    d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                    d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                    d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                    d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                    d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                    d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                    d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                    d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                    d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                    d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                    d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                    d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                    d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                    d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                    d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                    d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                    d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                    d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                    d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                    d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                    d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                    d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                    d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                    d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                    d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                    d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                    d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                    d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                    d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                    d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                    d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                    d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                    d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                    d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                    d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                    d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                    d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                    d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                    d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                    d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                    d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                    d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                    d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                    d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                    d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                    d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                    d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                    d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                    d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                    d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                    d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                    d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                    d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                    d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                    d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                    d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                    d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                    d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                    d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                    d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                    d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                    d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                    d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                    d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                    d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                    d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                    d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                    t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                    t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                    t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                    t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                    t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                    t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                    t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                    t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                    t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
                }

                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD21 = derivative3DX(t21,imgWidth,imgHeight,imgDepth,stepSize);
                dervXYD21 = derivative3DY(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD31 = derivative3DX(t31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXZD31 = derivative3DZ(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD12 = derivative3DY(t12,imgWidth,imgHeight,imgDepth,stepSize);
                dervYXD12 = derivative3DX(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD32 = derivative3DY(t32,imgWidth,imgHeight,imgDepth,stepSize);
                dervYZD32 = derivative3DZ(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD13 = derivative3DZ(t13,imgWidth,imgHeight,imgDepth,stepSize);
                dervZXD13 = derivative3DX(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD23 = derivative3DZ(t23,imgWidth,imgHeight,imgDepth,stepSize);
                dervZYD23 = derivative3DY(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == updatedPxlLocs[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStep*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
                }
                delete [] dervX;
                dervX = NULL;
                delete [] dervY;
                dervY = NULL;
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
                delete [] dervXY;
                dervXY = NULL;
                delete [] dervXZ;
                dervXZ = NULL;
                delete [] dervYZ;
                dervYZ = NULL;
                delete [] dervZZD33;
                dervZZD33 = NULL;
                delete [] dervXXD11;
                dervXXD11 = NULL;
                delete [] dervXYD21;
                dervXYD21 = NULL;
                delete [] dervYXD12;
                dervYXD12 = NULL;
                delete [] dervYD21;
                dervYD21 = NULL;
                delete [] dervZD31;
                dervZD31 = NULL;
                delete [] dervYYD22;
                dervYYD22 = NULL;
                delete [] dervXZD31;
                dervXZD31 = NULL;
                delete [] dervYZD32;
                dervYZD32 = NULL;
                delete [] dervZXD13;
                dervZXD13 = NULL;
                delete [] dervZYD23;
                dervZYD23 = NULL;
                delete [] dervXD12;
                dervXD12 = NULL;
                delete [] dervZD32;
                dervZD32 = NULL;
                delete [] dervXD13;
                dervXD13 = NULL;
                delete [] dervYD23;
                dervYD23 = NULL;
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
            }
    //        double l2normErrorOri = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
    //        printf("%lf\n", l2normErrorOri);
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    } else {
        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);    // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on z axis after x and y axis.

                //First derivatives of convolved mask
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.

                //Second derivatives
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
                dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
                dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xy.
                dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xz.
                dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
                dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. yz.


                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = scatImageArr[i];
                    }
                    tempImgArrayCurr[i] = scatImageArr[i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                    double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                    double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];
    //                double en[6][3][3], vn[3][3];

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

    //                    vn[0][0] = 1;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = 0;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 1;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 1;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = 0;
    //                    vn[0][1] = 0;
    //                    vn[0][2] = dervZConv[i]/normi;;

    //                    vn[1][0] = 0;
    //                    vn[1][1] = 0;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = 0;
    //                    vn[2][1] = 0;
    //                    vn[2][2] = 0;

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
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

    //                    vn[0][0] = dervXConv[i]/normi;
    //                    vn[0][1] = dervYConv[i]/normi;
    //                    vn[0][2] = dervZConv[i]/normi;

    //                    vn[1][0] = dervYConv[i]/xy_norm_i;
    //                    vn[1][1] = -dervXConv[i]/xy_norm_i;
    //                    vn[1][2] = 0;

    //                    vn[2][0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
    //                    vn[2][2] = -xy_norm_i_square/(xy_norm_i*normi);

    //                    if(v1[2]-vn[0][2] < 0) {
    //                        qDebug() << v1[2]-vn[0][2];
    //                    }

    //                    for(int i=0; i<3; i++) {
    //                        if(v1[i]-vn[0][i] < 0) {
    //                            qDebug() << v1[i]-vn[0][i];
    //                        }
    //                    }
                    }

    //                for(int i=0; i<3; i++) {
    //                    if(v1[i]-vn[0][i] < 0) {
    //                        qDebug() << v1[i]-vn[0][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v2[i]-vn[1][i] < 0) {
    //                        qDebug() << v2[i]-vn[1][i];
    //                    }
    //                }
    //                for(int i=0; i<3; i++) {
    //                    if(v3[i]-vn[2][i] < 0) {
    //                        qDebug() << v3[i]-vn[2][i];
    //                    }
    //                }


    //                m1 = green_diff(normi, 1);
    //                m1 = aubert_diff(norm_i_square, 1);
    //                m1 = li1(normi, 0.2);
    //                m1 = pm_diff(norm_i_square, 0.2);
    //                m1 = gr_diff(norm_i_square, 1);
    //                m1 = pm_diff2(norm_i_square, 0.2);

                    m1 = charbonnier_diff(norm_i_square, 0.561268);
                    m2 = 1;
                    m3 = 1;
                    m4 = sqrt(m1*m2);
                    //m4 = (m1+m2)/2;
                    //m4 = m1;
                    //m4 = 1;

                    m5 = sqrt(m1*m3);
    //                m5 = 0;

                    m6 = 1;
    //                m6 = 0;

//            qDebug() << "Ana";

    //                for(int k=0; k<3; k++) {
    //                    for(int i=0; i<3; i++) {
    //                        for(int j=0; j<3; j++) {
    //                            en[k][i][j] = v;
    //                        }
    //                    }
    //                }



    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e1[i][j] = v1[i]*v1[j];
    //                    }
    //                }
                    e1[0][0] = v1[0]*v1[0];
                    e1[0][1] = v1[0]*v1[1];
                    e1[0][2] = v1[0]*v1[2];
                    e1[1][0] = v1[1]*v1[0];
                    e1[1][1] = v1[1]*v1[1];
                    e1[1][2] = v1[1]*v1[2];
                    e1[2][0] = v1[2]*v1[0];
                    e1[2][1] = v1[2]*v1[1];
                    e1[2][2] = v1[2]*v1[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e2[i][j] = v2[i]*v2[j];
    //                    }
    //                }
                    e2[0][0] = v2[0]*v2[0];
                    e2[0][1] = v2[0]*v2[1];
                    e2[0][2] = v2[0]*v2[2];
                    e2[1][0] = v2[1]*v2[0];
                    e2[1][1] = v2[1]*v2[1];
                    e2[1][2] = v2[1]*v2[2];
                    e2[2][0] = v2[2]*v2[0];
                    e2[2][1] = v2[2]*v2[1];
                    e2[2][2] = v2[2]*v2[2];

    //                for(int i=0; i<3; i++) {
    //                    for(int j=0; j<3; j++) {
    //                        e3[i][j] = v3[i]*v3[j];
    //                    }
    //                }
                    e3[0][0] = v3[0]*v3[0];
                    e3[0][1] = v3[0]*v3[1];
                    e3[0][2] = v3[0]*v3[2];
                    e3[1][0] = v3[1]*v3[0];
                    e3[1][1] = v3[1]*v3[1];
                    e3[1][2] = v3[1]*v3[2];
                    e3[2][0] = v3[2]*v3[0];
                    e3[2][1] = v3[2]*v3[1];
                    e3[2][2] = v3[2]*v3[2];

                    e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                    e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                    e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                    e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                    e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                    e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                    e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                    e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                    e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                    e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                    e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                    e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                    e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                    e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                    e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                    e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                    e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                    e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                    e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                    e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                    e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                    e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                    e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                    e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                    e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                    e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                    e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                    d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                    d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                    d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                    d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                    d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                    d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                    d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                    d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                    d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                    d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                    d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                    d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                    d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                    d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                    d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                    d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                    d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                    d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                    d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                    d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                    d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                    d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                    d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                    d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                    d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                    d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                    d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                    d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                    d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                    d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                    d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                    d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                    d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                    d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                    d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                    d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                    d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                    d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                    d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                    d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                    d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                    d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                    d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                    d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                    d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                    d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                    d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                    d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                    d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                    d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                    d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                    d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                    d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                    d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                    d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                    d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                    d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                    d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                    d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                    d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                    d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                    d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                    d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                    d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                    d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                    d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                    d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                    d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                    d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                    d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                    d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                    d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                    d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                    d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                    d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                    d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                    d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                    d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                    d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                    d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                    d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                    t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                    t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                    t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                    t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                    t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                    t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                    t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                    t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                    t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
                }

                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD21 = derivative3DX(t21,imgWidth,imgHeight,imgDepth,stepSize);
                dervXYD21 = derivative3DY(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD31 = derivative3DX(t31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXZD31 = derivative3DZ(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD12 = derivative3DY(t12,imgWidth,imgHeight,imgDepth,stepSize);
                dervYXD12 = derivative3DX(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD32 = derivative3DY(t32,imgWidth,imgHeight,imgDepth,stepSize);
                dervYZD32 = derivative3DZ(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD13 = derivative3DZ(t13,imgWidth,imgHeight,imgDepth,stepSize);
                dervZXD13 = derivative3DX(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD23 = derivative3DZ(t23,imgWidth,imgHeight,imgDepth,stepSize);
                dervZYD23 = derivative3DY(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);

                alpha = (double)(4*n+2)/(2*n+3);

                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(i == randPxls[randArrTraceIndex]) {
                        randArrTraceIndex++;
                        continue;
                    }
                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStep*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
                }
                delete [] dervX;
                dervX = NULL;
                delete [] dervY;
                dervY = NULL;
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
                delete [] dervXY;
                dervXY = NULL;
                delete [] dervXZ;
                dervXZ = NULL;
                delete [] dervYZ;
                dervYZ = NULL;
                delete [] dervZZD33;
                dervZZD33 = NULL;
                delete [] dervXXD11;
                dervXXD11 = NULL;
                delete [] dervXYD21;
                dervXYD21 = NULL;
                delete [] dervYXD12;
                dervYXD12 = NULL;
                delete [] dervYD21;
                dervYD21 = NULL;
                delete [] dervZD31;
                dervZD31 = NULL;
                delete [] dervYYD22;
                dervYYD22 = NULL;
                delete [] dervXZD31;
                dervXZD31 = NULL;
                delete [] dervYZD32;
                dervYZD32 = NULL;
                delete [] dervZXD13;
                dervZXD13 = NULL;
                delete [] dervZYD23;
                dervZYD23 = NULL;
                delete [] dervXD12;
                dervXD12 = NULL;
                delete [] dervZD32;
                dervZD32 = NULL;
                delete [] dervXD13;
                dervXD13 = NULL;
                delete [] dervYD23;
                dervYD23 = NULL;
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
            }
    //        double l2normErrorOri = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
    //        printf("%lf\n", l2normErrorOri);
            l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
    }

    //Make nonnegative voxel values
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        if(scatImageArr[i] < 0) {
            scatImageArr[i] = 0;
        }
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] d1111;
    d1111 = NULL;
    delete [] d1112;
    d1112 = NULL;
    delete [] d1113;
    d1113 = NULL;
    delete [] d1121;
    d1121 = NULL;
    delete [] d1122;
    d1122 = NULL;
    delete [] d1123;
    d1123 = NULL;
    delete [] d1131;
    d1131 = NULL;
    delete [] d1132;
    d1132 = NULL;
    delete [] d1133;
    d1133 = NULL;

    delete [] d1211;
    d1211 = NULL;
    delete [] d1212;
    d1212 = NULL;
    delete [] d1213;
    d1213 = NULL;
    delete [] d1221;
    d1221 = NULL;
    delete [] d1222;
    d1222 = NULL;
    delete [] d1223;
    d1223 = NULL;
    delete [] d1231;
    d1231 = NULL;
    delete [] d1232;
    d1232 = NULL;
    delete [] d1233;
    d1233 = NULL;

    delete [] d1311;
    d1311 = NULL;
    delete [] d1312;
    d1312 = NULL;
    delete [] d1313;
    d1313 = NULL;
    delete [] d1321;
    d1321 = NULL;
    delete [] d1322;
    d1322 = NULL;
    delete [] d1323;
    d1323 = NULL;
    delete [] d1331;
    d1331 = NULL;
    delete [] d1332;
    d1332 = NULL;
    delete [] d1333;
    d1333 = NULL;


    delete [] d2111;
    d2111 = NULL;
    delete [] d2112;
    d2112 = NULL;
    delete [] d2113;
    d2113 = NULL;
    delete [] d2121;
    d2121 = NULL;
    delete [] d2122;
    d2122 = NULL;
    delete [] d2123;
    d2123 = NULL;
    delete [] d2131;
    d2131 = NULL;
    delete [] d2132;
    d2132 = NULL;
    delete [] d2133;
    d2133 = NULL;

    delete [] d2211;
    d2211 = NULL;
    delete [] d2212;
    d2212 = NULL;
    delete [] d2213;
    d2213 = NULL;
    delete [] d2221;
    d2221 = NULL;
    delete [] d2222;
    d2222 = NULL;
    delete [] d2223;
    d2223 = NULL;
    delete [] d2231;
    d2231 = NULL;
    delete [] d2232;
    d2232 = NULL;
    delete [] d2233;
    d2233 = NULL;

    delete [] d2311;
    d2311 = NULL;
    delete [] d2312;
    d2312 = NULL;
    delete [] d2313;
    d2313 = NULL;
    delete [] d2321;
    d2321 = NULL;
    delete [] d2322;
    d2322 = NULL;
    delete [] d2323;
    d2323 = NULL;
    delete [] d2331;
    d2331 = NULL;
    delete [] d2332;
    d2332 = NULL;
    delete [] d2333;
    d2333 = NULL;


    delete [] d3111;
    d3111 = NULL;
    delete [] d3112;
    d3112 = NULL;
    delete [] d3113;
    d3113 = NULL;
    delete [] d3121;
    d3121 = NULL;
    delete [] d3122;
    d3122 = NULL;
    delete [] d3123;
    d3123 = NULL;
    delete [] d3131;
    d3131 = NULL;
    delete [] d3132;
    d3132 = NULL;
    delete [] d3133;
    d3133 = NULL;

    delete [] d3211;
    d3211 = NULL;
    delete [] d3212;
    d3212 = NULL;
    delete [] d3213;
    d3213 = NULL;
    delete [] d3221;
    d3221 = NULL;
    delete [] d3222;
    d3222 = NULL;
    delete [] d3223;
    d3223 = NULL;
    delete [] d3231;
    d3231 = NULL;
    delete [] d3232;
    d3232 = NULL;
    delete [] d3233;
    d3233 = NULL;

    delete [] d3311;
    d3311 = NULL;
    delete [] d3312;
    d3312 = NULL;
    delete [] d3313;
    d3313 = NULL;
    delete [] d3321;
    d3321 = NULL;
    delete [] d3322;
    d3322 = NULL;
    delete [] d3323;
    d3323 = NULL;
    delete [] d3331;
    d3331 = NULL;
    delete [] d3332;
    d3332 = NULL;
    delete [] d3333;
    d3333 = NULL;

    delete [] t11;
    t11 = NULL;
    delete [] t12;
    t12 = NULL;
    delete [] t13;
    t13 = NULL;
    delete [] t21;
    t21 = NULL;
    delete [] t22;
    t22 = NULL;
    delete [] t23;
    t23 = NULL;
    delete [] t31;
    t31 = NULL;
    delete [] t32;
    t32 = NULL;
    delete [] t33;
    t33 = NULL;

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}



//***************************************************************************************************************************************************************************************************************************//
//4D image inpainting with EED
void eed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ) {
    clock_t begin = clock();

    int kernelSize = 5, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror=0, aaerror=0, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;
    qDebug() << gridSpcX << gridSpcY << gridSpcZ;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

//    imgTimeLen = 1;
    for(int t=0; t<imgTimeLen; t++) {
        qDebug() << t << ". volume";
        //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
        //    qDebug() << "MSE error: " << mserror << "\n";
        //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
        //    printf("Error: %lf\n", aaerror);
        l2normError = l2Norm(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "l2 norm error: " << l2normError << "\n";

        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                //First convolved derivatives
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

                //First derivatives
                dervX = derivative3DX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                dervY = derivative3DY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                dervZ = derivative3DZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                dervForX = dervForw3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);

    // //             Adaptive contrast parameter selection******************************
    //            double tv_norm = 0, contPar;
    //            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
    //                double norm_i_square, normi;

    //                norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i];
    //                normi = sqrt(norm_i_square);
    //                tv_norm = tv_norm + normi;
    //            }
    //            tv_norm = tv_norm/(imgWidth*imgHeight*imgDepth);
    //            contPar = 0.002*tv_norm;
    // //            qDebug() << tv_norm << contPar;
    // //            ********************************************************************

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = inpaintedImageArr[t][i];
                    }
                    tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
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

                    double contPar = 0.01;
                    d1 = charbonnier_diff(norm_i_square, contPar);
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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);
//                randArrTraceIndex = 0;
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(scatImageArr[t][i] == 0.0) {
//                        randArrTraceIndex++;
//                        continue;
//                        qDebug() << scatImageArr[t][i] << imageArr[t][i] << inpaintedImageArr[t][i];
                    } else {
                        inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                        tempImgArrayPrev[i] = tempImgArrayCurr[i];

//                        inpaintedImageArr[t][i] = 0;
                    }
//                    if(i == randPxls[randArrTraceIndex]) {
//                        randArrTraceIndex++;
//                        continue;
//                    }
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
            }
            l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
        mserror = mse(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "AAE error: " << aaerror << "\n";

        if(t+1<imgTimeLen) {
            l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
            qDebug() << "l2 norm error: " << l2normError << "\n";
        }
    }


    qDebug() << "Test (After): " << inpaintedImageArr[0][200772] << scatImageArr[0][200772] << imageArr[0][200772];



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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}



//4D image inpainting with FOEED
void foeed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int kernelSize = 5, N, randArrTraceIndex;
    float stepSize = 1;
    double sigma = 1.0, alpha, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervXX, *dervYY, *dervZZ, *dervXY, *dervYZ, *dervX, *dervY, *dervXZ, *dervXConv, *dervYConv, *dervZConv, *t11, *t12, *t13, *t21, *t22, *t23, *t31, *t32, *t33;
    double *d1111, *d1112, *d1113, *d1121, *d1122, *d1123, *d1131, *d1132, *d1133, *d1211, *d1212, *d1213, *d1221, *d1222, *d1223, *d1231, *d1232, *d1233, *d1311, *d1312, *d1313, *d1321, *d1322, *d1323, *d1331, *d1332, *d1333;
    double *d2111, *d2112, *d2113, *d2121, *d2122, *d2123, *d2131, *d2132, *d2133, *d2211, *d2212, *d2213, *d2221, *d2222, *d2223, *d2231, *d2232, *d2233, *d2311, *d2312, *d2313, *d2321, *d2322, *d2323, *d2331, *d2332, *d2333;
    double *d3111, *d3112, *d3113, *d3121, *d3122, *d3123, *d3131, *d3132, *d3133, *d3211, *d3212, *d3213, *d3221, *d3222, *d3223, *d3231, *d3232, *d3233, *d3311, *d3312, *d3313, *d3321, *d3322, *d3323, *d3331, *d3332, *d3333;
    double *dervXXD11, *dervYYD22, *dervZZD33, *dervYD21, *dervXYD21, *dervXD12, *dervYXD12, *dervZD31, *dervXZD31, *dervXD13, *dervZXD13, *dervYD23, *dervZYD23, *dervZD32, *dervYZD32;
    double *tempImgArrayCurr, *tempImgArrayPrev;

    qDebug() << "4D FOEED: " << "Kernel Size: " << kernelSize;

    //************************************************************
    //Memory allocation for 4th order diffusion tensor entries.
    d1111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    //************************************************************
//    double* te[3][3];
//    for(int i=0; i<3; i++) {
//        for(int j=0; j<3; j++) {
//            te[i][j] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
//        }
//    }

    t11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    //Memory allocation for temporary image array.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

//    imgTimeLen = 1;
    for(int t=0; t<imgTimeLen; t++) {
        qDebug() << t << ". volume";
        //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
        //    qDebug() << "MSE error: " << mserror << "\n";
        //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
        //    printf("Error: %lf\n", aaerror);
        l2normError = l2Norm(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "l2 norm error: " << l2normError << "\n";

        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.

                //First derivatives of convolved mask
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.

                //Second derivatives
                dervXX = derivative3DXX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
                dervX = derivative3DX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
                dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xy.
                dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xz.
                dervY = derivative3DY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
                dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. yz.


                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = inpaintedImageArr[t][i];
                    }
                    tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                    double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                    double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];
    //                double en[6][3][3], vn[3][3];

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

    //                m1 = green_diff(normi, 1);
    //                m1 = aubert_diff(norm_i_square, 1);
    //                m1 = li1(normi, 0.2);
    //                m1 = pm_diff(norm_i_square, 0.2);
    //                m1 = gr_diff(norm_i_square, 1);
    //                m1 = pm_diff2(norm_i_square, 0.2);
                    m1 = charbonnier_diff(norm_i_square, 0.01);
                    m2 = 1;
                    m3 = 1;
                    m4 = sqrt(m1*m2);
                    //m4 = (m1+m2)/2;
                    //m4 = m1;
                    //m4 = 1;
                    m5 = sqrt(m1*m3);
    //                m5 = 0;
                    m6 = 1;
    //                m6 = 0;

                    e1[0][0] = v1[0]*v1[0];
                    e1[0][1] = v1[0]*v1[1];
                    e1[0][2] = v1[0]*v1[2];
                    e1[1][0] = v1[1]*v1[0];
                    e1[1][1] = v1[1]*v1[1];
                    e1[1][2] = v1[1]*v1[2];
                    e1[2][0] = v1[2]*v1[0];
                    e1[2][1] = v1[2]*v1[1];
                    e1[2][2] = v1[2]*v1[2];

                    e2[0][0] = v2[0]*v2[0];
                    e2[0][1] = v2[0]*v2[1];
                    e2[0][2] = v2[0]*v2[2];
                    e2[1][0] = v2[1]*v2[0];
                    e2[1][1] = v2[1]*v2[1];
                    e2[1][2] = v2[1]*v2[2];
                    e2[2][0] = v2[2]*v2[0];
                    e2[2][1] = v2[2]*v2[1];
                    e2[2][2] = v2[2]*v2[2];

                    e3[0][0] = v3[0]*v3[0];
                    e3[0][1] = v3[0]*v3[1];
                    e3[0][2] = v3[0]*v3[2];
                    e3[1][0] = v3[1]*v3[0];
                    e3[1][1] = v3[1]*v3[1];
                    e3[1][2] = v3[1]*v3[2];
                    e3[2][0] = v3[2]*v3[0];
                    e3[2][1] = v3[2]*v3[1];
                    e3[2][2] = v3[2]*v3[2];

                    e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                    e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                    e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                    e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                    e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                    e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                    e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                    e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                    e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                    e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                    e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                    e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                    e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                    e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                    e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                    e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                    e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                    e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                    e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                    e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                    e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                    e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                    e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                    e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                    e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                    e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                    e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                    d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                    d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                    d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                    d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                    d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                    d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                    d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                    d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                    d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                    d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                    d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                    d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                    d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                    d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                    d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                    d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                    d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                    d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                    d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                    d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                    d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                    d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                    d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                    d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                    d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                    d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                    d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                    d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                    d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                    d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                    d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                    d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                    d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                    d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                    d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                    d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                    d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                    d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                    d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                    d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                    d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                    d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                    d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                    d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                    d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                    d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                    d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                    d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                    d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                    d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                    d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                    d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                    d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                    d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                    d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                    d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                    d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                    d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                    d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                    d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                    d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                    d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                    d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                    d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                    d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                    d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                    d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                    d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                    d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                    d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                    d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                    d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                    d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                    d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                    d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                    d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                    d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                    d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                    d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                    d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                    d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                    t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                    t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                    t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                    t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                    t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                    t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                    t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                    t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                    t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
                }

                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD21 = derivative3DX(t21,imgWidth,imgHeight,imgDepth,stepSize);
                dervXYD21 = derivative3DY(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD31 = derivative3DX(t31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXZD31 = derivative3DZ(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD12 = derivative3DY(t12,imgWidth,imgHeight,imgDepth,stepSize);
                dervYXD12 = derivative3DX(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD32 = derivative3DY(t32,imgWidth,imgHeight,imgDepth,stepSize);
                dervYZD32 = derivative3DZ(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD13 = derivative3DZ(t13,imgWidth,imgHeight,imgDepth,stepSize);
                dervZXD13 = derivative3DX(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD23 = derivative3DZ(t23,imgWidth,imgHeight,imgDepth,stepSize);
                dervZYD23 = derivative3DY(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);


                alpha = (double)(4*n+2)/(2*n+3);
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(scatImageArr[t][i] == 0.0) {
//                        qDebug() << scatImageArr[t][i] << imageArr[t][i] << inpaintedImageArr[t][i];
                    } else {
                        inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] - timeStepSize*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                        tempImgArrayPrev[i] = tempImgArrayCurr[i];
//                        inpaintedImageArr[t][i] = 0;
                    }
                }
                delete [] dervX;
                dervX = NULL;
                delete [] dervY;
                dervY = NULL;
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
                delete [] dervXY;
                dervXY = NULL;
                delete [] dervXZ;
                dervXZ = NULL;
                delete [] dervYZ;
                dervYZ = NULL;
                delete [] dervZZD33;
                dervZZD33 = NULL;
                delete [] dervXXD11;
                dervXXD11 = NULL;
                delete [] dervXYD21;
                dervXYD21 = NULL;
                delete [] dervYXD12;
                dervYXD12 = NULL;
                delete [] dervYD21;
                dervYD21 = NULL;
                delete [] dervZD31;
                dervZD31 = NULL;
                delete [] dervYYD22;
                dervYYD22 = NULL;
                delete [] dervXZD31;
                dervXZD31 = NULL;
                delete [] dervYZD32;
                dervYZD32 = NULL;
                delete [] dervZXD13;
                dervZXD13 = NULL;
                delete [] dervZYD23;
                dervZYD23 = NULL;
                delete [] dervXD12;
                dervXD12 = NULL;
                delete [] dervZD32;
                dervZD32 = NULL;
                delete [] dervXD13;
                dervXD13 = NULL;
                delete [] dervYD23;
                dervYD23 = NULL;
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
            }
            l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
        mserror = mse(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "AAE error: " << aaerror << "\n";

        if(t+1<imgTimeLen) {
            l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
            qDebug() << "l2 norm error: " << l2normError << "\n";
        }
    }

    delete [] d1111;
    d1111 = NULL;
    delete [] d1112;
    d1112 = NULL;
    delete [] d1113;
    d1113 = NULL;
    delete [] d1121;
    d1121 = NULL;
    delete [] d1122;
    d1122 = NULL;
    delete [] d1123;
    d1123 = NULL;
    delete [] d1131;
    d1131 = NULL;
    delete [] d1132;
    d1132 = NULL;
    delete [] d1133;
    d1133 = NULL;

    delete [] d1211;
    d1211 = NULL;
    delete [] d1212;
    d1212 = NULL;
    delete [] d1213;
    d1213 = NULL;
    delete [] d1221;
    d1221 = NULL;
    delete [] d1222;
    d1222 = NULL;
    delete [] d1223;
    d1223 = NULL;
    delete [] d1231;
    d1231 = NULL;
    delete [] d1232;
    d1232 = NULL;
    delete [] d1233;
    d1233 = NULL;

    delete [] d1311;
    d1311 = NULL;
    delete [] d1312;
    d1312 = NULL;
    delete [] d1313;
    d1313 = NULL;
    delete [] d1321;
    d1321 = NULL;
    delete [] d1322;
    d1322 = NULL;
    delete [] d1323;
    d1323 = NULL;
    delete [] d1331;
    d1331 = NULL;
    delete [] d1332;
    d1332 = NULL;
    delete [] d1333;
    d1333 = NULL;


    delete [] d2111;
    d2111 = NULL;
    delete [] d2112;
    d2112 = NULL;
    delete [] d2113;
    d2113 = NULL;
    delete [] d2121;
    d2121 = NULL;
    delete [] d2122;
    d2122 = NULL;
    delete [] d2123;
    d2123 = NULL;
    delete [] d2131;
    d2131 = NULL;
    delete [] d2132;
    d2132 = NULL;
    delete [] d2133;
    d2133 = NULL;

    delete [] d2211;
    d2211 = NULL;
    delete [] d2212;
    d2212 = NULL;
    delete [] d2213;
    d2213 = NULL;
    delete [] d2221;
    d2221 = NULL;
    delete [] d2222;
    d2222 = NULL;
    delete [] d2223;
    d2223 = NULL;
    delete [] d2231;
    d2231 = NULL;
    delete [] d2232;
    d2232 = NULL;
    delete [] d2233;
    d2233 = NULL;

    delete [] d2311;
    d2311 = NULL;
    delete [] d2312;
    d2312 = NULL;
    delete [] d2313;
    d2313 = NULL;
    delete [] d2321;
    d2321 = NULL;
    delete [] d2322;
    d2322 = NULL;
    delete [] d2323;
    d2323 = NULL;
    delete [] d2331;
    d2331 = NULL;
    delete [] d2332;
    d2332 = NULL;
    delete [] d2333;
    d2333 = NULL;


    delete [] d3111;
    d3111 = NULL;
    delete [] d3112;
    d3112 = NULL;
    delete [] d3113;
    d3113 = NULL;
    delete [] d3121;
    d3121 = NULL;
    delete [] d3122;
    d3122 = NULL;
    delete [] d3123;
    d3123 = NULL;
    delete [] d3131;
    d3131 = NULL;
    delete [] d3132;
    d3132 = NULL;
    delete [] d3133;
    d3133 = NULL;

    delete [] d3211;
    d3211 = NULL;
    delete [] d3212;
    d3212 = NULL;
    delete [] d3213;
    d3213 = NULL;
    delete [] d3221;
    d3221 = NULL;
    delete [] d3222;
    d3222 = NULL;
    delete [] d3223;
    d3223 = NULL;
    delete [] d3231;
    d3231 = NULL;
    delete [] d3232;
    d3232 = NULL;
    delete [] d3233;
    d3233 = NULL;

    delete [] d3311;
    d3311 = NULL;
    delete [] d3312;
    d3312 = NULL;
    delete [] d3313;
    d3313 = NULL;
    delete [] d3321;
    d3321 = NULL;
    delete [] d3322;
    d3322 = NULL;
    delete [] d3323;
    d3323 = NULL;
    delete [] d3331;
    d3331 = NULL;
    delete [] d3332;
    d3332 = NULL;
    delete [] d3333;
    d3333 = NULL;

    delete [] t11;
    t11 = NULL;
    delete [] t12;
    t12 = NULL;
    delete [] t13;
    t13 = NULL;
    delete [] t21;
    t21 = NULL;
    delete [] t22;
    t22 = NULL;
    delete [] t23;
    t23 = NULL;
    delete [] t31;
    t31 = NULL;
    delete [] t32;
    t32 = NULL;
    delete [] t33;
    t33 = NULL;

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//4D image inpainting with Linear Homogenous Diffusion
void linear_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int stepSize = 1, N;
    double mserror, aaerror, l2normError, alpha;
    double *dervYY, *dervXX, *dervZZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gridSpcX = gridSpcY = gridSpcZ = 1;
    qDebug() << gridSpcX << gridSpcY << gridSpcZ;

    for(int t=0; t<imgTimeLen; t++) {
        qDebug() << "Volume No: " << t;
        //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
        //    qDebug() << "MSE error: " << mserror << "\n";
        //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
        //    printf("Error: %lf\n", aaerror);
        l2normError = l2Norm(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "l2 norm error: " << l2normError << "\n";

        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                //First derivatives
                dervXX = derivative3DXX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.

                alpha = (double)(4*n+2)/(2*n+3);
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = inpaintedImageArr[t][i];
                    }
                    tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                    if(scatImageArr[t][i] == 0.0) {
//                        randArrTraceIndex++;
//                        continue;
//                        qDebug() << scatImageArr[t][i] << imageArr[t][i] << inpaintedImageArr[t][i];
                    } else {
                        inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] + timeStepSize*(dervXX[i]+dervYY[i]+dervZZ[i])) + (1-alpha)*tempImgArrayPrev[i];
                        tempImgArrayPrev[i] = tempImgArrayCurr[i];
                    }
                }
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
            }
            l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
        mserror = mse(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = aae(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "AAE error: " << aaerror << "\n";

        if(t+1<imgTimeLen) {
            l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
            qDebug() << "l2 norm error: " << l2normError << "\n";
        }
    }


    qDebug() << "Test (After): " << inpaintedImageArr[0][200772] << scatImageArr[0][200772] << imageArr[0][200772];

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//***************************************************************************************************************************************************************************************************************************//
//4D image signal dropout imputation with EED
void eed_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int kernelSize = 3, N;
//    float stepSize = 1;
    double sigma = 0.5, gausKernel[kernelSize], mserror=0, aaerror=0, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;
    qDebug() << gridSpcX << gridSpcY << gridSpcZ;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    vector<int> imgVolLine = corruptedVolumes(scatImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen);
//    vector<int> imgVolLine;
//    imgVolLine.push_back(2);

    for(int i=0; i<imgVolLine.size(); i++) {
        qDebug() << imgVolLine[i] << " ";
    }    

// //#######################################################################################################################
//    // Make regular grid mask
//    for(int v=0; v<imgVolLine.size(); v++) {
//        int t = imgVolLine[v];
//        qDebug() << "Volume No: " << t;

//        vector<int> imgSlcLine = corruptedSlices(scatImageArr[t], imgWidth, imgHeight, imgDepth);
// //        for(int i=0; i<imgSlcLine.size(); i++) {
// //            qDebug() << imgSlcLine[i] << " ";
// //        }

//        for(int s=0; s<imgSlcLine.size(); s++) {
//            makeRegGrid_on_slice(scatImageArr[t], imgWidth, imgHeight, imgSlcLine[s], 2, 0);
//        }
//    }
//    return;
// //#######################################################################################################################


    double tol_mserror = 0, tol_aaerror = 0;
    for(int v=0; v<imgVolLine.size(); v++) {
        int t = imgVolLine[v];
 //        for(int t=22; t<imgTimeLen; t++) {
            qDebug() << "Volume No: " << t;
            //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
            //    qDebug() << "MSE error: " << mserror << "\n";
            //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
            //    printf("Error: %lf\n", aaerror);
 //            l2normError = l2Norm(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
            l2normError = 1000;
            qDebug() << "l2 norm error: " << l2normError << "\n";

            double contPar = 0.75;
    //        if(t == 0) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 21) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 42) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 63) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 84) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else {
    //            contPar = 0.0005;            // B=700 or B=2000 optimal contrast parameter
    //        }

 //            //Contrast Paremeter estimation
 //            double quantile = 90.0;
 //            double estContPar;
 //            estContPar = quantlCriterPM4Inpainting4D(quantile, imageArr[t], scatImageArr[t], imgWidth, imgHeight, imgDepth, sigma, kernelSize);
 //            estContPar = estContPar/25;
 //            qDebug() << "Estimated Contrast Parameter: " << estContPar;

            int stopTime = 0;
            N = numSteps;
            while(l2normError > tol) {
    //            if(stopTime > 3) {
    //                break;
    //            }
                for(int n=0; n<N; n++) {
                    outConvX = convolution3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                    outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                    outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                    //First convolved derivatives
                    dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                    dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                    dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

                    //First derivatives
                    dervX = derivative3DX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                    dervY = derivative3DY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                    dervZ = derivative3DZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                    dervForX = dervForw3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                    dervForY = dervForw3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                    dervForZ = dervForw3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);
                    dervBacX = dervBack3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                    dervBacY = dervBack3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                    dervBacZ = dervBack3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);

                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(n == 0) {
                            tempImgArrayPrev[i] = inpaintedImageArr[t][i];
                        }
                        tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                        //Define eigenvectors and entries for the diffusion tensor.
                        double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
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

 //                        d1 = pm_diff2(norm_i_square, contPar);
                        d1 = charbonnier_diff(norm_i_square, contPar);
 //                        d1 = charbonnier_diff(norm_i_square, estContPar);
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

                    dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                    dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                    dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                    dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                    dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                    dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                    alpha = (double)(4*n+2)/(2*n+3);
                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(scatImageArr[t][i] == 0.0) {
    //                        randArrTraceIndex++;
    //                        continue;
    //                        qDebug() << scatImageArr[t][i] << imageArr[t][i] << inpaintedImageArr[t][i];
                        } else {
                            inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                            tempImgArrayPrev[i] = tempImgArrayCurr[i];
                        }
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
                }

                stopTime++;

                l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
                qDebug() << l2normError << "\n";
            }
    //        mserror = mse(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
    //        qDebug() << "MSE error: " << mserror << "\n";
    //        aaerror = aae(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
    //        qDebug() << "AAE error: " << aaerror << "\n";

            mserror = sub_mse(refImageArr[t],inpaintedImageArr[t], scatImageArr[t],imgWidth*imgHeight*imgDepth);
            qDebug() << "MSE error: " << mserror << "\n";
            tol_mserror = tol_mserror + mserror;
            aaerror = sub_aae(refImageArr[t],inpaintedImageArr[t], scatImageArr[t], imgWidth*imgHeight*imgDepth);
            qDebug() << "AAE error: " << aaerror << "\n";
            tol_aaerror = tol_aaerror + aaerror;

            if(t+1<imgTimeLen) {
                l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
                qDebug() << "l2 norm error: " << l2normError << "\n";
            }
 //        }
    }
    tol_aaerror = sub_aae_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average AAE error: " << tol_aaerror << "\n";
    tol_mserror = sub_mse_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average MSE error: " << tol_mserror << "\n";


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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//4D image signal dropout imputation with FOEED
void foeed_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int kernelSize = 3, N;
    float stepSize = 1;
    double sigma = 1.0, alpha, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervXX, *dervYY, *dervZZ, *dervXY, *dervYZ, *dervX, *dervY, *dervXZ, *dervXConv, *dervYConv, *dervZConv, *t11, *t12, *t13, *t21, *t22, *t23, *t31, *t32, *t33;
    double *d1111, *d1112, *d1113, *d1121, *d1122, *d1123, *d1131, *d1132, *d1133, *d1211, *d1212, *d1213, *d1221, *d1222, *d1223, *d1231, *d1232, *d1233, *d1311, *d1312, *d1313, *d1321, *d1322, *d1323, *d1331, *d1332, *d1333;
    double *d2111, *d2112, *d2113, *d2121, *d2122, *d2123, *d2131, *d2132, *d2133, *d2211, *d2212, *d2213, *d2221, *d2222, *d2223, *d2231, *d2232, *d2233, *d2311, *d2312, *d2313, *d2321, *d2322, *d2323, *d2331, *d2332, *d2333;
    double *d3111, *d3112, *d3113, *d3121, *d3122, *d3123, *d3131, *d3132, *d3133, *d3211, *d3212, *d3213, *d3221, *d3222, *d3223, *d3231, *d3232, *d3233, *d3311, *d3312, *d3313, *d3321, *d3322, *d3323, *d3331, *d3332, *d3333;
    double *dervXXD11, *dervYYD22, *dervZZD33, *dervYD21, *dervXYD21, *dervXD12, *dervYXD12, *dervZD31, *dervXZD31, *dervXD13, *dervZXD13, *dervYD23, *dervZYD23, *dervZD32, *dervYZD32;
    double *tempImgArrayCurr, *tempImgArrayPrev;

    qDebug() << "4D FOEED: " << "Kernel Size: " << kernelSize;

    //************************************************************
    //Memory allocation for 4th order diffusion tensor entries.
    d1111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    //************************************************************
//    double* te[3][3];
//    for(int i=0; i<3; i++) {
//        for(int j=0; j<3; j++) {
//            te[i][j] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
//        }
//    }

    t11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    //Memory allocation for temporary image array.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.


    vector<int> imgVolLine = corruptedVolumes(scatImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen);
    for(int i=0; i<imgVolLine.size(); i++) {
        qDebug() << imgVolLine[i] << " ";
    }

    double tol_mserror = 0, tol_aaerror = 0;

    for(int v=0; v<imgVolLine.size(); v++) {
        int t = imgVolLine[v];
//    imgTimeLen = 23;
//    for(int t=22; t<imgTimeLen; t++) {
        qDebug() << "Volume No: " << t;
        //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
        //    qDebug() << "MSE error: " << mserror << "\n";
        //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
        //    printf("Error: %lf\n", aaerror);
        l2normError = l2Norm(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "l2 norm error: " << l2normError << "\n";

        double contPar = 1.0;
//        if(t == 0) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 21) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 42) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 63) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 84) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else {
//            contPar = 0.0005;            // B=700 or B=2000 optimal contrast parameter
//        }

        N = numSteps;
        while(l2normError > tol) {
            for(int n=0; n<N; n++) {
                outConvX = convolution3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);    // Calculate the convolution on x axis.
                outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);               // Calculate the convolution on y axis after x axis.
                outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);             // Calculate the convolution on z axis after x and y axis.

                //First derivatives of convolved mask
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                                // Derivative of convolved image w.r.t. x.
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                                // Derivative of convolved image w.r.t. y.
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                                // Derivative of convolved image w.r.t. z.

                //Second derivatives
                dervXX = derivative3DXX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
                dervYY = derivative3DYY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
                dervZZ = derivative3DZZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
                dervX = derivative3DX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
                dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                         // Derivative w.r.t. xy.
                dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                         // Derivative w.r.t. xz.
                dervY = derivative3DY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
                dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                         // Derivative w.r.t. yz.

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[i] = inpaintedImageArr[t][i];
                    }
                    tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                    double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                    double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                    double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];

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

    //                m1 = green_diff(normi, 1);
    //                m1 = aubert_diff(norm_i_square, 1);
    //                m1 = li1(normi, 0.2);
    //                m1 = pm_diff(norm_i_square, 0.2);
    //                m1 = gr_diff(norm_i_square, 1);
    //                m1 = pm_diff2(norm_i_square, 0.2);
                    m1 = charbonnier_diff(norm_i_square, contPar);
                    m2 = 1;
                    m3 = 1;
                    m4 = sqrt(m1*m2);
                    //m4 = (m1+m2)/2;
                    //m4 = m1;
                    //m4 = 1;
                    m5 = sqrt(m1*m3);
    //                m5 = 0;
                    m6 = 1;
    //                m6 = 0;

                    e1[0][0] = v1[0]*v1[0];
                    e1[0][1] = v1[0]*v1[1];
                    e1[0][2] = v1[0]*v1[2];
                    e1[1][0] = v1[1]*v1[0];
                    e1[1][1] = v1[1]*v1[1];
                    e1[1][2] = v1[1]*v1[2];
                    e1[2][0] = v1[2]*v1[0];
                    e1[2][1] = v1[2]*v1[1];
                    e1[2][2] = v1[2]*v1[2];

                    e2[0][0] = v2[0]*v2[0];
                    e2[0][1] = v2[0]*v2[1];
                    e2[0][2] = v2[0]*v2[2];
                    e2[1][0] = v2[1]*v2[0];
                    e2[1][1] = v2[1]*v2[1];
                    e2[1][2] = v2[1]*v2[2];
                    e2[2][0] = v2[2]*v2[0];
                    e2[2][1] = v2[2]*v2[1];
                    e2[2][2] = v2[2]*v2[2];

                    e3[0][0] = v3[0]*v3[0];
                    e3[0][1] = v3[0]*v3[1];
                    e3[0][2] = v3[0]*v3[2];
                    e3[1][0] = v3[1]*v3[0];
                    e3[1][1] = v3[1]*v3[1];
                    e3[1][2] = v3[1]*v3[2];
                    e3[2][0] = v3[2]*v3[0];
                    e3[2][1] = v3[2]*v3[1];
                    e3[2][2] = v3[2]*v3[2];

                    e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                    e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                    e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                    e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                    e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                    e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                    e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                    e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                    e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                    e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                    e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                    e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                    e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                    e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                    e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                    e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                    e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                    e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                    e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                    e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                    e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                    e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                    e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                    e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                    e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                    e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                    e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                    d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                    d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                    d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                    d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                    d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                    d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                    d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                    d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                    d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                    d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                    d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                    d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                    d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                    d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                    d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                    d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                    d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                    d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                    d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                    d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                    d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                    d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                    d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                    d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                    d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                    d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                    d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                    d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                    d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                    d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                    d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                    d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                    d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                    d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                    d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                    d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                    d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                    d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                    d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                    d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                    d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                    d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                    d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                    d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                    d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                    d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                    d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                    d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                    d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                    d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                    d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                    d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                    d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                    d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                    d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                    d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                    d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                    d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                    d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                    d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                    d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                    d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                    d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                    d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                    d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                    d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                    d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                    d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                    d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                    d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                    d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                    d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                    d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                    d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                    d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                    d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                    d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                    d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                    d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                    d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                    d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                    t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                    t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                    t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                    t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                    t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                    t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                    t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                    t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                    t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
                }

                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD21 = derivative3DX(t21,imgWidth,imgHeight,imgDepth,stepSize);
                dervXYD21 = derivative3DY(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD31 = derivative3DX(t31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXZD31 = derivative3DZ(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD12 = derivative3DY(t12,imgWidth,imgHeight,imgDepth,stepSize);
                dervYXD12 = derivative3DX(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
                dervZD32 = derivative3DY(t32,imgWidth,imgHeight,imgDepth,stepSize);
                dervYZD32 = derivative3DZ(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
                dervXD13 = derivative3DZ(t13,imgWidth,imgHeight,imgDepth,stepSize);
                dervZXD13 = derivative3DX(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
                dervYD23 = derivative3DZ(t23,imgWidth,imgHeight,imgDepth,stepSize);
                dervZYD23 = derivative3DY(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);


                alpha = (double)(4*n+2)/(2*n+3);
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(scatImageArr[t][i] == 0.0) {
//                        qDebug() << scatImageArr[t][i] << imageArr[t][i] << inpaintedImageArr[t][i];
                    } else {
                        inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] - timeStepSize*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                        tempImgArrayPrev[i] = tempImgArrayCurr[i];
                    }
                }
                delete [] dervX;
                dervX = NULL;
                delete [] dervY;
                dervY = NULL;
                delete [] dervXX;
                dervXX = NULL;
                delete [] dervYY;
                dervYY = NULL;
                delete [] dervZZ;
                dervZZ = NULL;
                delete [] dervXY;
                dervXY = NULL;
                delete [] dervXZ;
                dervXZ = NULL;
                delete [] dervYZ;
                dervYZ = NULL;
                delete [] dervZZD33;
                dervZZD33 = NULL;
                delete [] dervXXD11;
                dervXXD11 = NULL;
                delete [] dervXYD21;
                dervXYD21 = NULL;
                delete [] dervYXD12;
                dervYXD12 = NULL;
                delete [] dervYD21;
                dervYD21 = NULL;
                delete [] dervZD31;
                dervZD31 = NULL;
                delete [] dervYYD22;
                dervYYD22 = NULL;
                delete [] dervXZD31;
                dervXZD31 = NULL;
                delete [] dervYZD32;
                dervYZD32 = NULL;
                delete [] dervZXD13;
                dervZXD13 = NULL;
                delete [] dervZYD23;
                dervZYD23 = NULL;
                delete [] dervXD12;
                dervXD12 = NULL;
                delete [] dervZD32;
                dervZD32 = NULL;
                delete [] dervXD13;
                dervXD13 = NULL;
                delete [] dervYD23;
                dervYD23 = NULL;
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
            }
            l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
            qDebug() << l2normError << "\n";
        }
//        mserror = mse(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
//        qDebug() << "MSE error: " << mserror << "\n";
//        aaerror = aae(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
//        qDebug() << "AAE error: " << aaerror << "\n";
        mserror = sub_mse(refImageArr[t],inpaintedImageArr[t], scatImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = sub_aae(refImageArr[t],inpaintedImageArr[t], scatImageArr[t], imgWidth*imgHeight*imgDepth);
        qDebug() << "AAE error: " << aaerror << "\n";


        if(t+1<imgTimeLen) {
            l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
            qDebug() << "l2 norm error: " << l2normError << "\n";
        }
//    }
    }
    tol_aaerror = sub_aae_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average AAE error: " << tol_aaerror << "\n";
    tol_mserror = sub_mse_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average MSE error: " << tol_mserror << "\n";

    delete [] d1111;
    d1111 = NULL;
    delete [] d1112;
    d1112 = NULL;
    delete [] d1113;
    d1113 = NULL;
    delete [] d1121;
    d1121 = NULL;
    delete [] d1122;
    d1122 = NULL;
    delete [] d1123;
    d1123 = NULL;
    delete [] d1131;
    d1131 = NULL;
    delete [] d1132;
    d1132 = NULL;
    delete [] d1133;
    d1133 = NULL;

    delete [] d1211;
    d1211 = NULL;
    delete [] d1212;
    d1212 = NULL;
    delete [] d1213;
    d1213 = NULL;
    delete [] d1221;
    d1221 = NULL;
    delete [] d1222;
    d1222 = NULL;
    delete [] d1223;
    d1223 = NULL;
    delete [] d1231;
    d1231 = NULL;
    delete [] d1232;
    d1232 = NULL;
    delete [] d1233;
    d1233 = NULL;

    delete [] d1311;
    d1311 = NULL;
    delete [] d1312;
    d1312 = NULL;
    delete [] d1313;
    d1313 = NULL;
    delete [] d1321;
    d1321 = NULL;
    delete [] d1322;
    d1322 = NULL;
    delete [] d1323;
    d1323 = NULL;
    delete [] d1331;
    d1331 = NULL;
    delete [] d1332;
    d1332 = NULL;
    delete [] d1333;
    d1333 = NULL;


    delete [] d2111;
    d2111 = NULL;
    delete [] d2112;
    d2112 = NULL;
    delete [] d2113;
    d2113 = NULL;
    delete [] d2121;
    d2121 = NULL;
    delete [] d2122;
    d2122 = NULL;
    delete [] d2123;
    d2123 = NULL;
    delete [] d2131;
    d2131 = NULL;
    delete [] d2132;
    d2132 = NULL;
    delete [] d2133;
    d2133 = NULL;

    delete [] d2211;
    d2211 = NULL;
    delete [] d2212;
    d2212 = NULL;
    delete [] d2213;
    d2213 = NULL;
    delete [] d2221;
    d2221 = NULL;
    delete [] d2222;
    d2222 = NULL;
    delete [] d2223;
    d2223 = NULL;
    delete [] d2231;
    d2231 = NULL;
    delete [] d2232;
    d2232 = NULL;
    delete [] d2233;
    d2233 = NULL;

    delete [] d2311;
    d2311 = NULL;
    delete [] d2312;
    d2312 = NULL;
    delete [] d2313;
    d2313 = NULL;
    delete [] d2321;
    d2321 = NULL;
    delete [] d2322;
    d2322 = NULL;
    delete [] d2323;
    d2323 = NULL;
    delete [] d2331;
    d2331 = NULL;
    delete [] d2332;
    d2332 = NULL;
    delete [] d2333;
    d2333 = NULL;


    delete [] d3111;
    d3111 = NULL;
    delete [] d3112;
    d3112 = NULL;
    delete [] d3113;
    d3113 = NULL;
    delete [] d3121;
    d3121 = NULL;
    delete [] d3122;
    d3122 = NULL;
    delete [] d3123;
    d3123 = NULL;
    delete [] d3131;
    d3131 = NULL;
    delete [] d3132;
    d3132 = NULL;
    delete [] d3133;
    d3133 = NULL;

    delete [] d3211;
    d3211 = NULL;
    delete [] d3212;
    d3212 = NULL;
    delete [] d3213;
    d3213 = NULL;
    delete [] d3221;
    d3221 = NULL;
    delete [] d3222;
    d3222 = NULL;
    delete [] d3223;
    d3223 = NULL;
    delete [] d3231;
    d3231 = NULL;
    delete [] d3232;
    d3232 = NULL;
    delete [] d3233;
    d3233 = NULL;

    delete [] d3311;
    d3311 = NULL;
    delete [] d3312;
    d3312 = NULL;
    delete [] d3313;
    d3313 = NULL;
    delete [] d3321;
    d3321 = NULL;
    delete [] d3322;
    d3322 = NULL;
    delete [] d3323;
    d3323 = NULL;
    delete [] d3331;
    d3331 = NULL;
    delete [] d3332;
    d3332 = NULL;
    delete [] d3333;
    d3333 = NULL;

    delete [] t11;
    t11 = NULL;
    delete [] t12;
    t12 = NULL;
    delete [] t13;
    t13 = NULL;
    delete [] t21;
    t21 = NULL;
    delete [] t22;
    t22 = NULL;
    delete [] t23;
    t23 = NULL;
    delete [] t31;
    t31 = NULL;
    delete [] t32;
    t32 = NULL;
    delete [] t33;
    t33 = NULL;

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//4D image signal dropout imputation with EED Explicit Scheme
void eed_4d_signalDropoutImputation_ExplicitScheme(float tol, float timeStepSize, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int kernelSize = 3;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror=0, aaerror=0, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;
    qDebug() << gridSpcX << gridSpcY << gridSpcZ;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.


    //Random Selection of 20 voxels from dropouts (valid only to the current simulated mask data!) ******************************************************
    int locLen = 20, volNum = 22, upperBound = imgWidth*imgHeight*imgDepth;
    double* randLocs = randVoxLocs(locLen, scatImageArr, volNum, upperBound);
    for(int i=0; i<locLen; i++) {
        int x,y,z;
        z = floor(randLocs[i]/(imgWidth*imgHeight));
        y = floor((randLocs[i] - z*imgWidth*imgHeight)/imgWidth);
        x = floor(randLocs[i] - y*imgWidth - z*imgWidth*imgHeight);
        qDebug() << randLocs[i] << x << y << z;
    }
    // **************************************************************************************************************************************************


    imgTimeLen = 23;
    for(int t=22; t<imgTimeLen; t++) {
        qDebug() << t << ". volume";
        //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
        //    qDebug() << "MSE error: " << mserror << "\n";
        //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
        //    printf("Error: %lf\n", aaerror);
//        l2normError = l2Norm(imageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
        l2normError = 1000;
        qDebug() << "l2 norm error: " << l2normError << "\n";

        double contPar = 0.25;
//        if(t == 0) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 21) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 42) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 63) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 84) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else {
//            contPar = 0.0005;            // B=700 or B=2000 optimal contrast parameter
//        }


        int l = 19;
        double minErr=abs(refImageArr[t][int(randLocs[l])]-inpaintedImageArr[t][int(randLocs[l])]), minErrIter=0;
        double initErr = minErr;

        int stopTime = 0;
        while(l2normError > tol) {
//            if(stopTime > 14) {
//                break;
//            }
            outConvX = convolution3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
            outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
            outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
            //First convolved derivatives
            dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
            dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
            dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

            //First derivatives
            dervX = derivative3DX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
            dervY = derivative3DY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
            dervZ = derivative3DZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

            dervForX = dervForw3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
            dervForY = dervForw3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
            dervForZ = dervForw3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);
            dervBacX = dervBack3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
            dervBacY = dervBack3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
            dervBacZ = dervBack3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);

            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                //Define eigenvectors and entries for the diffusion tensor.
                double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
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
                d1 = charbonnier_diff(norm_i_square, contPar);
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

            dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
            dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
            dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
            dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
            dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
            dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                if(scatImageArr[t][i] == 0.0) {
                } else {
                    inpaintedImageArr[t][i] = inpaintedImageArr[t][i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i]);
                }
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

            //*****************************************************************
//            mserror = sub_mse(refImageArr[t],inpaintedImageArr[t], scatImageArr[t],imgWidth*imgHeight*imgDepth);
//            qDebug() << "MSE error: " << mserror << "\n";
//            aaerror = sub_aae(refImageArr[t],inpaintedImageArr[t], scatImageArr[t], imgWidth*imgHeight*imgDepth);
//            qDebug() << "AAE error: " << aaerror << "\n";

//            int intPart = (int)aaerror;
//            double deciPart = (10*aaerror-10*intPart)/10;
//            string decPartStr = to_string(deciPart);
//            decPartStr = decPartStr.substr(2);
//            string randVoxReconsErrStr = to_string(intPart) + "," + decPartStr;
//            std::cout << randVoxReconsErrStr << std::endl;

//            qDebug() << "Number of iteration: " << stopTime;
            //*****************************************************************
            string randVoxReconsErrStr;
 //            for(int i=0; i<locLen; i++) {
 //                qDebug() << randLocs[i];
 //                randVoxReconsErrStr = randVoxReconsErrStr + " " + to_string(abs(refImageArr[t][int(randLocs[i])]-inpaintedImageArr[t][int(randLocs[i])]));
 //            }

            double number = abs(refImageArr[t][int(randLocs[l])]-inpaintedImageArr[t][int(randLocs[l])]);
            int intPart = (int)number;
            double deciPart = (10*number-10*intPart)/10;
            string decPartStr = to_string(deciPart);
             decPartStr = decPartStr.substr(2);
            randVoxReconsErrStr = to_string(intPart) + "," + decPartStr;
            std::cout << randVoxReconsErrStr << std::endl;
            if(minErr > number) {
                minErr = number;
                minErrIter = stopTime;
            }

            stopTime++;

            l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
//            qDebug() << l2normError << "\n";
        }

        qDebug() << "Optimal Error: " << minErr << "Optimal Error Iter. No: " << minErrIter << "Improvement: " << 100-minErr*100/initErr << "Initial Error: " << initErr;

        qDebug() << "Number of iterations: " << stopTime;

//        mserror = mse(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
//        qDebug() << "MSE error: " << mserror << "\n";
//        aaerror = aae(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
//        qDebug() << "AAE error: " << aaerror << "\n";

        mserror = sub_mse(refImageArr[t],inpaintedImageArr[t], scatImageArr[t],imgWidth*imgHeight*imgDepth);
        qDebug() << "MSE error: " << mserror << "\n";
        aaerror = sub_aae(refImageArr[t],inpaintedImageArr[t], scatImageArr[t], imgWidth*imgHeight*imgDepth);
        qDebug() << "AAE error: " << aaerror << "\n";

        if(t+1<imgTimeLen) {
            l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
//            qDebug() << "l2 norm error: " << l2normError << "\n";
        }
    }

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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}

//4D image signal dropout imputation with EED
void eed_with_qSpace_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int kernelSize = 3, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror=0, aaerror=0, l2normError, alpha;
    double *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double **tempImgArrayPrev, **tempImgArrayCurr;
    double **outConvX, **outConvXY, **outConvXYZ;
    double v1[3], v2[3], v3[3];
    double  *ssdXConv, *ssdYConv, *ssdZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;
    qDebug() << gridSpcX << gridSpcY << gridSpcZ;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double*[sizeof(double)*imgTimeLen];
    tempImgArrayPrev = new double*[sizeof(double)*imgTimeLen];
    for(int i = 0; i < imgTimeLen; ++i) {
        tempImgArrayCurr[i] = new double[sizeof(double)*imgWidth*imgHeight*imgDepth];
        tempImgArrayPrev[i] = new double[sizeof(double)*imgWidth*imgHeight*imgDepth];
    }

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    double contPar = 0.25;
//        if(t == 0) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 21) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 42) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 63) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else if(t == 84) {
//            contPar = 67.5;              // B=0 optimal contrast parameter
//        } else {
//            contPar = 0.0005;            // B=700 or B=2000 optimal contrast parameter
//        }

    l2normError = 1000;
//    l2normError = l2Norm_4D(imageArr,inpaintedImageArr,imgTimeLen,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    int stopTime = 0;
    N = numSteps;
    while(l2normError > tol) {
//            if(stopTime > 3) {
//                break;
//            }
        for(int n=0; n<N; n++) {
            outConvX = convolution3DX_of_4D(inpaintedImageArr, imgTimeLen, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
            outConvXY = convolution3DY_of_4D(outConvX, imgTimeLen, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);             // Calculate the convolution on y axis after x axis.
            outConvXYZ = convolution3DZ_of_4D(outConvXY, imgTimeLen, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);           // Calculate the convolution on z axis after x and y axis.

            imgTimeLen = 23;
            for(int t=22; t<imgTimeLen; t++) {
//                qDebug() << "Volume No: " << t;

                //First convolved derivatives
                ssdXConv = ssd_4DX_central_nonDropouts(outConvXYZ, scatImageArr, imgTimeLen, imgWidth, imgHeight, imgDepth);           // q-space profile similarity of convolved image w.r.t. x.
                ssdYConv = ssd_4DY_central_nonDropouts(outConvXYZ, scatImageArr, imgTimeLen, imgWidth, imgHeight, imgDepth);           // q-space profile similarity of convolved image w.r.t. y.
                ssdZConv = ssd_4DZ_central_nonDropouts(outConvXYZ, scatImageArr, imgTimeLen, imgWidth, imgHeight, imgDepth);           // q-space profile similarity of convolved image w.r.t. z.
//                ssdXConv = ssd_4DX_central(outConvXYZ, imgTimeLen, imgWidth, imgHeight, imgDepth);                            // q-space profile similarity of convolved image w.r.t. x.
//                ssdYConv = ssd_4DY_central(outConvXYZ, imgTimeLen, imgWidth, imgHeight, imgDepth);                            // q-space profile similarity of convolved image w.r.t. y.
//                ssdZConv = ssd_4DZ_central(outConvXYZ, imgTimeLen, imgWidth, imgHeight, imgDepth);                            // q-space profile similarity of convolved image w.r.t. z.

                //First derivatives
                dervX = derivative3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);                         // Derivative w.r.t. x.
                dervY = derivative3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);                         // Derivative w.r.t. y.
                dervZ = derivative3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);                         // Derivative w.r.t. z.

                dervForX = dervForw3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                dervForY = dervForw3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                dervForZ = dervForw3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);
                dervBacX = dervBack3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                dervBacY = dervBack3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                dervBacZ = dervBack3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);

                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(n == 0) {
                        tempImgArrayPrev[t][i] = inpaintedImageArr[t][i];
                    }
                    tempImgArrayCurr[t][i] = inpaintedImageArr[t][i];

                    //Define eigenvectors and entries for the diffusion tensor.
                    double norm_i_square = ssdXConv[i]*ssdXConv[i] + ssdYConv[i]*ssdYConv[i] + ssdZConv[i]*ssdZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
                    double xy_norm_i_square = ssdXConv[i]*ssdXConv[i] + ssdYConv[i]*ssdYConv[i], xy_norm_i = sqrt(xy_norm_i_square);

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
                        v1[2] = ssdZConv[i]/normi;

                        v2[0] = 0;
                        v2[1] = 0;
                        v2[2] = 0;

                        v3[0] = 0;
                        v3[1] = 0;
                        v3[2] = 0;
                    } else {
                        v1[0] = ssdXConv[i]/normi;
                        v1[1] = ssdYConv[i]/normi;
                        v1[2] = ssdZConv[i]/normi;

                        v2[0] = ssdYConv[i]/xy_norm_i;
                        v2[1] = -ssdXConv[i]/xy_norm_i;
                        v2[2] = 0;

                        v3[0] = (ssdXConv[i]*ssdZConv[i])/(xy_norm_i*normi);
                        v3[1] = (ssdYConv[i]*ssdZConv[i])/(xy_norm_i*normi);
                        v3[2] = -xy_norm_i_square/(xy_norm_i*normi);
                    }
                    d1 = charbonnier_diff(norm_i_square, contPar);
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

                dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                alpha = (double)(4*n+2)/(2*n+3);
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    if(scatImageArr[t][i] == 0.0) {
                    } else {
                        inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[t][i];
                        tempImgArrayPrev[t][i] = tempImgArrayCurr[t][i];
                    }
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
                delete [] ssdXConv;
                ssdXConv = NULL;
                delete [] ssdYConv;
                ssdYConv = NULL;
                delete [] ssdZConv;
                ssdZConv = NULL;
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
            }
            //Free each sub-array
            for(int i=0; i < imgTimeLen; ++i) {
                delete[] outConvX[i];
                delete[] outConvXY[i];
                delete[] outConvXYZ[i];
            }
            //Free the array of pointers
            delete [] outConvX;
            outConvX = NULL;
            delete [] outConvXY;
            outConvXY = NULL;
            delete [] outConvXYZ;
            outConvXYZ = NULL;
        }
        stopTime++;

//        l2normError = l2Norm_4D(tempImgArrayCurr, inpaintedImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
        l2normError = l2Norm(tempImgArrayCurr[22], inpaintedImageArr[22], imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";
    }
    mserror = sub_mse_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = sub_aae_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
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

    for(int i=0; i < imgTimeLen; ++i) {
        delete[] tempImgArrayCurr[i];
        delete[] tempImgArrayPrev[i];
    }
    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}


//4D image signal dropout imputation with EED
void eed_shore_inter_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
{
    clock_t begin = clock();

    int kernelSize = 3, N;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror=0, aaerror=0, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double *tempImgArrayPrev, *tempImgArrayCurr;
    double v1[3], v2[3], v3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    gridSpcX = gridSpcY = gridSpcZ = 1;
    qDebug() << gridSpcX << gridSpcY << gridSpcZ;

    //Memory allocation for temporary image arrays.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    vector<int> imgVolLine = corruptedVolumes(scatImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen);
//    vector<int> imgVolLine;
//    imgVolLine.push_back(2);

    for(int i=0; i<imgVolLine.size(); i++) {
        qDebug() << imgVolLine[i] << " ";
    }



    // *********** Begin: Run Python script ************************************************************
    // Initialize the Python Interpreter
    Py_Initialize();

//    QString python_path;
//    python_path = sys.path;
//    python_path<<"sys.path.append(\""<<path<<"/extensions/datastores/\")";
    PyRun_SimpleString("import sys");

    char buff[FILENAME_MAX]; //create string buffer to hold path
    getcwd( buff, FILENAME_MAX);
    string current_working_dir(buff);
    string shore_working_dir = current_working_dir + "/shore_imputation/";
    const char * shore_working_dir_c = shore_working_dir.c_str();
    qDebug() << shore_working_dir_c;
//    std::cout << shore_working_dir;

//    PyRun_SimpleString('shore_working_dir_c');
//    int a = system("python3 shore_working_dir_c");
//    qDebug() << "System output: " << a;

    // Finish the Python Interpreter
    Py_Finalize();
    // *********** End: Run Python script **************************************************************


    double tol_mserror = 0, tol_aaerror = 0;
    for(int v=0; v<imgVolLine.size(); v++) {
        int t = imgVolLine[v];
//        for(int t=22; t<imgTimeLen; t++) {
            qDebug() << "Volume No: " << t;
            //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
            //    qDebug() << "MSE error: " << mserror << "\n";
            //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
            //    printf("Error: %lf\n", aaerror);
//            l2normError = l2Norm(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
            l2normError = 1000;
            qDebug() << "l2 norm error: " << l2normError << "\n";

            double contPar = 0.75;
    //        if(t == 0) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 21) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 42) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 63) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else if(t == 84) {
    //            contPar = 67.5;              // B=0 optimal contrast parameter
    //        } else {
    //            contPar = 0.0005;            // B=700 or B=2000 optimal contrast parameter
    //        }

            int stopTime = 0;
            N = numSteps;
            while(l2normError > tol) {
    //            if(stopTime > 3) {
    //                break;
    //            }
                for(int n=0; n<N; n++) {
                    outConvX = convolution3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
                    outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
                    outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
                    //First convolved derivatives
                    dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
                    dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
                    dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

                    //First derivatives
                    dervX = derivative3DX(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
                    dervY = derivative3DY(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
                    dervZ = derivative3DZ(inpaintedImageArr[t],imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

                    dervForX = dervForw3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                    dervForY = dervForw3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                    dervForZ = dervForw3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);
                    dervBacX = dervBack3DX(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcX);
                    dervBacY = dervBack3DY(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcY);
                    dervBacZ = dervBack3DZ(inpaintedImageArr[t], imgWidth, imgHeight, imgDepth, gridSpcZ);

                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(n == 0) {
                            tempImgArrayPrev[i] = inpaintedImageArr[t][i];
                        }
                        tempImgArrayCurr[i] = inpaintedImageArr[t][i];

                        //Define eigenvectors and entries for the diffusion tensor.
                        double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
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

                        d1 = charbonnier_diff(norm_i_square, contPar);
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

                    dervXD12 = derivative3DX(d12,imgWidth,imgHeight,imgDepth,gridSpcX);
                    dervXD13 = derivative3DX(d13,imgWidth,imgHeight,imgDepth,gridSpcX);
                    dervYD21 = derivative3DY(d21,imgWidth,imgHeight,imgDepth,gridSpcY);
                    dervYD23 = derivative3DY(d23,imgWidth,imgHeight,imgDepth,gridSpcY);
                    dervZD31 = derivative3DZ(d31,imgWidth,imgHeight,imgDepth,gridSpcZ);
                    dervZD32 = derivative3DZ(d32,imgWidth,imgHeight,imgDepth,gridSpcZ);

                    alpha = (double)(4*n+2)/(2*n+3);
                    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                        if(scatImageArr[t][i] == 0.0) {
    //                        randArrTraceIndex++;
    //                        continue;
    //                        qDebug() << scatImageArr[t][i] << imageArr[t][i] << inpaintedImageArr[t][i];
                        } else {
                            inpaintedImageArr[t][i] = alpha*(inpaintedImageArr[t][i] + timeStepSize*((sumForX[i]*dervForX[i] - sumBacX[i]*dervBacX[i] + sumForY[i]*dervForY[i]  - sumBacY[i]*dervBacY[i] + sumForZ[i]*dervForZ[i] - sumBacZ[i]*dervBacZ[i])/2 + dervXD12[i] + dervXD13[i] + dervYD21[i] + dervYD23[i] + dervZD31[i] + dervZD32[i])) + (1-alpha)*tempImgArrayPrev[i];
                            tempImgArrayPrev[i] = tempImgArrayCurr[i];
                        }
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
                }

                stopTime++;

                l2normError = l2Norm(tempImgArrayCurr,inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
                qDebug() << l2normError << "\n";
            }
    //        mserror = mse(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
    //        qDebug() << "MSE error: " << mserror << "\n";
    //        aaerror = aae(refImageArr[t],inpaintedImageArr[t],imgWidth*imgHeight*imgDepth);
    //        qDebug() << "AAE error: " << aaerror << "\n";

            mserror = sub_mse(refImageArr[t],inpaintedImageArr[t], scatImageArr[t],imgWidth*imgHeight*imgDepth);
            qDebug() << "MSE error: " << mserror << "\n";
            tol_mserror = tol_mserror + mserror;
            aaerror = sub_aae(refImageArr[t],inpaintedImageArr[t], scatImageArr[t], imgWidth*imgHeight*imgDepth);
            qDebug() << "AAE error: " << aaerror << "\n";
            tol_aaerror = tol_aaerror + aaerror;

            if(t+1<imgTimeLen) {
                l2normError = l2Norm(imageArr[t+1],inpaintedImageArr[t+1],imgWidth*imgHeight*imgDepth);
                qDebug() << "l2 norm error: " << l2normError << "\n";
            }
//        }
    }
    tol_aaerror = sub_aae_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average AAE error: " << tol_aaerror << "\n";
    tol_mserror = sub_mse_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average MSE error: " << tol_mserror << "\n";


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

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}


//4D image signal dropout imputation with Simple Averaging Neighboring Slices
void average_4d_signalDropoutImputation(double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen)
{
    clock_t begin = clock();

    vector<int> imgVolLine = corruptedVolumes(scatImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen);
    for(int i=0; i<imgVolLine.size(); i++) {
        qDebug() << imgVolLine[i] << " ";
    }

    for(int v=0; v<imgVolLine.size(); v++) {
        int t = imgVolLine[v];
        qDebug() << "Volume No: " << t;

        vector<int> imgSlcLine = corruptedSlices(scatImageArr[t], imgWidth, imgHeight, imgDepth);
     //        for(int i=0; i<imgSlcLine.size(); i++) {
     //            qDebug() << imgSlcLine[i] << " ";
     //        }

        for(int s=0; s<imgSlcLine.size(); s++) {
            if(imgSlcLine[s] == 0) {
                int slc = imgSlcLine[s];
                for(int i=slc*imgWidth*imgHeight; i<(slc+1)*imgWidth*imgHeight; i++) {
                    inpaintedImageArr[t][i] = imageArr[t][i+imgWidth*imgHeight];
                }
            } else if(imgSlcLine[s] != imgDepth) {
                int slc = imgSlcLine[s];
                for(int i=slc*imgWidth*imgHeight; i<(slc+1)*imgWidth*imgHeight; i++) {
                    inpaintedImageArr[t][i] = imageArr[t][i-imgWidth*imgHeight];
                }
            } else {
                int slc = imgSlcLine[s];
                for(int i=slc*imgWidth*imgHeight; i<(slc+1)*imgWidth*imgHeight; i++) {
                    inpaintedImageArr[t][i] = (imageArr[t][i+imgWidth*imgHeight] + imageArr[t][i-imgWidth*imgHeight])/2;
                }
            }
        }
    }

    double tol_aaerror = sub_aae_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average AAE error: " << tol_aaerror << "\n";
    double tol_mserror = sub_mse_4D(refImageArr,inpaintedImageArr, scatImageArr, imgTimeLen, imgWidth*imgHeight*imgDepth);
    qDebug() << "Average MSE error: " << tol_mserror << "\n";

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
//~ Multi Thread Linear Homogenous Fourth Order Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void mulThread_linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, int maskLen)
{
    clock_t begin = clock();

    int randArrTraceIndex;
    float stepSize = 1;
    double *dervXX, *dervYY, *dervZZ, *dervXXXX, *dervYYYY, *dervZZZZ, *dervXXYY, *dervZZXX, *dervYYZZ;

    double mserror, aaerror, l2normError, alpha;
    double* tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    // a global instance of std::mutex to protect global variable
//    std::mutex myMutex;

//    ctpl::thread_pool t(9 /* two threads in the pool */);
//    std::vector<std::future<void>> results(9);

    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << l2normError << "\n";

    while(l2normError > tol) {
        for(int n=0; n<numSteps; n++) {
//            //Second and Fourth derivative
//            std::thread t1([&dervXX, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            std::thread t2([&dervYY, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            std::thread t3([&dervZZ, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            t1.join();
//            std::thread t4([&dervXXXX, dervXX, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervXXXX = derivative3DXX(dervXX,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            std::thread t5([&dervXXYY, dervXX, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervXXYY = derivative3DYY(dervXX,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            t2.join();
//            std::thread t6([&dervYYYY, dervYY, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervYYYY = derivative3DYY(dervYY,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            std::thread t7([&dervYYZZ, dervYY, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervYYZZ = derivative3DZZ(dervYY,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            t3.join();
//            std::thread t8([&dervZZZZ, dervZZ, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervZZZZ = derivative3DZZ(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            std::thread t9([&dervZZXX, dervZZ, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervZZXX = derivative3DXX(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            t4.join();
//            t5.join();
//            t6.join();
//            t7.join();
//            t8.join();
//            t9.join();


//            //Second and Fourth derivative
//            results[0] = t.push([&dervXX, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            results[1] = t.push([&dervYY, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            results[2] = t.push([&dervZZ, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            results[0].get();
//            t.push([&dervXXXX, dervXX, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervXXXX = derivative3DXX(dervXX,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            t.push([&dervXXYY, dervXX, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervXXYY = derivative3DYY(dervXX,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            results[1].get();
//            t.push([&dervYYYY, dervYY, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervYYYY = derivative3DYY(dervYY,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            t.push([&dervYYZZ, dervYY, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervYYZZ = derivative3DZZ(dervYY,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            results[2].get();
//            t.push([&dervZZZZ, dervZZ, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervZZZZ = derivative3DZZ(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);
//            });
//            t.push([&dervZZXX, dervZZ, imgWidth, imgHeight, imgDepth, stepSize](){
//                dervZZXX = derivative3DXX(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);
//            });

//            results[3].get();
//            results[4].get();
//            results[5].get();
//            results[6].get();
//            results[7].get();
//            results[8].get();


            dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. xx.
            dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. yy.
            dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. zz.
            dervXXXX = derivative3DXX(dervXX,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. xxxx.
            dervYYYY = derivative3DYY(dervYY,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. yyyy.
            dervZZZZ = derivative3DZZ(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. zzzz.
            dervXXYY = derivative3DYY(dervXX,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. xxyy.
            dervZZXX = derivative3DXX(dervZZ,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. zzxx.
            dervYYZZ = derivative3DZZ(dervYY,imgWidth,imgHeight,imgDepth,stepSize);                               // Derivative w.r.t. yyzz.

            alpha = (double)(4*n+2)/(2*n+3);

            randArrTraceIndex = 0;
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                if(n == 0) {
                    tempImgArrayPrev[i] = scatImageArr[i];
                }
                tempImgArrayCurr[i] = scatImageArr[i];

                if(i == randPxls[randArrTraceIndex]) {
                    randArrTraceIndex++;
                    continue;
                }
                scatImageArr[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
                tempImgArrayPrev[i] = tempImgArrayCurr[i];
            }

//            qDebug() << maskLen;
//            randArrTraceIndex = 0;
//            std::thread t1([=,&tempImgArrayPrev,&scatImageArr,&tempImgArrayCurr]()mutable{
//                randArrTraceIndex = 0;
//                qDebug() << "Thread 1 (begin): " << randArrTraceIndex << randPxls[randArrTraceIndex];
//                for(int i=0; i<randPxls[maskLen/4]; i++) {
//                    if(n == 0) {
//                        tempImgArrayPrev[i] = scatImageArr[i];
//                    }
//                    tempImgArrayCurr[i] = scatImageArr[i];

//                    if(i == randPxls[randArrTraceIndex]) {
//                        randArrTraceIndex++;
//                        continue;
//                    }
//                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
//                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
//                }
//                qDebug() << "Thread 1 (end): " << randArrTraceIndex << randPxls[randArrTraceIndex];
//            });

//            std::thread t2([=,&tempImgArrayPrev,&scatImageArr,&tempImgArrayCurr]()mutable{
//                randArrTraceIndex = maskLen/4;
//                qDebug() << "Thread 2 (begin): " << randArrTraceIndex << randPxls[randArrTraceIndex];
//                for(int i=randPxls[maskLen/4]; i<imgWidth*imgHeight*imgDepth; i++) {
//                    if(n == 0) {
//                        tempImgArrayPrev[i] = scatImageArr[i];
//                    }
//                    tempImgArrayCurr[i] = scatImageArr[i];

//                    if(i == randPxls[randArrTraceIndex]) {
//                        randArrTraceIndex++;
//                        continue;
//                    }
//                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
//                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
//                }
//                qDebug() << "Thread 2 (end): " << randArrTraceIndex-1 << randPxls[randArrTraceIndex-1];
//            });

//            std::thread t3([=,&randArrTraceIndex,&tempImgArrayPrev,&scatImageArr,&tempImgArrayCurr]()mutable{
//                for(int i=2*(imgWidth*imgHeight*imgDepth)/4; i<3*(imgWidth*imgHeight*imgDepth)/4; i++) {
//                    if(n == 0) {
//                        tempImgArrayPrev[i] = scatImageArr[i];
//                    }
//                    tempImgArrayCurr[i] = scatImageArr[i];

//                    if(i == randPxls[randArrTraceIndex]) {
//                        randArrTraceIndex++;
//                        continue;
//                    }
//                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
//                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
//                }
//            });

//            std::thread t4([=,&randArrTraceIndex,&tempImgArrayPrev,&scatImageArr,&tempImgArrayCurr]()mutable{
//                for(int i=3*(imgWidth*imgHeight*imgDepth)/4; i<imgWidth*imgHeight*imgDepth; i++) {
//                    if(n == 0) {
//                        tempImgArrayPrev[i] = scatImageArr[i];
//                    }
//                    tempImgArrayCurr[i] = scatImageArr[i];

//                    if(i == randPxls[randArrTraceIndex]) {
//                        randArrTraceIndex++;
//                        continue;
//                    }
//                    scatImageArr[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
//                    tempImgArrayPrev[i] = tempImgArrayCurr[i];
//                }
//            });

//            t1.join();
//            qDebug() << randArrTraceIndex;
//            t2.join();
//            qDebug() << randArrTraceIndex;
//            t3.join();
//            qDebug() << randArrTraceIndex;
//            t4.join();
//            qDebug() << randArrTraceIndex;

            delete [] dervXX;
            dervXX = NULL;
            delete [] dervYY;
            dervYY = NULL;
            delete [] dervZZ;
            dervZZ = NULL;
            delete [] dervXXYY;
            dervXXYY = NULL;
            delete dervYYZZ;
            dervYYZZ = NULL;
            delete [] dervZZXX;
            dervZZXX = NULL;
            delete [] dervXXXX;
            dervXXXX = NULL;
            delete [] dervYYYY;
            dervYYYY = NULL;
            delete [] dervZZZZ;
            dervZZZZ = NULL;
        }

        l2normError = l2Norm(scatImageArr,tempImgArrayCurr,imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE: " << aaerror << "\n";

    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;
    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
//~ Multi Thread Fourth Order Edge Enhancing Anisotropic Diffusion based Inpaiting by Fast Semi-Iterative Scheme with fourth order diffusion tensor
void mulThread_foeed_3d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth)
{
    clock_t begin = clock();

    int kernelSize = 3, N, randArrTraceIndex;
    float stepSize = 1;
    double sigma = 1.0, alpha, gausKernel[kernelSize], mserror, aaerror, l2normError;
    double *outConvX, *outConvXY, *outConvXYZ, *dervXX, *dervYY, *dervZZ, *dervXY, *dervYZ, *dervX, *dervY, *dervXZ, *dervXConv, *dervYConv, *dervZConv, *t11, *t12, *t13, *t21, *t22, *t23, *t31, *t32, *t33;
    double *d1111, *d1112, *d1113, *d1121, *d1122, *d1123, *d1131, *d1132, *d1133, *d1211, *d1212, *d1213, *d1221, *d1222, *d1223, *d1231, *d1232, *d1233, *d1311, *d1312, *d1313, *d1321, *d1322, *d1323, *d1331, *d1332, *d1333;
    double *d2111, *d2112, *d2113, *d2121, *d2122, *d2123, *d2131, *d2132, *d2133, *d2211, *d2212, *d2213, *d2221, *d2222, *d2223, *d2231, *d2232, *d2233, *d2311, *d2312, *d2313, *d2321, *d2322, *d2323, *d2331, *d2332, *d2333;
    double *d3111, *d3112, *d3113, *d3121, *d3122, *d3123, *d3131, *d3132, *d3133, *d3211, *d3212, *d3213, *d3221, *d3222, *d3223, *d3231, *d3232, *d3233, *d3311, *d3312, *d3313, *d3321, *d3322, *d3323, *d3331, *d3332, *d3333;
    double *dervXXD11, *dervYYD22, *dervZZD33, *dervYD21, *dervXYD21, *dervXD12, *dervYXD12, *dervZD31, *dervXZD31, *dervXD13, *dervZXD13, *dervYD23, *dervZYD23, *dervZD32, *dervYZD32;
    double *tempImgArrayCurr, *tempImgArrayPrev;

    //************************************************************
    //Memory allocation for 4th order diffusion tensor entries.
    d1111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d1311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d1333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d2311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d2333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3111 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3112 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3113 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3121 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3122 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3123 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3131 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3132 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3133 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3211 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3212 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3213 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3221 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3222 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3223 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3231 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3232 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3233 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    d3311 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3312 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3313 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3321 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3322 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3323 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3331 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3332 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    d3333 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    //************************************************************
    t11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    t33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    //Memory allocation for temporary image array.
    tempImgArrayCurr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    tempImgArrayPrev = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

//    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
//    qDebug() << "MSE error: " << mserror << "\n";
//    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
//    printf("Error: %lf\n", aaerror);
//    qDebug() << "Here";
    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "l2 norm error: " << l2normError << "\n";

    N = numSteps;
    while(l2normError > tol) {
        for(int n=0; n<N; n++) {
            outConvX = convolution3DX(scatImageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);    // Calculate the convolution on x axis.
            outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on y axis after x axis.
            outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);       // Calculate the convolution on z axis after x and y axis.

            //First derivatives of convolved mask
            std::thread t1([&dervXConv, outConvXYZ, imgWidth, imgHeight, imgDepth, stepSize](){
                dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);
            });
            std::thread t2([&dervYConv, outConvXYZ, imgWidth, imgHeight, imgDepth, stepSize](){
                dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);
            });
            std::thread t3([&dervZConv, outConvXYZ, imgWidth, imgHeight, imgDepth, stepSize](){
                dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);
            });

            //First derivatives of convolved mask
//            dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. x.
//            dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. y.
//            dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative of convolved image w.r.t. z.

            //Second derivatives
            std::thread t4([&dervXX, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
                dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
            });
            std::thread t5([&dervZZ, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
                dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
            });
            std::thread t6([&dervYY, scatImageArr, imgWidth, imgHeight, imgDepth, stepSize](){
                dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);
            });

            //Second derivatives
//            dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. xx.
//            dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. yy.
//            dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                         // Derivative w.r.t. zz.
            dervX = derivative3DX(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. x.
            dervXY = derivative3DY(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xy.
            dervXZ = derivative3DZ(dervX,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. xz.
            dervY = derivative3DY(scatImageArr,imgWidth,imgHeight,imgDepth,stepSize);                           // Derivative w.r.t. y.
            dervYZ = derivative3DZ(dervY,imgWidth,imgHeight,imgDepth,stepSize);                                 // Derivative w.r.t. yz.


            t1.join();
            t2.join();
            t3.join();
            t4.join();
            t5.join();
            t6.join();

            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                if(n == 0) {
                    tempImgArrayPrev[i] = scatImageArr[i];
                }
                tempImgArrayCurr[i] = scatImageArr[i];

                //Define eigenvectors and entries for the diffusion tensor.
                double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
                double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);
                double v1[3], v2[3], v3[3], m1, m2, m3, m4, m5, m6;
                double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];

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

//                m1 = green_diff(normi, 1);
//                m1 = aubert_diff(norm_i_square, 1);
//                m1 = li1(normi, 0.2);
//                m1 = pm_diff(norm_i_square, 0.2);
//                m1 = gr_diff(norm_i_square, 1);
//                m1 = pm_diff2(norm_i_square, 0.2);

                m1 = charbonnier_diff(norm_i_square, 0.1);
                m2 = 1;
                m3 = 1;
                m4 = sqrt(m1*m2);
                //m4 = (m1+m2)/2;
                //m4 = m1;
                //m4 = 1;
                m5 = sqrt(m1*m3);
                m6 = 1;

//                for(int i=0; i<3; i++) {
//                    for(int j=0; j<3; j++) {
//                        e1[i][j] = v1[i]*v1[j];
//                    }
//                }
                e1[0][0] = v1[0]*v1[0];
                e1[0][1] = v1[0]*v1[1];
                e1[0][2] = v1[0]*v1[2];
                e1[1][0] = v1[1]*v1[0];
                e1[1][1] = v1[1]*v1[1];
                e1[1][2] = v1[1]*v1[2];
                e1[2][0] = v1[2]*v1[0];
                e1[2][1] = v1[2]*v1[1];
                e1[2][2] = v1[2]*v1[2];

//                for(int i=0; i<3; i++) {
//                    for(int j=0; j<3; j++) {
//                        e2[i][j] = v2[i]*v2[j];
//                    }
//                }
                e2[0][0] = v2[0]*v2[0];
                e2[0][1] = v2[0]*v2[1];
                e2[0][2] = v2[0]*v2[2];
                e2[1][0] = v2[1]*v2[0];
                e2[1][1] = v2[1]*v2[1];
                e2[1][2] = v2[1]*v2[2];
                e2[2][0] = v2[2]*v2[0];
                e2[2][1] = v2[2]*v2[1];
                e2[2][2] = v2[2]*v2[2];

//                for(int i=0; i<3; i++) {
//                    for(int j=0; j<3; j++) {
//                        e3[i][j] = v3[i]*v3[j];
//                    }
//                }
                e3[0][0] = v3[0]*v3[0];
                e3[0][1] = v3[0]*v3[1];
                e3[0][2] = v3[0]*v3[2];
                e3[1][0] = v3[1]*v3[0];
                e3[1][1] = v3[1]*v3[1];
                e3[1][2] = v3[1]*v3[2];
                e3[2][0] = v3[2]*v3[0];
                e3[2][1] = v3[2]*v3[1];
                e3[2][2] = v3[2]*v3[2];

                e4[0][0] = (sqrt(2)*v1[0]*v2[0]);
                e4[0][1] = (v1[1]*v2[0] + v2[1]*v1[0])/sqrt(2);
                e4[0][2] = (v1[2]*v2[0] + v2[2]*v1[0])/sqrt(2);
                e4[1][0] = (v1[0]*v2[1] + v2[0]*v1[1])/sqrt(2);
                e4[1][1] = (sqrt(2)*v1[1]*v2[1]);
                e4[1][2] = (v1[2]*v2[1] + v2[2]*v1[1])/sqrt(2);
                e4[2][0] = (v1[0]*v2[2] + v2[0]*v1[2])/sqrt(2);
                e4[2][1] = (v1[1]*v2[2] + v2[1]*v1[2])/sqrt(2);
                e4[2][2] = (sqrt(2)*v1[2]*v2[2]);

                e5[0][0] = (sqrt(2)*v1[0]*v3[0]);
                e5[0][1] = (v1[1]*v3[0] + v3[1]*v1[0])/sqrt(2);
                e5[0][2] = (v1[2]*v3[0] + v3[2]*v1[0])/sqrt(2);
                e5[1][0] = (v1[0]*v3[1] + v3[0]*v1[1])/sqrt(2);
                e5[1][1] = (sqrt(2)*v1[1]*v3[1]);
                e5[1][2] = (v1[2]*v3[1] + v3[2]*v1[1])/sqrt(2);
                e5[2][0] = (v1[0]*v3[2] + v3[0]*v1[2])/sqrt(2);
                e5[2][1] = (v1[1]*v3[2] + v3[1]*v1[2])/sqrt(2);
                e5[2][2] = (sqrt(2)*v1[2]*v3[2]);

                e6[0][0] = (sqrt(2)*v2[0]*v3[0]);
                e6[0][1] = (v2[1]*v3[0] + v3[1]*v2[0])/sqrt(2);
                e6[0][2] = (v2[2]*v3[0] + v3[2]*v2[0])/sqrt(2);
                e6[1][0] = (v2[0]*v3[1] + v3[0]*v2[1])/sqrt(2);
                e6[1][1] = (sqrt(2)*v2[1]*v3[1]);
                e6[1][2] = (v2[2]*v3[1] + v3[2]*v2[1])/sqrt(2);
                e6[2][0] = (v2[0]*v3[2] + v3[0]*v2[2])/sqrt(2);
                e6[2][1] = (v2[1]*v3[2] + v3[1]*v2[2])/sqrt(2);
                e6[2][2] = (sqrt(2)*v2[2]*v3[2]);

                d1111[i] = m1*e1[0][0]*e1[0][0] + m2*e2[0][0]*e2[0][0] + m3*e3[0][0]*e3[0][0] + m4*e4[0][0]*e4[0][0] + m5*e5[0][0]*e5[0][0] + m6*e6[0][0]*e6[0][0];
                d1112[i] = m1*e1[0][0]*e1[0][1] + m2*e2[0][0]*e2[0][1] + m3*e3[0][0]*e3[0][1] + m4*e4[0][0]*e4[0][1] + m5*e5[0][0]*e5[0][1] + m6*e6[0][0]*e6[0][1];
                d1113[i] = m1*e1[0][0]*e1[0][2] + m2*e2[0][0]*e2[0][2] + m3*e3[0][0]*e3[0][2] + m4*e4[0][0]*e4[0][2] + m5*e5[0][0]*e5[0][2] + m6*e6[0][0]*e6[0][2];
                d1121[i] = m1*e1[0][0]*e1[1][0] + m2*e2[0][0]*e2[1][0] + m3*e3[0][0]*e3[1][0] + m4*e4[0][0]*e4[1][0] + m5*e5[0][0]*e5[1][0] + m6*e6[0][0]*e6[1][0];
                d1122[i] = m1*e1[0][0]*e1[1][1] + m2*e2[0][0]*e2[1][1] + m3*e3[0][0]*e3[1][1] + m4*e4[0][0]*e4[1][1] + m5*e5[0][0]*e5[1][1] + m6*e6[0][0]*e6[1][1];
                d1123[i] = m1*e1[0][0]*e1[1][2] + m2*e2[0][0]*e2[1][2] + m3*e3[0][0]*e3[1][2] + m4*e4[0][0]*e4[1][2] + m5*e5[0][0]*e5[1][2] + m6*e6[0][0]*e6[1][2];
                d1131[i] = m1*e1[0][0]*e1[2][0] + m2*e2[0][0]*e2[2][0] + m3*e3[0][0]*e3[2][0] + m4*e4[0][0]*e4[2][0] + m5*e5[0][0]*e5[2][0] + m6*e6[0][0]*e6[2][0];
                d1132[i] = m1*e1[0][0]*e1[2][1] + m2*e2[0][0]*e2[2][1] + m3*e3[0][0]*e3[2][1] + m4*e4[0][0]*e4[2][1] + m5*e5[0][0]*e5[2][1] + m6*e6[0][0]*e6[2][1];
                d1133[i] = m1*e1[0][0]*e1[2][2] + m2*e2[0][0]*e2[2][2] + m3*e3[0][0]*e3[2][2] + m4*e4[0][0]*e4[2][2] + m5*e5[0][0]*e5[2][2] + m6*e6[0][0]*e6[2][2];

                d1211[i] = m1*e1[0][1]*e1[0][0] + m2*e2[0][1]*e2[0][0] + m3*e3[0][1]*e3[0][0] + m4*e4[0][1]*e4[0][0] + m5*e5[0][1]*e5[0][0] + m6*e6[0][1]*e6[0][0];
                d1212[i] = m1*e1[0][1]*e1[0][1] + m2*e2[0][1]*e2[0][1] + m3*e3[0][1]*e3[0][1] + m4*e4[0][1]*e4[0][1] + m5*e5[0][1]*e5[0][1] + m6*e6[0][1]*e6[0][1];
                d1213[i] = m1*e1[0][1]*e1[0][2] + m2*e2[0][1]*e2[0][2] + m3*e3[0][1]*e3[0][2] + m4*e4[0][1]*e4[0][2] + m5*e5[0][1]*e5[0][2] + m6*e6[0][1]*e6[0][2];
                d1221[i] = m1*e1[0][1]*e1[1][0] + m2*e2[0][1]*e2[1][0] + m3*e3[0][1]*e3[1][0] + m4*e4[0][1]*e4[1][0] + m5*e5[0][1]*e5[1][0] + m6*e6[0][1]*e6[1][0];
                d1222[i] = m1*e1[0][1]*e1[1][1] + m2*e2[0][1]*e2[1][1] + m3*e3[0][1]*e3[1][1] + m4*e4[0][1]*e4[1][1] + m5*e5[0][1]*e5[1][1] + m6*e6[0][1]*e6[1][1];
                d1223[i] = m1*e1[0][1]*e1[1][2] + m2*e2[0][1]*e2[1][2] + m3*e3[0][1]*e3[1][2] + m4*e4[0][1]*e4[1][2] + m5*e5[0][1]*e5[1][2] + m6*e6[0][1]*e6[1][2];
                d1231[i] = m1*e1[0][1]*e1[2][0] + m2*e2[0][1]*e2[2][0] + m3*e3[0][1]*e3[2][0] + m4*e4[0][1]*e4[2][0] + m5*e5[0][1]*e5[2][0] + m6*e6[0][1]*e6[2][0];
                d1232[i] = m1*e1[0][1]*e1[2][1] + m2*e2[0][1]*e2[2][1] + m3*e3[0][1]*e3[2][1] + m4*e4[0][1]*e4[2][1] + m5*e5[0][1]*e5[2][1] + m6*e6[0][1]*e6[2][1];
                d1233[i] = m1*e1[0][1]*e1[2][2] + m2*e2[0][1]*e2[2][2] + m3*e3[0][1]*e3[2][2] + m4*e4[0][1]*e4[2][2] + m5*e5[0][1]*e5[2][2] + m6*e6[0][1]*e6[2][2];

                d1311[i] = m1*e1[0][2]*e1[0][0] + m2*e2[0][2]*e2[0][0] + m3*e3[0][2]*e3[0][0] + m4*e4[0][2]*e4[0][0] + m5*e5[0][2]*e5[0][0] + m6*e6[0][2]*e6[0][0];
                d1312[i] = m1*e1[0][2]*e1[0][1] + m2*e2[0][2]*e2[0][1] + m3*e3[0][2]*e3[0][1] + m4*e4[0][2]*e4[0][1] + m5*e5[0][2]*e5[0][1] + m6*e6[0][2]*e6[0][1];
                d1313[i] = m1*e1[0][2]*e1[0][2] + m2*e2[0][2]*e2[0][2] + m3*e3[0][2]*e3[0][2] + m4*e4[0][2]*e4[0][2] + m5*e5[0][2]*e5[0][2] + m6*e6[0][2]*e6[0][2];
                d1321[i] = m1*e1[0][2]*e1[1][0] + m2*e2[0][2]*e2[1][0] + m3*e3[0][2]*e3[1][0] + m4*e4[0][2]*e4[1][0] + m5*e5[0][2]*e5[1][0] + m6*e6[0][2]*e6[1][0];
                d1322[i] = m1*e1[0][2]*e1[1][1] + m2*e2[0][2]*e2[1][1] + m3*e3[0][2]*e3[1][1] + m4*e4[0][2]*e4[1][1] + m5*e5[0][2]*e5[1][1] + m6*e6[0][2]*e6[1][1];
                d1323[i] = m1*e1[0][2]*e1[1][2] + m2*e2[0][2]*e2[1][2] + m3*e3[0][2]*e3[1][2] + m4*e4[0][2]*e4[1][2] + m5*e5[0][2]*e5[1][2] + m6*e6[0][2]*e6[1][2];
                d1331[i] = m1*e1[0][2]*e1[2][0] + m2*e2[0][2]*e2[2][0] + m3*e3[0][2]*e3[2][0] + m4*e4[0][2]*e4[2][0] + m5*e5[0][2]*e5[2][0] + m6*e6[0][2]*e6[2][0];
                d1332[i] = m1*e1[0][2]*e1[2][1] + m2*e2[0][2]*e2[2][1] + m3*e3[0][2]*e3[2][1] + m4*e4[0][2]*e4[2][1] + m5*e5[0][2]*e5[2][1] + m6*e6[0][2]*e6[2][1];
                d1333[i] = m1*e1[0][2]*e1[2][2] + m2*e2[0][2]*e2[2][2] + m3*e3[0][2]*e3[2][2] + m4*e4[0][2]*e4[2][2] + m5*e5[0][2]*e5[2][2] + m6*e6[0][2]*e6[2][2];

                d2111[i] = m1*e1[1][0]*e1[0][0] + m2*e2[1][0]*e2[0][0] + m3*e3[1][0]*e3[0][0] + m4*e4[1][0]*e4[0][0] + m5*e5[1][0]*e5[0][0] + m6*e6[1][0]*e6[0][0];
                d2112[i] = m1*e1[1][0]*e1[0][1] + m2*e2[1][0]*e2[0][1] + m3*e3[1][0]*e3[0][1] + m4*e4[1][0]*e4[0][1] + m5*e5[1][0]*e5[0][1] + m6*e6[1][0]*e6[0][1];
                d2113[i] = m1*e1[1][0]*e1[0][2] + m2*e2[1][0]*e2[0][2] + m3*e3[1][0]*e3[0][2] + m4*e4[1][0]*e4[0][2] + m5*e5[1][0]*e5[0][2] + m6*e6[1][0]*e6[0][2];
                d2121[i] = m1*e1[1][0]*e1[1][0] + m2*e2[1][0]*e2[1][0] + m3*e3[1][0]*e3[1][0] + m4*e4[1][0]*e4[1][0] + m5*e5[1][0]*e5[1][0] + m6*e6[1][0]*e6[1][0];
                d2122[i] = m1*e1[1][0]*e1[1][1] + m2*e2[1][0]*e2[1][1] + m3*e3[1][0]*e3[1][1] + m4*e4[1][0]*e4[1][1] + m5*e5[1][0]*e5[1][1] + m6*e6[1][0]*e6[1][1];
                d2123[i] = m1*e1[1][0]*e1[1][2] + m2*e2[1][0]*e2[1][2] + m3*e3[1][0]*e3[1][2] + m4*e4[1][0]*e4[1][2] + m5*e5[1][0]*e5[1][2] + m6*e6[1][0]*e6[1][2];
                d2131[i] = m1*e1[1][0]*e1[2][0] + m2*e2[1][0]*e2[2][0] + m3*e3[1][0]*e3[2][0] + m4*e4[1][0]*e4[2][0] + m5*e5[1][0]*e5[2][0] + m6*e6[1][0]*e6[2][0];
                d2132[i] = m1*e1[1][0]*e1[2][1] + m2*e2[1][0]*e2[2][1] + m3*e3[1][0]*e3[2][1] + m4*e4[1][0]*e4[2][1] + m5*e5[1][0]*e5[2][1] + m6*e6[1][0]*e6[2][1];
                d2133[i] = m1*e1[1][0]*e1[2][2] + m2*e2[1][0]*e2[2][2] + m3*e3[1][0]*e3[2][2] + m4*e4[1][0]*e4[2][2] + m5*e5[1][0]*e5[2][2] + m6*e6[1][0]*e6[2][2];

                d2211[i] = m1*e1[1][1]*e1[0][0] + m2*e2[1][1]*e2[0][0] + m3*e3[1][1]*e3[0][0] + m4*e4[1][1]*e4[0][0] + m5*e5[1][1]*e5[0][0] + m6*e6[1][1]*e6[0][0];
                d2212[i] = m1*e1[1][1]*e1[0][1] + m2*e2[1][1]*e2[0][1] + m3*e3[1][1]*e3[0][1] + m4*e4[1][1]*e4[0][1] + m5*e5[1][1]*e5[0][1] + m6*e6[1][1]*e6[0][1];
                d2213[i] = m1*e1[1][1]*e1[0][2] + m2*e2[1][1]*e2[0][2] + m3*e3[1][1]*e3[0][2] + m4*e4[1][1]*e4[0][2] + m5*e5[1][1]*e5[0][2] + m6*e6[1][1]*e6[0][2];
                d2221[i] = m1*e1[1][1]*e1[1][0] + m2*e2[1][1]*e2[1][0] + m3*e3[1][1]*e3[1][0] + m4*e4[1][1]*e4[1][0] + m5*e5[1][1]*e5[1][0] + m6*e6[1][1]*e6[1][0];
                d2222[i] = m1*e1[1][1]*e1[1][1] + m2*e2[1][1]*e2[1][1] + m3*e3[1][1]*e3[1][1] + m4*e4[1][1]*e4[1][1] + m5*e5[1][1]*e5[1][1] + m6*e6[1][1]*e6[1][1];
                d2223[i] = m1*e1[1][1]*e1[1][2] + m2*e2[1][1]*e2[1][2] + m3*e3[1][1]*e3[1][2] + m4*e4[1][1]*e4[1][2] + m5*e5[1][1]*e5[1][2] + m6*e6[1][1]*e6[1][2];
                d2231[i] = m1*e1[1][1]*e1[2][0] + m2*e2[1][1]*e2[2][0] + m3*e3[1][1]*e3[2][0] + m4*e4[1][1]*e4[2][0] + m5*e5[1][1]*e5[2][0] + m6*e6[1][1]*e6[2][0];
                d2232[i] = m1*e1[1][1]*e1[2][1] + m2*e2[1][1]*e2[2][1] + m3*e3[1][1]*e3[2][1] + m4*e4[1][1]*e4[2][1] + m5*e5[1][1]*e5[2][1] + m6*e6[1][1]*e6[2][1];
                d2233[i] = m1*e1[1][1]*e1[2][2] + m2*e2[1][1]*e2[2][2] + m3*e3[1][1]*e3[2][2] + m4*e4[1][1]*e4[2][2] + m5*e5[1][1]*e5[2][2] + m6*e6[1][1]*e6[2][2];

                d2311[i] = m1*e1[1][2]*e1[0][0] + m2*e2[1][2]*e2[0][0] + m3*e3[1][2]*e3[0][0] + m4*e4[1][2]*e4[0][0] + m5*e5[1][2]*e5[0][0] + m6*e6[1][2]*e6[0][0];
                d2312[i] = m1*e1[1][2]*e1[0][1] + m2*e2[1][2]*e2[0][1] + m3*e3[1][2]*e3[0][1] + m4*e4[1][2]*e4[0][1] + m5*e5[1][2]*e5[0][1] + m6*e6[1][2]*e6[0][1];
                d2313[i] = m1*e1[1][2]*e1[0][2] + m2*e2[1][2]*e2[0][2] + m3*e3[1][2]*e3[0][2] + m4*e4[1][2]*e4[0][2] + m5*e5[1][2]*e5[0][2] + m6*e6[1][2]*e6[0][2];
                d2321[i] = m1*e1[1][2]*e1[1][0] + m2*e2[1][2]*e2[1][0] + m3*e3[1][2]*e3[1][0] + m4*e4[1][2]*e4[1][0] + m5*e5[1][2]*e5[1][0] + m6*e6[1][2]*e6[1][0];
                d2322[i] = m1*e1[1][2]*e1[1][1] + m2*e2[1][2]*e2[1][1] + m3*e3[1][2]*e3[1][1] + m4*e4[1][2]*e4[1][1] + m5*e5[1][2]*e5[1][1] + m6*e6[1][2]*e6[1][1];
                d2323[i] = m1*e1[1][2]*e1[1][2] + m2*e2[1][2]*e2[1][2] + m3*e3[1][2]*e3[1][2] + m4*e4[1][2]*e4[1][2] + m5*e5[1][2]*e5[1][2] + m6*e6[1][2]*e6[1][2];
                d2331[i] = m1*e1[1][2]*e1[2][0] + m2*e2[1][2]*e2[2][0] + m3*e3[1][2]*e3[2][0] + m4*e4[1][2]*e4[2][0] + m5*e5[1][2]*e5[2][0] + m6*e6[1][2]*e6[2][0];
                d2332[i] = m1*e1[1][2]*e1[2][1] + m2*e2[1][2]*e2[2][1] + m3*e3[1][2]*e3[2][1] + m4*e4[1][2]*e4[2][1] + m5*e5[1][2]*e5[2][1] + m6*e6[1][2]*e6[2][1];
                d2333[i] = m1*e1[1][2]*e1[2][2] + m2*e2[1][2]*e2[2][2] + m3*e3[1][2]*e3[2][2] + m4*e4[1][2]*e4[2][2] + m5*e5[1][2]*e5[2][2] + m6*e6[1][2]*e6[2][2];


                d3111[i] = m1*e1[2][0]*e1[0][0] + m2*e2[2][0]*e2[0][0] + m3*e3[2][0]*e3[0][0] + m4*e4[2][0]*e4[0][0] + m5*e5[2][0]*e5[0][0] + m6*e6[2][0]*e6[0][0];
                d3112[i] = m1*e1[2][0]*e1[0][1] + m2*e2[2][0]*e2[0][1] + m3*e3[2][0]*e3[0][1] + m4*e4[2][0]*e4[0][1] + m5*e5[2][0]*e5[0][1] + m6*e6[2][0]*e6[0][1];
                d3113[i] = m1*e1[2][0]*e1[0][2] + m2*e2[2][0]*e2[0][2] + m3*e3[2][0]*e3[0][2] + m4*e4[2][0]*e4[0][2] + m5*e5[2][0]*e5[0][2] + m6*e6[2][0]*e6[0][2];
                d3121[i] = m1*e1[2][0]*e1[1][0] + m2*e2[2][0]*e2[1][0] + m3*e3[2][0]*e3[1][0] + m4*e4[2][0]*e4[1][0] + m5*e5[2][0]*e5[1][0] + m6*e6[2][0]*e6[1][0];
                d3122[i] = m1*e1[2][0]*e1[1][1] + m2*e2[2][0]*e2[1][1] + m3*e3[2][0]*e3[1][1] + m4*e4[2][0]*e4[1][1] + m5*e5[2][0]*e5[1][1] + m6*e6[2][0]*e6[1][1];
                d3123[i] = m1*e1[2][0]*e1[1][2] + m2*e2[2][0]*e2[1][2] + m3*e3[2][0]*e3[1][2] + m4*e4[2][0]*e4[1][2] + m5*e5[2][0]*e5[1][2] + m6*e6[2][0]*e6[1][2];
                d3131[i] = m1*e1[2][0]*e1[2][0] + m2*e2[2][0]*e2[2][0] + m3*e3[2][0]*e3[2][0] + m4*e4[2][0]*e4[2][0] + m5*e5[2][0]*e5[2][0] + m6*e6[2][0]*e6[2][0];
                d3132[i] = m1*e1[2][0]*e1[2][1] + m2*e2[2][0]*e2[2][1] + m3*e3[2][0]*e3[2][1] + m4*e4[2][0]*e4[2][1] + m5*e5[2][0]*e5[2][1] + m6*e6[2][0]*e6[2][1];
                d3133[i] = m1*e1[2][0]*e1[2][2] + m2*e2[2][0]*e2[2][2] + m3*e3[2][0]*e3[2][2] + m4*e4[2][0]*e4[2][2] + m5*e5[2][0]*e5[2][2] + m6*e6[2][0]*e6[2][2];

                d3211[i] = m1*e1[2][1]*e1[0][0] + m2*e2[2][1]*e2[0][0] + m3*e3[2][1]*e3[0][0] + m4*e4[2][1]*e4[0][0] + m5*e5[2][1]*e5[0][0] + m6*e6[2][1]*e6[0][0];
                d3212[i] = m1*e1[2][1]*e1[0][1] + m2*e2[2][1]*e2[0][1] + m3*e3[2][1]*e3[0][1] + m4*e4[2][1]*e4[0][1] + m5*e5[2][1]*e5[0][1] + m6*e6[2][1]*e6[0][1];
                d3213[i] = m1*e1[2][1]*e1[0][2] + m2*e2[2][1]*e2[0][2] + m3*e3[2][1]*e3[0][2] + m4*e4[2][1]*e4[0][2] + m5*e5[2][1]*e5[0][2] + m6*e6[2][1]*e6[0][2];
                d3221[i] = m1*e1[2][1]*e1[1][0] + m2*e2[2][1]*e2[1][0] + m3*e3[2][1]*e3[1][0] + m4*e4[2][1]*e4[1][0] + m5*e5[2][1]*e5[1][0] + m6*e6[2][1]*e6[1][0];
                d3222[i] = m1*e1[2][1]*e1[1][1] + m2*e2[2][1]*e2[1][1] + m3*e3[2][1]*e3[1][1] + m4*e4[2][1]*e4[1][1] + m5*e5[2][1]*e5[1][1] + m6*e6[2][1]*e6[1][1];
                d3223[i] = m1*e1[2][1]*e1[1][2] + m2*e2[2][1]*e2[1][2] + m3*e3[2][1]*e3[1][2] + m4*e4[2][1]*e4[1][2] + m5*e5[2][1]*e5[1][2] + m6*e6[2][1]*e6[1][2];
                d3231[i] = m1*e1[2][1]*e1[2][0] + m2*e2[2][1]*e2[2][0] + m3*e3[2][1]*e3[2][0] + m4*e4[2][1]*e4[2][0] + m5*e5[2][1]*e5[2][0] + m6*e6[2][1]*e6[2][0];
                d3232[i] = m1*e1[2][1]*e1[2][1] + m2*e2[2][1]*e2[2][1] + m3*e3[2][1]*e3[2][1] + m4*e4[2][1]*e4[2][1] + m5*e5[2][1]*e5[2][1] + m6*e6[2][1]*e6[2][1];
                d3233[i] = m1*e1[2][1]*e1[2][2] + m2*e2[2][1]*e2[2][2] + m3*e3[2][1]*e3[2][2] + m4*e4[2][1]*e4[2][2] + m5*e5[2][1]*e5[2][2] + m6*e6[2][1]*e6[2][2];

                d3311[i] = m1*e1[2][2]*e1[0][0] + m2*e2[2][2]*e2[0][0] + m3*e3[2][2]*e3[0][0] + m4*e4[2][2]*e4[0][0] + m5*e5[2][2]*e5[0][0] + m6*e6[2][2]*e6[0][0];
                d3312[i] = m1*e1[2][2]*e1[0][1] + m2*e2[2][2]*e2[0][1] + m3*e3[2][2]*e3[0][1] + m4*e4[2][2]*e4[0][1] + m5*e5[2][2]*e5[0][1] + m6*e6[2][2]*e6[0][1];
                d3313[i] = m1*e1[2][2]*e1[0][2] + m2*e2[2][2]*e2[0][2] + m3*e3[2][2]*e3[0][2] + m4*e4[2][2]*e4[0][2] + m5*e5[2][2]*e5[0][2] + m6*e6[2][2]*e6[0][2];
                d3321[i] = m1*e1[2][2]*e1[1][0] + m2*e2[2][2]*e2[1][0] + m3*e3[2][2]*e3[1][0] + m4*e4[2][2]*e4[1][0] + m5*e5[2][2]*e5[1][0] + m6*e6[2][2]*e6[1][0];
                d3322[i] = m1*e1[2][2]*e1[1][1] + m2*e2[2][2]*e2[1][1] + m3*e3[2][2]*e3[1][1] + m4*e4[2][2]*e4[1][1] + m5*e5[2][2]*e5[1][1] + m6*e6[2][2]*e6[1][1];
                d3323[i] = m1*e1[2][2]*e1[1][2] + m2*e2[2][2]*e2[1][2] + m3*e3[2][2]*e3[1][2] + m4*e4[2][2]*e4[1][2] + m5*e5[2][2]*e5[1][2] + m6*e6[2][2]*e6[1][2];
                d3331[i] = m1*e1[2][2]*e1[2][0] + m2*e2[2][2]*e2[2][0] + m3*e3[2][2]*e3[2][0] + m4*e4[2][2]*e4[2][0] + m5*e5[2][2]*e5[2][0] + m6*e6[2][2]*e6[2][0];
                d3332[i] = m1*e1[2][2]*e1[2][1] + m2*e2[2][2]*e2[2][1] + m3*e3[2][2]*e3[2][1] + m4*e4[2][2]*e4[2][1] + m5*e5[2][2]*e5[2][1] + m6*e6[2][2]*e6[2][1];
                d3333[i] = m1*e1[2][2]*e1[2][2] + m2*e2[2][2]*e2[2][2] + m3*e3[2][2]*e3[2][2] + m4*e4[2][2]*e4[2][2] + m5*e5[2][2]*e5[2][2] + m6*e6[2][2]*e6[2][2];

                t11[i] = d1111[i]*dervXX[i] + d1112[i]*dervXY[i] + d1113[i]*dervXZ[i] + d1121[i]*dervXY[i] + d1122[i]*dervYY[i] + d1123[i]*dervYZ[i] + d1131[i]*dervXZ[i] + d1132[i]*dervYZ[i] + d1133[i]*dervZZ[i];
                t12[i] = d1211[i]*dervXX[i] + d1212[i]*dervXY[i] + d1213[i]*dervXZ[i] + d1221[i]*dervXY[i] + d1222[i]*dervYY[i] + d1223[i]*dervYZ[i] + d1231[i]*dervXZ[i] + d1232[i]*dervYZ[i] + d1233[i]*dervZZ[i];
                t13[i] = d1311[i]*dervXX[i] + d1312[i]*dervXY[i] + d1313[i]*dervXZ[i] + d1321[i]*dervXY[i] + d1322[i]*dervYY[i] + d1323[i]*dervYZ[i] + d1331[i]*dervXZ[i] + d1332[i]*dervYZ[i] + d1333[i]*dervZZ[i];

                t21[i] = d2111[i]*dervXX[i] + d2112[i]*dervXY[i] + d2113[i]*dervXZ[i] + d2121[i]*dervXY[i] + d2122[i]*dervYY[i] + d2123[i]*dervYZ[i] + d2131[i]*dervXZ[i] + d2132[i]*dervYZ[i] + d2133[i]*dervZZ[i];
                t22[i] = d2211[i]*dervXX[i] + d2212[i]*dervXY[i] + d2213[i]*dervXZ[i] + d2221[i]*dervXY[i] + d2222[i]*dervYY[i] + d2223[i]*dervYZ[i] + d2231[i]*dervXZ[i] + d2232[i]*dervYZ[i] + d2233[i]*dervZZ[i];
                t23[i] = d2311[i]*dervXX[i] + d2312[i]*dervXY[i] + d2313[i]*dervXZ[i] + d2321[i]*dervXY[i] + d2322[i]*dervYY[i] + d2323[i]*dervYZ[i] + d2331[i]*dervXZ[i] + d2332[i]*dervYZ[i] + d2333[i]*dervZZ[i];

                t31[i] = d3111[i]*dervXX[i] + d3112[i]*dervXY[i] + d3113[i]*dervXZ[i] + d3121[i]*dervXY[i] + d3122[i]*dervYY[i] + d3123[i]*dervYZ[i] + d3131[i]*dervXZ[i] + d3132[i]*dervYZ[i] + d3133[i]*dervZZ[i];
                t32[i] = d3211[i]*dervXX[i] + d3212[i]*dervXY[i] + d3213[i]*dervXZ[i] + d3221[i]*dervXY[i] + d3222[i]*dervYY[i] + d3223[i]*dervYZ[i] + d3231[i]*dervXZ[i] + d3232[i]*dervYZ[i] + d3233[i]*dervZZ[i];
                t33[i] = d3311[i]*dervXX[i] + d3312[i]*dervXY[i] + d3313[i]*dervXZ[i] + d3321[i]*dervXY[i] + d3322[i]*dervYY[i] + d3323[i]*dervYZ[i] + d3331[i]*dervXZ[i] + d3332[i]*dervYZ[i] + d3333[i]*dervZZ[i];
            }

            std::thread t7([&dervXXD11, t11, imgWidth, imgHeight, imgDepth, stepSize](){
                dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
            });
            std::thread t8([&dervYYD22, t22, imgWidth, imgHeight, imgDepth, stepSize](){
                dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
            });
            std::thread t9([&dervZZD33, t33, imgWidth, imgHeight, imgDepth, stepSize](){
                dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
            });

//            dervXXD11 = derivative3DXX(t11,imgWidth,imgHeight,imgDepth,stepSize);
//            dervYYD22 = derivative3DYY(t22,imgWidth,imgHeight,imgDepth,stepSize);
//            dervZZD33 = derivative3DZZ(t33,imgWidth,imgHeight,imgDepth,stepSize);
            dervYD21 = derivative3DY(t21,imgWidth,imgHeight,imgDepth,stepSize);
            dervXYD21 = derivative3DX(dervYD21,imgWidth,imgHeight,imgDepth,stepSize);
            dervZD31 = derivative3DZ(t31,imgWidth,imgHeight,imgDepth,stepSize);
            dervXZD31 = derivative3DX(dervZD31,imgWidth,imgHeight,imgDepth,stepSize);
            dervXD12 = derivative3DX(t12,imgWidth,imgHeight,imgDepth,stepSize);
            dervYXD12 = derivative3DY(dervXD12,imgWidth,imgHeight,imgDepth,stepSize);
            dervZD32 = derivative3DZ(t32,imgWidth,imgHeight,imgDepth,stepSize);
            dervYZD32 = derivative3DY(dervZD32,imgWidth,imgHeight,imgDepth,stepSize);
            dervXD13 = derivative3DX(t13,imgWidth,imgHeight,imgDepth,stepSize);
            dervZXD13 = derivative3DZ(dervXD13,imgWidth,imgHeight,imgDepth,stepSize);
            dervYD23 = derivative3DY(t23,imgWidth,imgHeight,imgDepth,stepSize);
            dervZYD23 = derivative3DZ(dervYD23,imgWidth,imgHeight,imgDepth,stepSize);

            t7.join();
            t8.join();
            t9.join();

            alpha = (double)(4*n+2)/(2*n+3);

            randArrTraceIndex = 0;
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                if(i == randPxls[randArrTraceIndex]) {
                    randArrTraceIndex++;
                    continue;
                }
                scatImageArr[i] = alpha*(scatImageArr[i] - timeStep*(dervXXD11[i] + dervXYD21[i] + dervXZD31[i] + dervYXD12[i] + dervYYD22[i] + dervYZD32[i] + dervZXD13[i] + dervZYD23[i] + dervZZD33[i])) + (1-alpha)*tempImgArrayPrev[i];
                tempImgArrayPrev[i] = tempImgArrayCurr[i];
            }
            delete [] dervX;
            dervX = NULL;
            delete [] dervY;
            dervY = NULL;
            delete [] dervXX;
            dervXX = NULL;
            delete [] dervYY;
            dervYY = NULL;
            delete [] dervZZ;
            dervZZ = NULL;
            delete [] dervXY;
            dervXY = NULL;
            delete [] dervXZ;
            dervXZ = NULL;
            delete [] dervYZ;
            dervYZ = NULL;
            delete [] dervZZD33;
            dervZZD33 = NULL;
            delete [] dervXXD11;
            dervXXD11 = NULL;
            delete [] dervXYD21;
            dervXYD21 = NULL;
            delete [] dervYXD12;
            dervYXD12 = NULL;
            delete [] dervYD21;
            dervYD21 = NULL;
            delete [] dervZD31;
            dervZD31 = NULL;
            delete [] dervYYD22;
            dervYYD22 = NULL;
            delete [] dervXZD31;
            dervXZD31 = NULL;
            delete [] dervYZD32;
            dervYZD32 = NULL;
            delete [] dervZXD13;
            dervZXD13 = NULL;
            delete [] dervZYD23;
            dervZYD23 = NULL;
            delete [] dervXD12;
            dervXD12 = NULL;
            delete [] dervZD32;
            dervZD32 = NULL;
            delete [] dervXD13;
            dervXD13 = NULL;
            delete [] dervYD23;
            dervYD23 = NULL;
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
        }
//        double l2normErrorOri = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight);
//        printf("%lf\n", l2normErrorOri);
        l2normError = l2Norm(tempImgArrayCurr,scatImageArr,imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE error: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE error: " << aaerror << "\n";

    delete [] d1111;
    d1111 = NULL;
    delete [] d1112;
    d1112 = NULL;
    delete [] d1113;
    d1113 = NULL;
    delete [] d1121;
    d1121 = NULL;
    delete [] d1122;
    d1122 = NULL;
    delete [] d1123;
    d1123 = NULL;
    delete [] d1131;
    d1131 = NULL;
    delete [] d1132;
    d1132 = NULL;
    delete [] d1133;
    d1133 = NULL;

    delete [] d1211;
    d1211 = NULL;
    delete [] d1212;
    d1212 = NULL;
    delete [] d1213;
    d1213 = NULL;
    delete [] d1221;
    d1221 = NULL;
    delete [] d1222;
    d1222 = NULL;
    delete [] d1223;
    d1223 = NULL;
    delete [] d1231;
    d1231 = NULL;
    delete [] d1232;
    d1232 = NULL;
    delete [] d1233;
    d1233 = NULL;

    delete [] d1311;
    d1311 = NULL;
    delete [] d1312;
    d1312 = NULL;
    delete [] d1313;
    d1313 = NULL;
    delete [] d1321;
    d1321 = NULL;
    delete [] d1322;
    d1322 = NULL;
    delete [] d1323;
    d1323 = NULL;
    delete [] d1331;
    d1331 = NULL;
    delete [] d1332;
    d1332 = NULL;
    delete [] d1333;
    d1333 = NULL;


    delete [] d2111;
    d2111 = NULL;
    delete [] d2112;
    d2112 = NULL;
    delete [] d2113;
    d2113 = NULL;
    delete [] d2121;
    d2121 = NULL;
    delete [] d2122;
    d2122 = NULL;
    delete [] d2123;
    d2123 = NULL;
    delete [] d2131;
    d2131 = NULL;
    delete [] d2132;
    d2132 = NULL;
    delete [] d2133;
    d2133 = NULL;

    delete [] d2211;
    d2211 = NULL;
    delete [] d2212;
    d2212 = NULL;
    delete [] d2213;
    d2213 = NULL;
    delete [] d2221;
    d2221 = NULL;
    delete [] d2222;
    d2222 = NULL;
    delete [] d2223;
    d2223 = NULL;
    delete [] d2231;
    d2231 = NULL;
    delete [] d2232;
    d2232 = NULL;
    delete [] d2233;
    d2233 = NULL;

    delete [] d2311;
    d2311 = NULL;
    delete [] d2312;
    d2312 = NULL;
    delete [] d2313;
    d2313 = NULL;
    delete [] d2321;
    d2321 = NULL;
    delete [] d2322;
    d2322 = NULL;
    delete [] d2323;
    d2323 = NULL;
    delete [] d2331;
    d2331 = NULL;
    delete [] d2332;
    d2332 = NULL;
    delete [] d2333;
    d2333 = NULL;


    delete [] d3111;
    d3111 = NULL;
    delete [] d3112;
    d3112 = NULL;
    delete [] d3113;
    d3113 = NULL;
    delete [] d3121;
    d3121 = NULL;
    delete [] d3122;
    d3122 = NULL;
    delete [] d3123;
    d3123 = NULL;
    delete [] d3131;
    d3131 = NULL;
    delete [] d3132;
    d3132 = NULL;
    delete [] d3133;
    d3133 = NULL;

    delete [] d3211;
    d3211 = NULL;
    delete [] d3212;
    d3212 = NULL;
    delete [] d3213;
    d3213 = NULL;
    delete [] d3221;
    d3221 = NULL;
    delete [] d3222;
    d3222 = NULL;
    delete [] d3223;
    d3223 = NULL;
    delete [] d3231;
    d3231 = NULL;
    delete [] d3232;
    d3232 = NULL;
    delete [] d3233;
    d3233 = NULL;

    delete [] d3311;
    d3311 = NULL;
    delete [] d3312;
    d3312 = NULL;
    delete [] d3313;
    d3313 = NULL;
    delete [] d3321;
    d3321 = NULL;
    delete [] d3322;
    d3322 = NULL;
    delete [] d3323;
    d3323 = NULL;
    delete [] d3331;
    d3331 = NULL;
    delete [] d3332;
    d3332 = NULL;
    delete [] d3333;
    d3333 = NULL;

    delete [] t11;
    t11 = NULL;
    delete [] t12;
    t12 = NULL;
    delete [] t13;
    t13 = NULL;
    delete [] t21;
    t21 = NULL;
    delete [] t22;
    t22 = NULL;
    delete [] t23;
    t23 = NULL;
    delete [] t31;
    t31 = NULL;
    delete [] t32;
    t32 = NULL;
    delete [] t33;
    t33 = NULL;

    delete [] tempImgArrayCurr;
    tempImgArrayCurr = NULL;
    delete [] tempImgArrayPrev;
    tempImgArrayPrev = NULL;

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
//~ Recoursive Multi Thread Linear Homogenous Fourth Order Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void RecurMulThread_linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth)
{
    clock_t begin = clock();

    float stepSize = 1;
    double mserror, aaerror, l2normError;
    double* scatImageArrCurr;

    l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << l2normError << "\n";

    while(l2normError > tol) {
        scatImageArrCurr = FSIstepRecur_linfod_3d_inpainting_FSI(numSteps,timeStepSize,scatImageArr,randPxls,imgWidth,imgHeight,imgDepth,stepSize);
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            scatImageArr[i] = scatImageArrCurr[i];
        }

        delete [] scatImageArrCurr;
        scatImageArrCurr = NULL;

        l2normError = l2Norm(scatImageArr,imageArr,imgWidth*imgHeight*imgDepth);
        qDebug() << l2normError << "\n";
    }

    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "MSE: " << mserror << "\n";
    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
    qDebug() << "AAE: " << aaerror << "\n";

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
//***************************************************************************************************************************************************************************************************************************//
double* FSIstepRecur_linfod_3d_inpainting_FSI(int n, float timeStepSize, double* scatImageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float dervStepSize)
{
    std::lock_guard<std::mutex> lockGuard(myMutex);
    int randArrTraceIndex;
//    float stepSize = 1;
    double *dervXX, *dervYY, *dervZZ, *dervXXXX, *dervYYYY, *dervZZZZ, *dervXXYY, *dervZZXX, *dervYYZZ;
    double alpha;

    if(n == 0) {
        double* outArray = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            outArray[i] = scatImageArr[i];
        }

        return outArray;
    } else if(n == 1) {
        double* outArray = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

        dervXX = derivative3DXX(scatImageArr,imgWidth,imgHeight,imgDepth,dervStepSize);                           // Derivative w.r.t. xx.
        dervYY = derivative3DYY(scatImageArr,imgWidth,imgHeight,imgDepth,dervStepSize);                           // Derivative w.r.t. yy.
        dervZZ = derivative3DZZ(scatImageArr,imgWidth,imgHeight,imgDepth,dervStepSize);                           // Derivative w.r.t. zz.
        dervXXXX = derivative3DXX(dervXX,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. xxxx.
        dervYYYY = derivative3DYY(dervYY,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. yyyy.
        dervZZZZ = derivative3DZZ(dervZZ,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. zzzz.
        dervXXYY = derivative3DYY(dervXX,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. xxyy.
        dervZZXX = derivative3DXX(dervZZ,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. zzxx.
        dervYYZZ = derivative3DZZ(dervYY,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. yyzz.

        alpha = (double)(4*n+2)/(2*n+3);

        randArrTraceIndex = 0;
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[randArrTraceIndex]) {
                randArrTraceIndex++;
                continue;
            }
            outArray[i] = alpha*(scatImageArr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*scatImageArr[i] ;
        }
        delete [] dervXX;
        dervXX = NULL;
        delete [] dervYY;
        dervYY = NULL;
        delete [] dervZZ;
        dervZZ = NULL;
        delete [] dervXXYY;
        dervXXYY = NULL;
        delete dervYYZZ;
        dervYYZZ = NULL;
        delete [] dervZZXX;
        dervZZXX = NULL;
        delete [] dervXXXX;
        dervXXXX = NULL;
        delete [] dervYYYY;
        dervYYYY = NULL;
        delete [] dervZZZZ;
        dervZZZZ = NULL;

        return outArray;
    } else {
        double* outArray = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
//        double* tempImgArrayCurr = FSIstepRecur_linfod_3d_inpainting_FSI(n-1,timeStepSize,scatImageArr,randPxls,imgWidth,imgHeight,imgDepth,dervStepSize);
//        double* tempImgArrayPrev = FSIstepRecur_linfod_3d_inpainting_FSI(n-2,timeStepSize,scatImageArr,randPxls,imgWidth,imgHeight,imgDepth,dervStepSize);

        double *tempImgArrayCurr, *tempImgArrayPrev;
        std::thread t1([=,&tempImgArrayCurr](){tempImgArrayCurr = FSIstepRecur_linfod_3d_inpainting_FSI(n-1,timeStepSize,scatImageArr,randPxls,imgWidth,imgHeight,imgDepth,dervStepSize);});
        std::thread t2([=,&tempImgArrayPrev](){tempImgArrayPrev = FSIstepRecur_linfod_3d_inpainting_FSI(n-2,timeStepSize,scatImageArr,randPxls,imgWidth,imgHeight,imgDepth,dervStepSize);});
        t1.join();
        t2.join();

        dervXX = derivative3DXX(tempImgArrayCurr,imgWidth,imgHeight,imgDepth,dervStepSize);                           // Derivative w.r.t. xx.
        dervYY = derivative3DYY(tempImgArrayCurr,imgWidth,imgHeight,imgDepth,dervStepSize);                           // Derivative w.r.t. yy.
        dervZZ = derivative3DZZ(tempImgArrayCurr,imgWidth,imgHeight,imgDepth,dervStepSize);                           // Derivative w.r.t. zz.
        dervXXXX = derivative3DXX(dervXX,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. xxxx.
        dervYYYY = derivative3DYY(dervYY,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. yyyy.
        dervZZZZ = derivative3DZZ(dervZZ,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. zzzz.
        dervXXYY = derivative3DYY(dervXX,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. xxyy.
        dervZZXX = derivative3DXX(dervZZ,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. zzxx.
        dervYYZZ = derivative3DZZ(dervYY,imgWidth,imgHeight,imgDepth,dervStepSize);                               // Derivative w.r.t. yyzz.

        alpha = (double)(4*n+2)/(2*n+3);

        randArrTraceIndex = 0;
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(i == randPxls[randArrTraceIndex]) {
                randArrTraceIndex++;
                continue;
            }
            outArray[i] = alpha*(tempImgArrayCurr[i] - timeStepSize*(dervXXXX[i] + dervYYYY[i] + dervZZZZ[i] + 2*dervXXYY[i] + 2*dervZZXX[i] + 2*dervYYZZ[i])) + (1-alpha)*tempImgArrayPrev[i] ;
        }
        delete [] dervXX;
        dervXX = NULL;
        delete [] dervYY;
        dervYY = NULL;
        delete [] dervZZ;
        dervZZ = NULL;
        delete [] dervXXYY;
        dervXXYY = NULL;
        delete dervYYZZ;
        dervYYZZ = NULL;
        delete [] dervZZXX;
        dervZZXX = NULL;
        delete [] dervXXXX;
        dervXXXX = NULL;
        delete [] dervYYYY;
        dervYYYY = NULL;
        delete [] dervZZZZ;
        dervZZZZ = NULL;
        delete [] tempImgArrayPrev;
        tempImgArrayPrev = NULL;
        delete [] tempImgArrayCurr;
        tempImgArrayCurr = NULL;

        return outArray;
    }
}

//***************************************************************************************************************************************************************************************************************************//
void eed_3d_tensor(double** strucTensor, double** diffTensor, double* imageArr, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ) {
    clock_t begin = clock();

    int kernelSize = 3, randArrTraceIndex;
//    float stepSize = 1;
    double sigma = 1, gausKernel[kernelSize], mserror, aaerror, l2normError, alpha;
    double *outConvX, *outConvXY, *outConvXYZ, *dervY, *dervX, *dervZ, *dervXConv, *dervYConv, *dervZConv;
    double *dervYD21, *dervXD12, *dervXD13, *dervYD23, *dervZD31, *dervZD32, *dervForX, *dervBacX, *dervForY, *dervBacY, *sumForX, *sumBacX, *sumForY, *sumBacY, *sumForZ, *sumBacZ, *dervForZ, *dervBacZ;
    double v1[3], v2[3], v3[3], sv1[3], sv2[3], sv3[3];
//    double  *dervForXConv, *dervForYConv, *dervForZConv, *dervBacXConv, *dervBacYConv, *dervBacZConv;

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

    //Memory allocation for structural tensor entries dij.
    double* sd11 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd12 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd13 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd21 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd22 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd23 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd31 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd32 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
    double* sd33 = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

    gridSpcX = gridSpcY = gridSpcZ = 1;

    qDebug() << gridSpcX;

    gaussian1D_kernel(gausKernel,kernelSize,sigma);                             // Calculate the Gaussian kernel with sigma variance and size of kernelSize.

    //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
    //    qDebug() << "MSE error: " << mserror << "\n";
    //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
    //    printf("Error: %lf\n", aaerror);

    outConvX = convolution3DX(imageArr, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);     // Calculate the convolution on x axis.
    outConvXY = convolution3DY(outConvX, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);        // Calculate the convolution on y axis after x axis.
    outConvXYZ = convolution3DZ(outConvXY, imgWidth, imgHeight, imgDepth, gausKernel, kernelSize);      // Calculate the convolution on z axis after x and y axis.
    //First convolved derivatives
    dervXConv = derivative3DX(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcX);                         // Derivative of convolved image w.r.t. x.
    dervYConv = derivative3DY(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcY);                         // Derivative of convolved image w.r.t. y.
    dervZConv = derivative3DZ(outConvXYZ,imgWidth,imgHeight,imgDepth,gridSpcZ);                         // Derivative of convolved image w.r.t. z.

    //First derivatives
    dervX = derivative3DX(imageArr,imgWidth,imgHeight,imgDepth,gridSpcX);                           // Derivative w.r.t. x.
    dervY = derivative3DY(imageArr,imgWidth,imgHeight,imgDepth,gridSpcY);                           // Derivative w.r.t. y.
    dervZ = derivative3DZ(imageArr,imgWidth,imgHeight,imgDepth,gridSpcZ);                           // Derivative w.r.t. z.

    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        //Define eigenvectors and entries for the diffusion tensor.
        double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square), d1, d2, d3;
        double xy_norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i], xy_norm_i = sqrt(xy_norm_i_square);

    //                                    qDebug() << "v1" << v1[0] << v1[1] << v1[2];
    //                                    qDebug() << "v2" << v2[0] << v2[1] << v2[2];
    //                                    qDebug() << "v3" << v3[0] << v3[1] << v3[2];
    //                                    qDebug() << "V12" << v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
    //                                    qDebug() << "V13" << v1[0]*v3[0]+v1[1]*v3[1]+v1[2]*v3[2];
    //                                    qDebug() << "V32" << v3[0]*v2[0]+v3[1]*v2[1]+v3[2]*v2[2];
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


    //                double sm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
    //                if(sm == 1) {
    //                    qDebug() << "One";
    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    //                } else {
    //                    qDebug() << "Not One";
    //                    qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
    //                }
    //                qDebug() << v2[0]*v2[0] << v2[1]*v2[1] << v2[2]*v2[2] << "Sum: " << v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];
    //                qDebug() << v3[0]*v3[0] << v3[1]*v3[1] << v3[2]*v3[2] << "Sum: " << v3[0]*v3[0] + v3[1]*v3[1] + v3[2]*v3[2];

    //                double tempVal = std::max(dervForXConv[i]*dervBacXConv[i],0.0) + std::max(dervForYConv[i]*dervBacYConv[i],0.0) + std::max(dervForZConv[i]*dervBacZConv[i],0.0);
    //                d1 = charbonnier_diff(tempVal, 0.1);

        double contPar = 0.561268;
        d1 = charbonnier_diff(norm_i_square, contPar);
        d2 = 1;
        d3 = 1;

//                    //Take inverse of (v1,v2,v3) matrix
//                    double mat[3][3], invmat[3][3], determinant=0, g[3];
//                    mat[0][0] = v1[0];
//                    mat[1][0] = v1[1];
//                    mat[2][0] = v1[2];
//                    mat[0][1] = v2[0];
//                    mat[1][1] = v2[1];
//                    mat[2][1] = v2[2];
//                    mat[0][2] = v3[0];
//                    mat[1][2] = v3[1];
//                    mat[2][2] = v3[2];
//                    for(int i = 0; i < 3; i++) {
//                        determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
//                    }
//                    for(int i = 0; i < 3; i++){
//                        for(int j = 0; j < 3; j++)
//                          invmat[i][j] = ((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/determinant;
//                    }
//                    if(i==79468) {
//                        qDebug() << "Smoothed norm: " << norm_i_square << "Soothed XY norm: " << xy_norm_i_square << "Smoothed Z norm: " << dervZConv[i];

//                        for(int k = 0; k < 3; k++){
//                            qDebug() << "v1: " << v1[k];
//                            qDebug() << "v2: " << v2[k];
//                            qDebug() << "v3: " << v3[k];
//                        }
//                        for(int k = 0; k < 3; k++){
//                            qDebug() << "\n";
//                            for(int j = 0; j < 3; j++)
//                                qDebug() << mat[k][j];
//                        }
//                        qDebug() << "Inverse: " << i;
//                        for(int k = 0; k < 3; k++){
//                            for(int j = 0; j < 3; j++)
//                                qDebug() << invmat[k][j] << "\t";
//                                qDebug() << "\n";
//                        }
//                    }
//                    g[0] = invmat[0][0]*dervX[i] + invmat[0][1]*dervY[i] + invmat[0][2]*dervZ[i];
//                    g[1] = invmat[1][0]*dervX[i] + invmat[1][1]*dervY[i] + invmat[1][2]*dervZ[i];
//                    g[2] = invmat[2][0]*dervX[i] + invmat[2][1]*dervY[i] + invmat[2][2]*dervZ[i];
//                    if(i==79468) {
//                        qDebug() << "Coefficients: " << g[0] << g[1] << g[2];
//                    }
//                    d1 = charbonnier_diff(norm_i_square, estContPar);
//                    if(g[1] == 0) {
//                        d2 = 1;
//                    } else {
//                        d2 = 1/g[1];
//                    }
//                    if(g[2] == 0) {
//                        d3 = 1;
//                    } else {
//                        d3 = 1/g[2];
//                    }
//                    if(norm_i_square == 0 || xy_norm_i_square == 0) {
//                        d2 = 1;
//                        d3 = 1;
//                    } else {
//                        if(g[1] == 0) {
//                            d2 = 1;
//                        } else {
//                            d2 = 1/g[1];
//                        }
//                        if(g[2] == 0) {
//                            d3 = 1;
//                        } else {
//                            d3 = 1/g[2];
//                        }
// //                        d2 = 1/g[1];
// //                        d3 = 1/g[2];
//                    }
//                    if(i==79468) {
//                        qDebug() << "Tensor eigenvalues: " << d1 << d2 << d3;
//                    }

        //Diffusion tensor entries definition with eigenvalues d1, d2, d3 and orthgonal orthonormal eigenvectors v1, v2 and v3.
        d11[i] = d1*v1[0]*v1[0] + d2*v2[0]*v2[0] + d3*v3[0]*v3[0];
        d12[i] = d1*v1[0]*v1[1] + d2*v2[0]*v2[1] + d3*v3[0]*v3[1];
        d13[i] = d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2];
        d21[i] = d1*v1[1]*v1[0] + d2*v2[1]*v2[0] + d3*v3[1]*v3[0];
        d22[i] = d1*v1[1]*v1[1] + d2*v2[1]*v2[1] + d3*v3[1]*v3[1];
        d23[i] = d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2];
        d31[i] = d1*v1[0]*v1[2] + d2*v2[0]*v2[2] + d3*v3[0]*v3[2];
        d32[i] = d1*v1[1]*v1[2] + d2*v2[1]*v2[2] + d3*v3[1]*v3[2];
        d33[i] = d1*v1[2]*v1[2] + d2*v2[2]*v2[2] + d3*v3[2]*v3[2];

        if(i == 24324) {
            qDebug() << "Diffusion Tensor: ";
            qDebug() << d11[i] << d12[i] << d13[i];
            qDebug() << d21[i] << d22[i] << d23[i];
            qDebug() << d31[i] << d32[i] << d33[i];
            qDebug() << "\n";
        }
    }
    QDir dir("./img-outputs");
    if (!dir.exists())
        dir.mkpath(".");
    QFile file("./img-outputs/eed_diff_tensors");
    if(!file.open(QFile::WriteOnly | QFile::Text)) {
        qDebug() << "Couldn't open the file!";
    }
    QByteArray temp;
    for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
        double ind;
        char buf[20];

        ind = d11[i];
        ::sprintf(buf, "%lf", ind);
        temp.append(buf);
        temp.append("\n");

            ind = d12[i];
            ::sprintf(buf, "%lf", ind);
            temp.append(buf);
            temp.append("\n");

            ind = d13[i];
            ::sprintf(buf, "%lf", ind);
            temp.append(buf);
            temp.append("\n");

            ind = d22[i];
            ::sprintf(buf, "%lf", ind);
            temp.append(buf);
            temp.append("\n");

            ind = d23[i];
            ::sprintf(buf, "%lf", ind);
            temp.append(buf);
            temp.append("\n");

            ind = d33[i];
            ::sprintf(buf, "%lf", ind);
            temp.append(buf);
            temp.append("\n");
    }
    file.write(temp);
    file.flush();
    file.close();

    delete [] dervX;
    dervX = NULL;
    delete [] dervY;
    dervY = NULL;
    delete [] dervZ;
    dervZ = NULL;
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

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("Total time spent: %lf\n", time_spent);
}
