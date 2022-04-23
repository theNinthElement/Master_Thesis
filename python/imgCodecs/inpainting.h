#ifndef INPAINTING_H
#define INPAINTING_H

//**********************************************************************************************************************************************************************
//~ Edge Enhancing Diffusion-based Inpaiting by Explicit Scheme
void eed_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight);

//~ Edge Enhancing Diffusion-based Inpaiting by Fast Explicit Scheme (Steps)
void eed_inpaintingFED_steps(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight);

//~ Edge Enhancing Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void eed_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight);
//**********************************************************************************************************************************************************************
void eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask);
//**********************************************************************************************************************************************************************
void spatial_dti_eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask, double* spatialDti);
//**********************************************************************************************************************************************************************
void dti_eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask, double **eigVals, double ***eigVecs);
//***************************************************************************************************************************************************************************************************************************//
void linear_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask);
//***************************************************************************************************************************************************************************************************************************//
void eed_3d_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth);
//***************************************************************************************************************************************************************************************************************************//
//~ Linear Homogenous Fourth Order Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth);
//***************************************************************************************************************************************************************************************************************************//
//~ Fourth Order Edge Enhancing Anisotropic Diffusion based Inpaiting by Fast Semi-Iterative Scheme with fourth order diffusion tensor
void foeed_3d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask);


//***************************************************************************************************************************************************************************************************************************//
void eed_3d_to_4d_old_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, double* refImageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask);

//***************************************************************************************************************************************************************************************************************************//
void eed_3d_to_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, double* refImageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask);

//***************************************************************************************************************************************************************************************************************************//
void fluctuating_DT_eed_3d_to_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, double* refImageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask, float beta);

//***************************************************************************************************************************************************************************************************************************//
void eed_3d_with_temporal_inpainting_FSI(float tol, float timeStepSize, float temporalDiffusion, int numSteps, double* scatImageArr, double* imageArr, double* refImageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask, float beta);

//***************************************************************************************************************************************************************************************************************************//
void foeed_3d_with_temporal_data_inpainting_FSI(float tol, float timeStep, float temporalDiffusion, int numSteps, double* scatImageArr, double* imageArr, double* refImageArr, double* refImageArr2, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask);

void st_eed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, int** randPxls, int imgWidth, int imgHeight, int imgDepth, int imgVolume, float gridSpcX, float gridSpcY, float gridSpcZ, float gridSpct, bool zeros2mask);

void st_eed_4d_to_4D_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** refImageArr, int** randPxls, int imgWidth, int imgHeight, int imgDepth, int imgVolume, float gridSpcX, float gridSpcY, float gridSpcZ, float gridSpct, bool zeros2mask);

//~ Fourth Order Edge Enhancing Anisotropic Diffusion based Inpaiting by Fast Semi-Iterative Scheme with fourth order diffusion tensor
void st_foeed_4d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, int imgVolumes, bool zeros2mask);

//***************************************************************************************************************************************************************************************************************************//
//4D image inpainting with EED
void eed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//4D image inpainting with FOEED
void foeed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//4D image inpainting with Linear Homogenous Diffusion
void linear_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//***************************************************************************************************************************************************************************************************************************//
//4D image signal dropout imputation with EED
void eed_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr,int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//4D image signal dropout imputation with FOEED
void foeed_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);


//4D image signal dropout imputation with EED Explicit Scheme
void eed_4d_signalDropoutImputation_ExplicitScheme(float tol, float timeStepSize, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//4D image signal dropout imputation with EED
void eed_with_qSpace_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//4D image signal dropout imputation with EED
void eed_shore_inter_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ);

//4D image signal dropout imputation with Simple Averaging Neighboring Slices
void average_4d_signalDropoutImputation(double** scatImageArr, double** imageArr, double** inpaintedImageArr, double** refImageArr, int imgWidth, int imgHeight, int imgDepth, int imgTimeLen);
//***************************************************************************************************************************************************************************************************************************//
//~ Multi Thread Linear Homogenous Fourth Order Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void mulThread_linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, int maskLen);
//***************************************************************************************************************************************************************************************************************************//
//~ Multi Thread Fourth Order Edge Enhancing Anisotropic Diffusion based Inpaiting by Fast Semi-Iterative Scheme with fourth order diffusion tensor
void mulThread_foeed_3d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth);

//~ Multi Thread Linear Homogenous Fourth Order Diffusion-based Inpaiting by Fast Semi-Iterative Scheme
void RecurMulThread_linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth);
//***************************************************************************************************************************************************************************************************************************//
double* FSIstepRecur_linfod_3d_inpainting_FSI(int n, float timeStepSize, double* scatImageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, float dervStepSize);

//***************************************************************************************************************************************************************************************************************************//
void eed_3d_tensor(double** strucTensor, double** diffTensor, double* imageArr, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ);


#endif // INPAINTING_H
