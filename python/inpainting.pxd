cdef extern from "inpainting.cpp":
    void eed_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr,
                                    int* randPxls, int imgWidth, int imgHeight)
    void eed_inpaintingFED_steps(float tol, float timeStepSize, int numSteps, double* scatImageArr,
                                        double* imageArr, int* randPxls, int imgWidth, int imgHeight)
    void eed_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr,
                                        int* randPxls, int imgWidth, int imgHeight)
    void eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr,
                                            double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth,
                                            float gridSpcX, float gridSpcY, float gridSpcZ, bool zeros2mask)
    void spatial_dti_eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr,
                                                int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY,
                                                float gridSpcZ, bool zeros2mask, double* spatialDti)
    void dti_eed_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr, double* imageArr,
                                    int* randPxls, int imgWidth, int imgHeight, int imgDepth, float gridSpcX, float gridSpcY,
                                    float gridSpcZ, bool zeros2mask, double **eigVals, double ***eigVecs)
    void linear_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr,
                                       double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask)
    void eed_3d_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr,
                                int* randPxls, int imgWidth, int imgHeight, int imgDepth)
    void linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr,
                                        double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth)
    void foeed_3d_inpainting_FSI(float tol, float timeStep, int numSteps, double* scatImageArr,
                                        double* imageArr, int* randPxls, int imgWidth, int imgHeight, int imgDepth, bool zeros2mask)
    void eed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr,
                                        double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth,
                                        int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
    void foeed_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr, double** imageArr,
                                        double** inpaintedImageArr, int imgWidth, int imgHeight, int imgDepth,
                                        int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
    void linear_4d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr,
                                                double** imageArr, double** inpaintedImageArr, int imgWidth,
                                                int imgHeight, int imgDepth, int imgTimeLen, float gridSpcX,
                                                float gridSpcY, float gridSpcZ)
    void eed_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr,
                                                double** imageArr, double** inpaintedImageArr, double** refImageArr,
                                                int imgWidth, int imgHeight, int imgDepth, int imgTimeLen,
                                                float gridSpcX, float gridSpcY, float gridSpcZ)
    void foeed_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps, double** scatImageArr,
                                                double** imageArr, double** inpaintedImageArr, double** refImageArr,
                                                int imgWidth, int imgHeight, int imgDepth, int imgTimeLen,
                                                float gridSpcX, float gridSpcY, float gridSpcZ)
    void eed_4d_signalDropoutImputation_ExplicitScheme(float tol, float timeStepSize, double** scatImageArr,
                                                double** imageArr, double** inpaintedImageArr, double** refImageArr,
                                                int imgWidth, int imgHeight, int imgDepth, int imgTimeLen,
                                                float gridSpcX, float gridSpcY, float gridSpcZ)
    void eed_with_qSpace_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps,
                                                double** scatImageArr, double** imageArr, double** inpaintedImageArr,
                                                double** refImageArr, int imgWidth, int imgHeight, int imgDepth,
                                                int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
    void eed_with_qSpace_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps,
                                                double** scatImageArr, double** imageArr, double** inpaintedImageArr,
                                                double** refImageArr, int imgWidth, int imgHeight, int imgDepth,
                                                int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
    void eed_shore_inter_4d_signalDropoutImputation_FSI(float tol, float timeStepSize, int numSteps,
                                                double** scatImageArr, double** imageArr, double** inpaintedImageArr,
                                                double** refImageArr, int imgWidth, int imgHeight, int imgDepth,
                                                int imgTimeLen, float gridSpcX, float gridSpcY, float gridSpcZ)
    void average_4d_signalDropoutImputation(double** scatImageArr, double** imageArr, double** inpaintedImageArr,
                                                double** refImageArr, int imgWidth, int imgHeight, int imgDepth,
                                                int imgTimeLen)
    void mulThread_linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps, double* scatImageArr,
                                                double* imageArr, int* randPxls, int imgWidth, int imgHeight,
                                                int imgDepth, int maskLen)
    void mulThread_foeed_3d_inpainting_FSI(float tol, float timeStep, int numSteps,
                                                double* scatImageArr, double* imageArr, int* randPxls,
                                                int imgWidth, int imgHeight, int imgDepth)
    void RecurMulThread_linfod_3d_inpainting_FSI(float tol, float timeStepSize, int numSteps,
                                                double* scatImageArr, double* imageArr, int* randPxls,
                                                int imgWidth, int imgHeight, int imgDepth)
    double* FSIstepRecur_linfod_3d_inpainting_FSI(int n, float timeStepSize, double* scatImageArr,
                                                int* randPxls, int imgWidth, int imgHeight, int imgDepth,
                                                float dervStepSize)
    void eed_3d_tensor(double** strucTensor, double** diffTensor, double* imageArr, int imgWidth, int imgHeight,
                                                int imgDepth, float gridSpcX, float gridSpcY, float gridSpcZ)


