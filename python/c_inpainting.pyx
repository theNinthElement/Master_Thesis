from inpaiting cimport inpaitning


def c_eed_inpainting(float tol, float timeStepSize, double* scatImageArr, double* imageArr,
                                    int* randPxls, int imgWidth, int imgHeight)
    return inpaitning.eed_inpainting(tol, timeStepSize, &scatImageArr, &imageArr, &randPxls, imgWidth, imgHeight)

