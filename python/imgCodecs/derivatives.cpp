#include <QDebug>

//Spatial second central derivative on x axis.
double* derivativeXX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize) {
    int len = inputArHeight*inputArWidth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=1; i<len-1; i++) {
        outAr[i] = (inputAr[i+1] - 2*inputAr[i] + inputAr[i-1])/(stepSize*stepSize);
    }
    for(int k=0; k<inputArHeight; k++) {
//        outAr[k*inputArWidth] = (2*inputAr[k*inputArWidth+1] - 2*inputAr[k*inputArWidth])/(stepSize*stepSize);
//        outAr[(k+1)*inputArWidth-1] = (2*inputAr[(k+1)*inputArWidth-2] - 2*inputAr[(k+1)*inputArWidth-1])/(stepSize*stepSize);
        outAr[k*inputArWidth] = (inputAr[k*inputArWidth+1] - inputAr[k*inputArWidth])/(stepSize*stepSize);
        outAr[(k+1)*inputArWidth-1] = (inputAr[(k+1)*inputArWidth-2] - inputAr[(k+1)*inputArWidth-1])/(stepSize*stepSize);
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128.128)
}

//Spatial second central derivative on y axis.
double* derivativeYY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize) {
    int len = inputArHeight*inputArWidth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=inputArWidth; i<len-inputArWidth; i++) {
        outAr[i] = (inputAr[i+inputArWidth] - 2*inputAr[i] + inputAr[i-inputArWidth])/(stepSize*stepSize);
    }
    for(int k=0; k<inputArWidth; k++) {
//        outAr[k] = (2*inputAr[k+inputArWidth] - 2*inputAr[k] )/(stepSize*stepSize);
//        outAr[(inputArHeight-1)*inputArWidth+k] = (2*inputAr[(inputArHeight-2)*inputArWidth+k] - 2*inputAr[(inputArHeight-1)*inputArWidth+k] )/(stepSize*stepSize);
        outAr[k] = (inputAr[k+inputArWidth] - inputAr[k] )/(stepSize*stepSize);
        outAr[(inputArHeight-1)*inputArWidth+k] = (inputAr[(inputArHeight-2)*inputArWidth+k] - inputAr[(inputArHeight-1)*inputArWidth+k] )/(stepSize*stepSize);
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128.128)
}

//Spatial central derivative on x axis.
double* derivativeX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize)
{
    int len = inputArHeight*inputArWidth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=1; i<len-1; i++) {
        outAr[i] = (inputAr[i+1]-inputAr[i-1])/(2*stepSize);
    }
    for(int k=0; k<inputArHeight; k++) {
//        outAr[k*inputArWidth] = 0;
//        outAr[(k+1)*inputArWidth-1] = 0;
        outAr[k*inputArWidth] = (inputAr[k*inputArWidth+1]-inputAr[k*inputArWidth])/(2*stepSize);
        outAr[(k+1)*inputArWidth-1] = (inputAr[(k+1)*inputArWidth-1]-inputAr[(k+1)*inputArWidth-2])/(2*stepSize);
//        outAr[(k+1)*inputArWidth-1] = (inputAr[(k+1)*inputArWidth-2]-inputAr[(k+1)*inputArWidth-1])/(2*stepSize);
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128,128) for derivatives
}

//Spatial central derivative on y axis.
double* derivativeY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize) {
    int len = inputArHeight*inputArWidth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=inputArWidth; i<len-inputArWidth; i++) {
        outAr[i] = (inputAr[i+inputArWidth]-inputAr[i-inputArWidth])/(2*stepSize);
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    }
    for(int k=0; k<inputArWidth; k++) {
//        outAr[k] = 0;
//        outAr[(inputArHeight-1)*inputArWidth+k] = 0;
        outAr[k] = (inputAr[k+inputArWidth]-inputAr[k])/(2*stepSize);
        outAr[(inputArHeight-1)*inputArWidth+k] = (inputAr[(inputArHeight-1)*inputArWidth+k]-inputAr[(inputArHeight-2)*inputArWidth+k])/(2*stepSize);
//        outAr[(inputArHeight-1)*inputArWidth+k] = (inputAr[(inputArHeight-2)*inputArWidth+k]-inputAr[(inputArHeight-1)*inputArWidth+k])/(2*stepSize);
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);

    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//Spatial forward derivative on x axis.
double* dervForwX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len; i++)
    {
        if((i+1) % inputArWidth == 0) {
            outAr[i] = 0;
        } else {
            outAr[i] = (inputAr[i+1]-inputAr[i])/stepSize;
        }
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);

    return outAr;
}

//Spatial backward derivative on x axis.
double* dervBackX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len; i++)
    {
        if(i % inputArWidth == 0) {
            outAr[i] = 0;
        } else {
            outAr[i] = (inputAr[i]-inputAr[i-1])/stepSize;
        }
    }
//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;
}

//Spatial forward derivative on y axis.
double* dervForwY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len-inputArWidth; i++)
        outAr[i] = (inputAr[i+inputArWidth]-inputAr[i])/stepSize;
    for(int i=inputArWidth*(inputArHeight-1); i<inputArWidth*inputArHeight; i++)
        outAr[i] = 0;

//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    return outAr;
}

//Spatial backward derivative on y axis.
double* dervBackY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize)
{
    int len = inputArHeight*inputArWidth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<inputArWidth; i++)
        outAr[i] = 0;
    for(int i=inputArWidth; i<len; i++)
        outAr[i] = (inputAr[i]-inputAr[i-inputArWidth])/stepSize;

//    for(int i=0; i<len; i++)
//        printf("outAr[%d] = %lf\n", i, outAr[i]);

    return outAr;
}

//3D Spatial central derivative on x axis.
double* derivative3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=1; i<(slen-1); i++) {
            outAr[i+dlen] = (inputAr[i+1+dlen]-inputAr[i-1+dlen])/(2*stepSize);
        }
        for(int k=0; k<inputArHeight; k++) {
    //        outAr[k*inputArWidth] = 0;
    //        outAr[(k+1)*inputArWidth-1] = 0;
            outAr[k*inputArWidth+dlen] = (inputAr[k*inputArWidth+1+dlen]-inputAr[k*inputArWidth+dlen])/(2*stepSize);
            outAr[(k+1)*inputArWidth-1+dlen] = (inputAr[(k+1)*inputArWidth-1+dlen]-inputAr[(k+1)*inputArWidth-2+dlen])/(2*stepSize);
//            outAr[(k+1)*inputArWidth-1+dlen] = (inputAr[(k+1)*inputArWidth-2+dlen]-inputAr[(k+1)*inputArWidth-1+dlen])/(2*stepSize);
        }
    }
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128,128) for derivatives
}

//4D Spatial - Temporal central derivative on x axis.
double** derivative4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double)*inputArVolume];
    for (int i = 0; i < inputArVolume; i++)
    {
        *outAr = new double[sizeof (double) * len];
    }

    for (int v = 0; v < inputArVolume; v++){

        for(int d=0; d<inputArDepth; d++) {
            int dlen = d*inputArWidth*inputArHeight;
            int slen = inputArWidth*inputArHeight;

            for(int i=1; i<(slen-1); i++) {
                outAr[v][i+dlen] = (inputAr[v][i+1+dlen]-inputAr[v][i-1+dlen])/(2*stepSize);
            }
            for(int k=0; k<inputArHeight; k++) {
        //        outAr[k*inputArWidth] = 0;
        //        outAr[(k+1)*inputArWidth-1] = 0;
                outAr[v][k*inputArWidth+dlen] = (inputAr[v][k*inputArWidth+1+dlen]-inputAr[v][k*inputArWidth+dlen])/(2*stepSize);
                outAr[v][(k+1)*inputArWidth-1+dlen] = (inputAr[v][(k+1)*inputArWidth-1+dlen]-inputAr[v][(k+1)*inputArWidth-2+dlen])/(2*stepSize);
    //            outAr[(k+1)*inputArWidth-1+dlen] = (inputAr[(k+1)*inputArWidth-2+dlen]-inputAr[(k+1)*inputArWidth-1+dlen])/(2*stepSize);
            }
        }
    }
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128,128) for derivatives
}

//3D Spatial central derivative on y axis.
double* derivative3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=inputArWidth; i<slen-inputArWidth; i++) {
            outAr[i+dlen] = (inputAr[i+inputArWidth+dlen]-inputAr[i-inputArWidth+dlen])/(2*stepSize);
    //        printf("outAr[%d] = %lf\n", i, outAr[i]);
        }
        for(int k=0; k<inputArWidth; k++) {
    //        outAr[k] = 0;
    //        outAr[(inputArHeight-1)*inputArWidth+k] = 0;
            outAr[k+dlen] = (inputAr[k+inputArWidth+dlen]-inputAr[k+dlen])/(2*stepSize);
            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = (inputAr[(inputArHeight-1)*inputArWidth+k+dlen]-inputAr[(inputArHeight-2)*inputArWidth+k+dlen])/(2*stepSize);
//            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = (inputAr[(inputArHeight-2)*inputArWidth+k+dlen]-inputAr[(inputArHeight-1)*inputArWidth+k+dlen])/(2*stepSize);
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//4D Spatial-temporal central derivative on y axis.
double** derivative4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double* [sizeof(double) * inputArVolume];
    for (int i = 0; i < inputArVolume; i++) {
        *outAr = new double[sizeof(double) * len];
    }

    for (int v = 0; v < inputArVolume; v++) {

        for(int d=0; d<inputArDepth; d++) {
            int dlen = d*inputArWidth*inputArHeight;
            int slen = inputArWidth*inputArHeight;

            for(int i=inputArWidth; i<slen-inputArWidth; i++) {
                outAr[v][i+dlen] = (inputAr[v][i+inputArWidth+dlen]-inputAr[v][i-inputArWidth+dlen])/(2*stepSize);
        //        printf("outAr[%d] = %lf\n", i, outAr[i]);
            }
            for(int k=0; k<inputArWidth; k++) {
        //        outAr[k] = 0;
        //        outAr[(inputArHeight-1)*inputArWidth+k] = 0;
                outAr[v][k+dlen] = (inputAr[v][k+inputArWidth+dlen]-inputAr[v][k+dlen])/(2*stepSize);
                outAr[v][(inputArHeight-1)*inputArWidth+k+dlen] = (inputAr[v][(inputArHeight-1)*inputArWidth+k+dlen]-inputAr[v][(inputArHeight-2)*inputArWidth+k+dlen])/(2*stepSize);
    //            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = (inputAr[(inputArHeight-2)*inputArWidth+k+dlen]-inputAr[(inputArHeight-1)*inputArWidth+k+dlen])/(2*stepSize);
            }
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//3D Spatial central derivative on z axis.
double* derivative3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    int slen = inputArWidth*inputArHeight;

    for(int i=slen; i<len-slen; i++) {
        outAr[i] = (inputAr[i+slen]-inputAr[i-slen])/(2*stepSize);
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    }

    for(int h=0; h<inputArHeight; h++) {
        int hlen = h*inputArWidth;
        for(int k=0; k<inputArWidth; k++) {
    //        outAr[k] = 0;
    //        outAr[(inputArHeight-1)*inputArWidth+k] = 0;
            outAr[k+hlen] = (inputAr[k+slen+hlen]-inputAr[k+hlen])/(2*stepSize);
            outAr[len-slen+k+hlen] = (inputAr[len-slen+k+hlen]-inputAr[len-2*slen+k+hlen])/(2*stepSize);
//            outAr[len-slen+k+hlen] = (inputAr[len-2*slen+k+hlen]-inputAr[len-slen+k+hlen])/(2*stepSize);
        }
    }

    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//4D Spatial-temporal central derivative on z axis.
double** derivative4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i < inputArVolume; i++)
    {
        *outAr = new double[sizeof(double) * len];
    }

    for ( int v = 0; v < inputArVolume; v ++)
    {
        int slen = inputArWidth*inputArHeight;

        for(int i=slen; i<len-slen; i++) {
            outAr[v][i] = (inputAr[v][i+slen]-inputAr[v][i-slen])/(2*stepSize);
        }

        for(int h=0; h<inputArHeight; h++) {
            int hlen = h*inputArWidth;
            for(int k=0; k<inputArWidth; k++) {
                outAr[v][k+hlen] = (inputAr[k+slen+hlen]-inputAr[k+hlen])/(2*stepSize);
                outAr[v][len-slen+k+hlen] = (inputAr[len-slen+k+hlen]-inputAr[len-2*slen+k+hlen])/(2*stepSize);
            }
        }
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}

//4D Spatial-temporal central derivative on z axis.
double** derivative4DT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i < inputArVolume; i++)
    {
        *outAr = new double[sizeof(double) * len];
    }
    int v = 0;
    for( int i = 0; i<len; i++){

        v = 0;

        outAr[v][i] = (inputAr[v+1][i] - inputAr[v][i])/(2*stepSize);

        for (v = 1; v < inputArVolume-1; v++)
        {

            outAr[v][i] = (inputAr[v+1][i] - inputAr[v-1][i])/(2*stepSize);
        }

        outAr[v][i] = (inputAr[v][i] - inputAr[v-1][i])/(2*stepSize);
    }
    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}


//3D Spatial forward derivative on x axis.
double* dervForw3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<(d+1)*slen; i++) {
            if((i+1) % inputArWidth == 0) {
                outAr[i] = 0;
            } else {
                outAr[i] = (inputAr[i+1]-inputAr[i])/stepSize;
            }
        }
    }
    return outAr;
}

//4D Spatial forward derivative on x axis.
double** dervForw4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;


    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v = 0; v < inputArVolume; v++)
    {

        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<(d+1)*slen; i++) {
                if((i+1) % inputArWidth == 0) {
                    outAr[v][i] = 0;
                } else {
                    outAr[v][i] = (inputAr[v][i+1]-inputAr[v][i])/stepSize;
                }
            }
        }
    }
    return outAr;
}

//3D Spatial forward derivative on y axis.
double* dervForw3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<slen-inputArWidth+d*slen; i++)
            outAr[i] = (inputAr[i+inputArWidth]-inputAr[i])/stepSize;
        for(int i=inputArWidth*(inputArHeight-1)+d*slen; i<(d+1)*slen; i++)
            outAr[i] = 0;
    }
    return outAr;
}

//4D Spatial forward derivative on y axis.
double** dervForw4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v = 0; v < inputArVolume; v++)
    {
        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<slen-inputArWidth+d*slen; i++)
                outAr[v][i] = (inputAr[v][i+inputArWidth]-inputAr[v][i])/stepSize;
            for(int i=inputArWidth*(inputArHeight-1)+d*slen; i<(d+1)*slen; i++)
                outAr[v][i] = 0;
        }
    }
    return outAr;
}


//3D Spatial forward derivative on z axis.
double* dervForw3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=0; i<len-inputArWidth*inputArHeight; i++) {
        outAr[i] = (inputAr[i+inputArWidth*inputArHeight]-inputAr[i])/stepSize;
    }

    for(int i=len-inputArWidth*inputArHeight; i<len; i++) {
        outAr[i] = 0;
    }

    return outAr;
}

//4D Spatial forward derivative on z axis.
double** dervForw4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v = 0; v < inputArVolume; v++)
    {
        for(int i=0; i<len-inputArWidth*inputArHeight; i++)
        {
            outAr[v][i] = (inputAr[v][i+inputArWidth*inputArHeight]-inputAr[v][i])/stepSize;
        }

        for(int i=len-inputArWidth*inputArHeight; i<len; i++) {
            outAr[v][i] = 0;
        }
    }
    return outAr;
}

//4D Spatial forward derivative on t axis.
double** dervForw4DT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }
    \
    int v=0;

    for(int i=0; i<len-inputArWidth*inputArHeight; i++) {
        {
            for (v=0; v<inputArVolume-1; v++)
            {
                outAr[v][i] = (inputAr[v+1][i+len]-inputAr[v][i])/stepSize;
            }

            if(v==inputArVolume)
                outAr[v][i] = 0;

        }
    }
    return outAr;
}

//3D Spatial backward derivative on x axis.
double* dervBack3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<(d+1)*slen; i++) {
            if(i % inputArWidth == 0) {
                outAr[i] = 0;
            } else {
                outAr[i] = (inputAr[i]-inputAr[i-1])/stepSize;
            }
        }
    }
    return outAr;
}

//4D Spatial backward derivative on x axis.
double** dervBack4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v = 0; v < inputArVolume; v++) {

        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<(d+1)*slen; i++) {
                if(i % inputArWidth == 0) {
                    outAr[v][i] = 0;
                } else {
                    outAr[v][i] = (inputAr[v][i]-inputAr[v][i-1])/stepSize;
                }
            }
        }
    }
    return outAr;
}

//3D Spatial backward derivative on y axis.
double* dervBack3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int slen = inputArWidth*inputArHeight;

        for(int i=d*slen; i<d*slen+inputArWidth; i++)
            outAr[i] = 0;
        for(int i=inputArWidth+d*slen; i<(d+1)*slen; i++)
            outAr[i] = (inputAr[i]-inputAr[i-inputArWidth])/stepSize;
    }
    return outAr;
}

//4D Spatial backward derivative on y axis.
double** dervBack4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v = 0; v < inputArVolume; v++) {
        for(int d=0; d<inputArDepth; d++) {
            int slen = inputArWidth*inputArHeight;

            for(int i=d*slen; i<d*slen+inputArWidth; i++)
                outAr[v][i] = 0;
            for(int i=inputArWidth+d*slen; i<(d+1)*slen; i++)
                outAr[v][i] = (inputAr[v][i]-inputAr[v][i-inputArWidth])/stepSize;
        }
    }
    return outAr;
}

//3D Spatial backward derivative on z axis.
double* dervBack3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int i=inputArWidth*inputArHeight; i<len; i++) {
        outAr[i] = (inputAr[i]-inputAr[i-inputArWidth*inputArHeight])/stepSize;
    }

    for(int i=0; i<inputArWidth*inputArHeight; i++) {
        outAr[i] = 0;
    }

    return outAr;
}

//4D Spatial backward derivative on z axis.
double** dervBack4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth,int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v = 0; v < inputArVolume; v++) {
        for(int i=inputArWidth*inputArHeight; i<len; i++) {
            outAr[v][i] = (inputAr[v][i]-inputAr[v][i-inputArWidth*inputArHeight])/stepSize;
        }

        for(int i=0; i<inputArWidth*inputArHeight; i++) {
            outAr[v][i] = 0;
        }
    }

    return outAr;
}

//4D Spatial backward derivative on t axis.
double** dervBack4DT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth,int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;

    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    int v = 0;
    for(int i=0; i<len; i++) {

        outAr[0][i] = 0;

        for (v = 1; v < inputArVolume; v++) {
            outAr[v][i] = (inputAr[v][i]-inputAr[v-1][i])/stepSize;
        }
    }

    return outAr;
}

//3D Spatial second central derivative on x axis.
double* derivative3DXX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=1; i<(slen-1); i++) {
            outAr[i+dlen] = (inputAr[i+1+dlen] - 2*inputAr[i+dlen] + inputAr[i-1+dlen])/(stepSize*stepSize);
        }
        for(int k=0; k<inputArHeight; k++) {
    //        outAr[k*inputArWidth] = 0;
    //        outAr[(k+1)*inputArWidth-1] = 0;
            outAr[k*inputArWidth+dlen] = (inputAr[k*inputArWidth+1+dlen]-inputAr[k*inputArWidth+dlen])/(stepSize*stepSize);
//            outAr[(k+1)*inputArWidth-1+dlen] = (inputAr[(k+1)*inputArWidth-1+dlen]-inputAr[(k+1)*inputArWidth-2+dlen])/(stepSize*stepSize);
            outAr[(k+1)*inputArWidth-1+dlen] = (inputAr[(k+1)*inputArWidth-2+dlen]-inputAr[(k+1)*inputArWidth-1+dlen])/(stepSize*stepSize);
        }
    }

    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128.128)
}

//4D Spatial second central derivative on x axis.
double** derivative4DXX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for ( int v=0; v< inputArVolume; v++) {

        for(int d=0; d<inputArDepth; d++) {
            int dlen = d*inputArWidth*inputArHeight;
            int slen = inputArWidth*inputArHeight;

            for(int i=1; i<(slen-1); i++) {
                outAr[v][i+dlen] = (inputAr[v][i+1+dlen] - 2*inputAr[v][i+dlen] + inputAr[v][i-1+dlen])/(stepSize*stepSize);
            }
            for(int k=0; k<inputArHeight; k++) {
        //        outAr[k*inputArWidth] = 0;
        //        outAr[(k+1)*inputArWidth-1] = 0;
                outAr[v][k*inputArWidth+dlen] = (inputAr[v][k*inputArWidth+1+dlen]-inputAr[v][k*inputArWidth+dlen])/(stepSize*stepSize);
    //            outAr[(k+1)*inputArWidth-1+dlen] = (inputAr[(k+1)*inputArWidth-1+dlen]-inputAr[(k+1)*inputArWidth-2+dlen])/(stepSize*stepSize);
                outAr[v][(k+1)*inputArWidth-1+dlen] = (inputAr[v][(k+1)*inputArWidth-2+dlen]-inputAr[v][(k+1)*inputArWidth-1+dlen])/(stepSize*stepSize);
            }
        }
    }
    return outAr;   //when it is displayed, array should be scaled to (0,255) from (-128.128)
}

//3D Spatial second central derivative on y axis.
double* derivative3DYY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    for(int d=0; d<inputArDepth; d++) {
        int dlen = d*inputArWidth*inputArHeight;
        int slen = inputArWidth*inputArHeight;

        for(int i=inputArWidth; i<slen-inputArWidth; i++) {
//            qDebug() << "Index(-1): " << i-inputArWidth+dlen << "Value: " << inputAr[i-inputArWidth+dlen];
//            qDebug() << "Index(0): " << i+dlen << "Value: " << inputAr[i+dlen];
//            qDebug() << "Index(1): " << i+inputArWidth+dlen << "Value: " << inputAr[i+inputArWidth+dlen];

            outAr[i+dlen] = (inputAr[i+inputArWidth+dlen] - 2*inputAr[i+dlen] + inputAr[i-inputArWidth+dlen])/(stepSize*stepSize);
//            qDebug() << "Index(): " << i << "Derv. Value: " << outAr[i];
        }
        for(int k=0; k<inputArWidth; k++) {
//        outAr[k] = (2*inputAr[k+inputArWidth] - 2*inputAr[k] )/(stepSize*stepSize);
//        outAr[(inputArHeight-1)*inputArWidth+k] = (2*inputAr[(inputArHeight-2)*inputArWidth+k] - 2*inputAr[(inputArHeight-1)*inputArWidth+k] )/(stepSize*stepSize);
            outAr[k+dlen] = (inputAr[k+inputArWidth+dlen] - inputAr[k+dlen] )/(stepSize*stepSize);
            outAr[(inputArHeight-1)*inputArWidth+k+dlen] = (inputAr[(inputArHeight-2)*inputArWidth+k+dlen] - inputAr[(inputArHeight-1)*inputArWidth+k+dlen])/(stepSize*stepSize);
        }
    }

    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128.128)
}

//4D Spatial second central derivative on y axis.
double** derivative4DYY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for ( int v=0; v< inputArVolume; v++) {
        for(int d=0; d<inputArDepth; d++) {
            int dlen = d*inputArWidth*inputArHeight;
            int slen = inputArWidth*inputArHeight;

            for(int i=inputArWidth; i<slen-inputArWidth; i++) {

                outAr[v][i+dlen] = (inputAr[v][i+inputArWidth+dlen] - 2*inputAr[v][i+dlen] + inputAr[v][i-inputArWidth+dlen])/(stepSize*stepSize);
            }
            for(int k=0; k<inputArWidth; k++) {
                outAr[v][k+dlen] = (inputAr[v][k+inputArWidth+dlen] - inputAr[v][k+dlen] )/(stepSize*stepSize);
                outAr[v][(inputArHeight-1)*inputArWidth+k+dlen] = (inputAr[v][(inputArHeight-2)*inputArWidth+k+dlen] -
                        inputAr[v][(inputArHeight-1)*inputArWidth+k+dlen])/(stepSize*stepSize);
            }
        }
    }

    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128.128)
}

//3D Spatial central derivative on z axis.
double* derivative3DZZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double* outAr = new double[sizeof(double) * len];

    int slen = inputArWidth*inputArHeight;

    for(int i=slen; i<len-slen; i++) {
        outAr[i] = (inputAr[i-slen]-2*inputAr[i]+inputAr[i+slen])/(stepSize*stepSize);
//        printf("outAr[%d] = %lf\n", i, outAr[i]);
    }

    for(int h=0; h<inputArHeight; h++) {
        int hlen = h*inputArWidth;
        for(int k=0; k<inputArWidth; k++) {
    //        outAr[k] = 0;
    //        outAr[(inputArHeight-1)*inputArWidth+k] = 0;
            outAr[k+hlen] = (inputAr[k+slen+hlen]-inputAr[k+hlen])/(stepSize*stepSize);
            outAr[len-slen+k+hlen] = (inputAr[len-2*slen+k+hlen]-inputAr[len-slen+k+hlen])/(stepSize*stepSize);
        }
    }

    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}


//4D Spatial central derivative on z axis.
double** derivative4DZZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    for (int v=0; v < inputArVolume; v++) {

        int slen = inputArWidth*inputArHeight;

        for(int i=slen; i<len-slen; i++) {
            outAr[v][i] = (inputAr[v][i-slen]-2*inputAr[v][i]+inputAr[v][i+slen])/(stepSize*stepSize);
    //        printf("outAr[%d] = %lf\n", i, outAr[i]);
        }

        for(int h=0; h<inputArHeight; h++) {
            int hlen = h*inputArWidth;
            for(int k=0; k<inputArWidth; k++) {
        //        outAr[k] = 0;
        //        outAr[(inputArHeight-1)*inputArWidth+k] = 0;
                outAr[v][k+hlen] = (inputAr[v][k+slen+hlen]-inputAr[v][k+hlen])/(stepSize*stepSize);
                outAr[v][len-slen+k+hlen] = (inputAr[v][len-2*slen+k+hlen]-inputAr[v][len-slen+k+hlen])/(stepSize*stepSize);
            }
        }
    }

    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}


//4D Spatial central derivative on t axis.
double** derivative4DTT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize)
{
    int len = inputArHeight*inputArWidth*inputArDepth;
    // Decleartion and memory allocation
    double** outAr = new double*[sizeof(double) * inputArVolume];
    for (int i = 0; i<inputArVolume; i++)
    {
        outAr[i] = new double[sizeof(double)* len];
    }

    int v = 0;

    for (int i =0; i < len; i++)
    {
        v = 0;

        outAr[v][i] = (inputAr[v+1][i] - inputAr[v][i])/(stepSize*stepSize);

        for (v = 1; v < inputArVolume-1; v++) {
            outAr[v][i] = (inputAr[v-1][i] - 2*inputAr[v][i] + inputAr[v+1][i]) / (stepSize*stepSize);
        }

//        for ( int v = 1; v < inputArVolume; v++) {
            outAr[v][i] = (inputAr[v][i] - inputAr[v-1][i])/(stepSize*stepSize);
//        }

    }


    return outAr; //when it is displayed, array should be scaled to (0,255) from (-128,128) for dervatives.
}
