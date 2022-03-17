#ifndef DERIVATIVES_H
#define DERIVATIVES_H

//Spatial second central derivative on x axis.
double* derivativeXX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial second central derivative on y axis.
double* derivativeYY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial central derivative on x axis.
double* derivativeX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial central derivative on y axis.
double* derivativeY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial forward derivative on x axis.
double* dervForwX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial backward derivative on x axis.
double* dervBackX(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial forward derivative on y axis.
double* dervForwY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//Spatial backward derivative on y axis.
double* dervBackY(double* inputAr, int inputArWidth, int inputArHeight, int stepSize);

//3D Spatial central derivative on x axis.
double* derivative3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//4D Spatial-temporal central derivative on x axis.
double** derivative4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//3D Spatial central derivative on y axis.
double* derivative3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//4D Spatial central derivative on x axis.
double** derivative4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//3D Spatial central derivative on z axis.
double* derivative3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//4D Spatial central derivative on t axis.
double** derivative4DT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial central derivative on z axis.
double** derivative4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//3D Spatial forward derivative on x axis.
double* dervForw3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial forward derivative on y axis.
double* dervForw3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial forward derivative on z axis.
double* dervForw3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial backward derivative on x axis.
double* dervBack3DX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial backward derivative on y axis.
double* dervBack3DY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial backward derivative on z axis.
double* dervBack3DZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial second central derivative on x axis.
double* derivative3DXX(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial second central derivative on y axis.
double* derivative3DYY(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//3D Spatial central derivative on z axis.
double* derivative3DZZ(double* inputAr, int inputArWidth, int inputArHeight, int inputArDepth, float stepSize);

//4D Spatial forward derivative on x axis.
double** dervForw4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial forward derivative on y axis.
double** dervForw4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial forward derivative on z axis.
double** dervForw4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial forward derivative on t axis.
double** dervForw4DT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial backward derivative on x axis.
double** dervBack4DX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial backward derivative on y axis.
double** dervBack4DY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial backward derivative on z axis.
double** dervBack4DZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial backward derivative on t axis.
double** dervBack4DT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial second central derivative on x axis.
double** derivative4DXX(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial second central derivative on y axis.
double** derivative4DYY(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial central derivative on z axis.
double** derivative4DZZ(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);

//4D Spatial central derivative on t axis.
double** derivative4DTT(double** inputAr, int inputArWidth, int inputArHeight, int inputArDepth, int inputArVolume, float stepSize);


#endif // DERIVATIVES_H
