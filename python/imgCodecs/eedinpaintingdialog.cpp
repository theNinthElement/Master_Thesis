#include "eedinpaintingdialog.h"
#include "ui_eedinpaintingdialog.h"
#include "inpainting.h"
#include "supplementary_functions.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>

#include <QDebug>
#include "derivatives.h"
//#include "CannyEdgeDetector.h"
#include <opencv2/imgproc.hpp>
#include <QChar>
#include <QDoubleValidator>


//using namespace cv;

EEDInpaintingDialog::EEDInpaintingDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EEDInpaintingDialog)
{
    ui->setupUi(this);

    // The MSE tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_eedInpaintMSETol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The Tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_eedInpaintTol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The TimeStepSize lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_eedInpaintTimeStep->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The FED Inner Cycle lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_eedInpaintInnerFEDCycSize->setValidator(new QIntValidator(0, 999999999, this));
    // The FSI Inner Loop lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_EEDInpaintFSILoopSoze->setValidator(new QIntValidator(0, 999999999, this));
    // The FSI Inner Loop lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_EED4DBetaValue->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    ui->lineEdit_EEDTemporalDiffusion->setValidator(new QDoubleValidator(0, 999999999, 10, this));


    //Make the displayed images center of the label
    ui->label_eedInpaintImg->setAlignment(Qt::AlignCenter);
    ui->label_eedInpaintMaskImg->setAlignment(Qt::AlignCenter);

    //Make FSI is chosen by default and hidden MRE tolerance &  line edits
    ui->radioButton_eedInpaintFSI->setChecked(true);
    ui->label_eedInpaintMSETol->setVisible(false);
    ui->lineEdit_eedInpaintMSETol->setVisible(false);
    ui->label_eedInpaintInnerFEDCycSize->setVisible(false);
    ui->lineEdit_eedInpaintInnerFEDCycSize->setVisible(false);

    //Set default value
    ui->lineEdit_eedInpaintTimeStep->setPlaceholderText("0.1");
    ui->lineEdit_eedInpaintFileName->setPlaceholderText("inpaintedImage");
    ui->lineEdit_EEDInpaintFSILoopSoze->setPlaceholderText("20");
    ui->lineEdit_eedInpaintInnerFEDCycSize->setPlaceholderText("20");
    ui->lineEdit_eedInpaintTol->setPlaceholderText("0.00001");
    ui->lineEdit_eedInpaintMSETol->setPlaceholderText("100");
    ui->lineEdit_EED4DBetaValue->setPlaceholderText("0.1");
    ui->lineEdit_EEDTemporalDiffusion->setPlaceholderText("0.1");
}

EEDInpaintingDialog::~EEDInpaintingDialog()
{
    delete ui;
}

void EEDInpaintingDialog::on_radioButton_eedInpaintExplitScheme_clicked()
{
    //Hide the Number of Inner Loop steps and MSE tolerance line edits and labels
    ui->lineEdit_EEDInpaintFSILoopSoze->setVisible(false);
    ui->label_eedInpaintFSILoopSize->setVisible(false);
    ui->lineEdit_eedInpaintInnerFEDCycSize->setVisible(false);
    ui->label_eedInpaintInnerFEDCycSize->setVisible(false);
    ui->lineEdit_eedInpaintMSETol->setVisible(false);
    ui->label_eedInpaintMSETol->setVisible(false);
    ui->lineEdit_EED4DBetaValue->setVisible(false);
    ui->label_EEDBetaValue->setVisible(false);
    ui->lineEdit_EEDTemporalDiffusion->setVisible(false);
    ui->label_eedInpaintFSITemporalDiffusion->setVisible(false);


    ui->lineEdit_eedInpaintTol->setVisible(true);
    ui->label_eedInpaintTol->setVisible(true);
}

void EEDInpaintingDialog::on_radioButton_eedInpaintFED_clicked()
{
    ui->lineEdit_eedInpaintTol->setVisible(false);
    ui->label_eedInpaintTol->setVisible(false);
    ui->lineEdit_EEDInpaintFSILoopSoze->setVisible(false);
    ui->label_eedInpaintFSILoopSize->setVisible(false);

    ui->lineEdit_eedInpaintMSETol->setVisible(true);
    ui->label_eedInpaintMSETol->setVisible(true);
    ui->lineEdit_eedInpaintInnerFEDCycSize->setVisible(true);
    ui->label_eedInpaintInnerFEDCycSize->setVisible(true);
    ui->lineEdit_EED4DBetaValue->setVisible(false);
    ui->label_EEDBetaValue->setVisible(false);
    ui->lineEdit_EEDTemporalDiffusion->setVisible(false);
    ui->label_eedInpaintFSITemporalDiffusion->setVisible(false);
}

void EEDInpaintingDialog::on_radioButton_eedInpaintFSI_clicked()
{
    ui->lineEdit_eedInpaintInnerFEDCycSize->setVisible(false);
    ui->label_eedInpaintInnerFEDCycSize->setVisible(false);
    ui->lineEdit_eedInpaintMSETol->setVisible(false);
    ui->label_eedInpaintMSETol->setVisible(false);

    ui->lineEdit_eedInpaintTol->setVisible(true);
    ui->label_eedInpaintTol->setVisible(true);
    ui->lineEdit_EEDInpaintFSILoopSoze->setVisible(true);
    ui->label_eedInpaintFSILoopSize->setVisible(true);
    ui->lineEdit_EED4DBetaValue->setVisible(true);
    ui->label_EEDBetaValue->setVisible(true);
    ui->lineEdit_EEDTemporalDiffusion->setVisible(true);
    ui->label_eedInpaintFSITemporalDiffusion->setVisible(true);
}

void EEDInpaintingDialog::on_pushButton_eedInpaintUplImg_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input = nifti_image_read(fin, 1);
    if(!nim_input) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_input->dim[1];
    int imgHeight = nim_input->dim[2];
//    int imgDepth = nim_input->dim[3];
    int imgTimeLen = nim_input->dim[4];
//    int dataDim = ui->comboBox_randMask->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input);
    // Create image and set to the label
    imageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*imageObject);
    ui->label_eedInpaintImg->setPixmap(image.scaled(ui->label_eedInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);
}

void EEDInpaintingDialog::on_pushButton_eedInpaintUplRefImg_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input_ref = nifti_image_read(fin, 1);
    if(!nim_input_ref) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_input_ref->dim[1];
    int imgHeight = nim_input_ref->dim[2];
//    int imgDepth = nim_input->dim[3];
    int imgTimeLen = nim_input_ref->dim[4];
//    int dataDim = ui->comboBox_randMask->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_ref);
    // Create image and set to the label
    imageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*imageObject);
    ui->label_eedInpaintImg->setPixmap(image.scaled(ui->label_eedInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);
}


void EEDInpaintingDialog::on_pushButton_eedInpaintUplMask_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    bool loadEigenInfo = false;
    loadEigenInfo = ui->checkBox_load_eigenInfo->isChecked();

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input_mask = nifti_image_read(fin, 1);
    if(!nim_input_mask) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    //Writing pixels locations to randPxl array
    int strLen = imageDataPath.length();
    int maxIndx = 0;
    for(int i=0; i<strLen; i++) {
        int currIndx = imageDataPath.indexOf("/", i);
        if(currIndx > maxIndx) {
            maxIndx = currIndx;
        }
    }
    QString imageDataPath_subString = imageDataPath.mid(0,maxIndx);
    imageDataPath_subString = imageDataPath_subString + "/pixel_locations";


//    QFile file("./img-outputs/masks/pixel_locations");
    QFile file(imageDataPath_subString);
    if(!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Error", "File(voxel locations) is NOT open.");
        return;
    }
    QString arr= file.readAll();
    randPxlStrArr = arr.split('\n');
    file.close();
    //***************************************************
    int imgWidth = nim_input_mask->dim[1];
    int imgHeight = nim_input_mask->dim[2];
//    int imgDepth = nim_input_mask->dim[3];
    int imgTimeLen = nim_input_mask->dim[4];
//    int dataDim = ui->comboBox_randMask->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_mask);
    // Create image and set to the label
    maskImageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*maskImageObject);
    ui->label_eedInpaintMaskImg->setPixmap(image.scaled(ui->label_eedInpaintMaskImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);


    //To read eigen information
    if(loadEigenInfo == true) {
        QString imageDataPath_subString = imageDataPath.mid(0,maxIndx);
        imageDataPath_subString = imageDataPath_subString + "/eigVals.txt";
        QFile file(imageDataPath_subString);
        if (!file.open(QIODevice::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, "Error", "File(eigenvalues) is NOT open.");
            return;
        }
        QString arr= file.readAll();
        eigVals = arr.split('\n');
        file.close();


        imageDataPath_subString = imageDataPath.mid(0,maxIndx);
        imageDataPath_subString = imageDataPath_subString + "/eigVecs.txt";
        qDebug() << imageDataPath_subString;
        QFile file1(imageDataPath_subString);
        if (!file1.open(QIODevice::ReadOnly | QFile::Text)) {
            QMessageBox::warning(this, "Error", "File(eigenvectors) is NOT open.");
            return;
        }
        arr= file1.readAll();
        eigVecs = arr.split('\n');
        file1.close();

        qDebug() << "eigVals size: " << eigVals.size();
        qDebug() << "eigVecs size: " << eigVecs.size();

        //***************************************************
        // Read input dataset, including data
        imageDataPath_subString = imageDataPath.mid(0,maxIndx);
        imageDataPath_subString = imageDataPath_subString + "/predicted_registered.nii";
        QByteArray finTmp = imageDataPath_subString.toLocal8Bit();
        fin = finTmp.data();
        nim_dti = nifti_image_read(fin, 1);
        if(!nim_dti) {
            fprintf(stderr,"Failed to read DTI NIfTI image data from '%s'\n", fin);
            return;
        }
        //***************************************************
    }
}

void EEDInpaintingDialog::on_pushButton_eedInpaintRun_clicked()
{
    int innerLoopSize, innerCycleSize;
    float timeStep, tol, temporalDiffusion = 100;
    bool zeros2mask = false;
    bool loadEigenInfo = false;
    bool loadRefImage = false;
    double **eigenVals, ***eigenVecs;

    float beta = 0;


    loadEigenInfo = ui->checkBox_load_eigenInfo->isChecked();
    qDebug() << "Eingen Information loaded: " << loadEigenInfo;

    if(ui->checkBox_eedIntp_zeros2mask->isChecked()) {
        zeros2mask = true;
    }
    qDebug() << "Are zeros included: " << zeros2mask;

    if(ui->checkBox_EED3Dto4D->isChecked()) {
        loadRefImage = true;
    }
    qDebug() << "Ref Image to be used " << zeros2mask;


    if(ui->radioButton_eedInpaintExplitScheme->isChecked()) {
        if(ui->lineEdit_eedInpaintTol->text().isEmpty() || ui->lineEdit_eedInpaintTimeStep->text().isEmpty() || ui->label_eedInpaintImg->pixmap()==0 //.isNull()
                || ui->label_eedInpaintMaskImg->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter (time) step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_eedInpaintTimeStep->text()).toDouble();
            tol = (ui->lineEdit_eedInpaintTol->text()).toDouble();

            double *imageArr, *scatImageArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], percentage = randPxlStrArr.size()-1, imgDepth = nim_input->dim[3];

            percentage = (100*percentage)/(imgHeight*imgWidth*imgDepth);
            qDebug() << percentage;
            QString strPercentage = QString::number(percentage+1);
            qDebug() << strPercentage;

            // Decleartion and memory allocation
            int* randPxls = new int[sizeof(int) * (randPxlStrArr.size()-1)];

            for(int i=0; i<randPxlStrArr.size()-1; i++) {
                QString temp = randPxlStrArr[i];
                randPxls[i] = temp.toInt();
            }

//            double transMat[4][4] = nim_input->qto_ijk;
//            qDebug() << transMat[0];

            imageArr = nii_to_array(nim_input);
            scatImageArr = nii_to_array(nim_input_mask);

            if(ui->checkBox_EEDInpaintMonitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for Explicit Scheme has not been implemented yet :(");
                return;
            } else {
                if(imgDepth == 1) {
                    eed_inpainting(tol, timeStep, scatImageArr, imageArr, randPxls, imgWidth, imgHeight);
                } else {
                    qDebug() << "3D EED Explicit Scheme\n";
                    eed_3d_inpainting(tol, timeStep, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
                }
            }

            array_to_nii(nim_input,scatImageArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_eedInpaintImg->setPixmap(image.scaled(ui->label_eedInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    } else if(ui->radioButton_eedInpaintFED->isChecked()) {
        qDebug() << "FED Scheme\n";

        if(ui->lineEdit_eedInpaintTimeStep->text().isEmpty() || ui->lineEdit_eedInpaintMSETol->text().isEmpty() || ui->lineEdit_eedInpaintInnerFEDCycSize->text().isEmpty() || ui->label_eedInpaintImg->pixmap()==0 //.isNull()
                || ui->label_eedInpaintMaskImg->pixmap()==0 ) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_eedInpaintTimeStep->text()).toDouble();
            tol = (ui->lineEdit_eedInpaintMSETol->text()).toDouble();
            innerCycleSize = (ui->lineEdit_eedInpaintInnerFEDCycSize->text()).toDouble();

            double *imageArr, *scatImageArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

            // Decleartion and memory allocation
            int* randPxls = new int[sizeof(int) * (randPxlStrArr.size()-1)];

            for(int i=0; i<randPxlStrArr.size()-1; i++) {
                QString temp = randPxlStrArr[i];
                randPxls[i] = temp.toInt();
            }
            imageArr = nii_to_array(nim_input);
            scatImageArr = nii_to_array(nim_input_mask);

            if(ui->checkBox_EEDInpaintMonitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for FED Scheme has not been implemented yet :(");
                return;
            } else {
                eed_inpaintingFED_steps(tol, timeStep, innerCycleSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight);
            }

            array_to_nii(nim_input,scatImageArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_eedInpaintImg->setPixmap(image.scaled(ui->label_eedInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    } else if(ui->radioButton_eedInpaintFSI->isChecked()) {
        qDebug() << "FSI Scheme\n";

        if(!ui->lineEdit_EED4DBetaValue->text().isEmpty()) {
            beta = (ui->lineEdit_EED4DBetaValue->text()).toFloat();
        }

        if(!ui->lineEdit_EEDTemporalDiffusion->text().isEmpty()) {
            temporalDiffusion = (ui->lineEdit_EEDTemporalDiffusion->text()).toFloat();
        }


        if(ui->lineEdit_eedInpaintTimeStep->text().isEmpty() || ui->lineEdit_eedInpaintTol->text().isEmpty() || ui->lineEdit_EEDInpaintFSILoopSoze->text().isEmpty() || ui->label_eedInpaintImg->pixmap()==0 //.isNull()
                || ui->label_eedInpaintMaskImg->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_eedInpaintTimeStep->text()).toDouble();
            tol = (ui->lineEdit_eedInpaintTol->text()).toDouble();
            innerLoopSize = (ui->lineEdit_EEDInpaintFSILoopSoze->text()).toDouble();

            double *imageArr, *scatImageArr, *dtiImageArr, *refImageArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];
            float gridSpcX = nim_input->dx, gridSpcY = nim_input->dy, gridSpcZ = nim_input->dz, gridSpcT = nim_input->dt;

            // Decleartion and memory allocation
            int* randPxls = new int[sizeof(int) * (randPxlStrArr.size()-1)];
            for(int i=0; i<randPxlStrArr.size()-1; i++) {
                QString temp = randPxlStrArr[i];
                randPxls[i] = temp.toInt();
            }

            if(loadEigenInfo) {
                eigenVals = new double*[sizeof(double) * (eigVals.size()/3)];
                for(int i=0; i<(eigVals.size()/3); ++i) {
                    eigenVals[i] = new double[sizeof(double)*3];
                }

                eigenVecs = new double**[sizeof(double) * (eigVecs.size()/9)];
                for(int i=0; i<(eigVecs.size()/9); ++i) {
                    eigenVecs[i] = new double*[sizeof(double)*3];
                    for(int j=0; j<3; j++) {
                        eigenVecs[i][j] = new double[sizeof(double)*3];
                    }
                }
            }
            for(int i=0; i<eigVals.size()/3; i++) {
                for(int j=0; j<3; j++) {
                    QString temp = eigVals[3*i+j];
                    eigenVals[i][j] = temp.toDouble();
//                    qDebug() << i << j << eigenVals[i][j];
                }
            }

            qDebug() << eigVecs.size()/9 << eigVecs.size();
            for(int i=0; i<eigVecs.size()/9; i++) {
                for(int j=0; j<3; j++) {
                    for(int k=0; k<3; k++) {
                        QString temp = eigVecs[9*i+3*j+k];
                        eigenVecs[i][j][k] = temp.toDouble();
//                        qDebug() << i << j << k << eigenVecs[i][j][k];
                    }
                }
            }


            imageArr = nii_to_array(nim_input);
            scatImageArr = nii_to_array(nim_input_mask);
            refImageArr = nii_to_array(nim_input_ref);
//            if(nim_dti) {
//                dtiImageArr = nii_to_array(nim_dti);
//            }

            if(ui->checkBox_EEDInpaintMonitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for FSI Scheme has not been implemented yet :(");
                return;
            } else {
                if(imgDepth == 1) {
                    eed_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight);
                } else {
                    qDebug() << gridSpcX << gridSpcY << gridSpcZ << gridSpcT;
                    qDebug() << imgWidth << imgHeight << imgDepth;

                    if(loadEigenInfo) {
//                        qDebug() << "Spatial DTI EED: \n";
//                        spatial_dti_eed_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask, dtiImageArr);

                        qDebug() << "DTI EED: \n";
                        dti_eed_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask, eigenVals, eigenVecs);
                    } else if(loadRefImage) {
                        if(temporalDiffusion!=100)
                        {
                            qDebug() << "We are at feed_3d_with_temporal_inpainting_FSI";
                            eed_3d_with_temporal_inpainting_FSI(tol, timeStep, temporalDiffusion, innerLoopSize, scatImageArr, imageArr, refImageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask, beta);

                        } else if(beta>0)
                        {
                            qDebug() << "We are at fluctuating_DT_eed_3d_to_4d_inpainting_FSI";
                            fluctuating_DT_eed_3d_to_4d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, refImageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask, beta);

                        } else {
                            qDebug() << "We are at eed_3d_to_4d_inpainting_FSI";
//                          eed_3d_to_4d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, refImageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask);
                            eed_3d_to_4d_old_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, refImageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask);
                        }
                    }
                    else {
                        qDebug() << "3D EED_FSI\n";

//                        double **sturTen, **diffTen;
//                        eed_3d_tensor(sturTen, diffTen, imageArr, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ);

                        eed_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, zeros2mask);

                    }

//                    linfod_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
//                    qDebug() << "Linear Homogenous\n";
//                    linear_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, zeros2mask);
//                    int maskLen = randPxlStrArr.size()-1;
//                    mulThread_linfod_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, maskLen);
//                    mulThread_foeed_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
//                    RecurMulThread_linfod_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
                }
            }

//            //Slice MSE and AAE***************************************************************************
//            for(int slcNum=0; slcNum<imgDepth; slcNum++) {
//                double sum = 0, mse = 1, diff = 0, aae = 1;
//                for(int i=(slcNum)*imgWidth*imgHeight; i<(slcNum+1)*imgWidth*imgHeight; i++) {
//                    sum = sum + (imageArr[i]-scatImageArr[i])*(imageArr[i]-scatImageArr[i]);
//                    diff = diff + abs(imageArr[i]-scatImageArr[i]);
//                }
//                mse = sum/(imgWidth*imgHeight);
//                aae = diff/(imgWidth*imgHeight);

//                qDebug() << slcNum << mse << aae;
//            }


//            array_to_nii(nim_input,imageArr);
            array_to_nii(nim_input,scatImageArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_eedInpaintImg->setPixmap(image.scaled(ui->label_eedInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    }
}

void EEDInpaintingDialog::on_pushButton_eedInpaintSaveImg_clicked()
{
    bool isSaved;

    if(ui->label_eedInpaintImg->pixmap()==0 //.isNull()
            || ui->lineEdit_eedInpaintFileName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and name is not given");
    } else {
        QDir dir("./img-outputs/edge_enhancing_diffusion_based_inpainting");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/edge_enhancing_diffusion_based_inpainting/";
        int percentage = randPxlStrArr.size()-1, imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

        percentage = (100*percentage)/(imgHeight*imgWidth*imgDepth);
        qDebug() << "Known data percentage: " << percentage;
        QString strPercentage = QString::number(percentage+1);
        qDebug() << "Known data percentage: " << strPercentage;

        if(ui->radioButton_eedInpaintFSI->isChecked()) {
            nameStr = nameStr + ui->lineEdit_eedInpaintFileName->text() + "FSI" + "_innerLoopSize" + (ui->lineEdit_EEDInpaintFSILoopSoze->text()) +  "_timeStepSize" + (ui->lineEdit_eedInpaintTimeStep->text()) + "_tol" + (ui->lineEdit_eedInpaintTol->text()) + "_" + strPercentage + "percent_of_pixels" + ".png";
        } else if(ui->radioButton_eedInpaintExplitScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_eedInpaintFileName->text() + "_timeStepSize" + (ui->lineEdit_eedInpaintTimeStep->text()) + "_tol" + (ui->lineEdit_eedInpaintTol->text()) + "_" + strPercentage + "percent_of_pixels" + ".png";
        } else if(ui->radioButton_eedInpaintFED->isChecked()) {
            nameStr = nameStr + ui->lineEdit_eedInpaintFileName->text() + "FED" + "_innerCycleSize" + (ui->lineEdit_eedInpaintInnerFEDCycSize->text()) + "_timeStepSize" + (ui->lineEdit_eedInpaintTimeStep->text()) + "_tol" + (ui->lineEdit_eedInpaintMSETol->text()) + "_" + strPercentage + "percent_of_pixels" + ".png";
        }

        //***************************************************
        //Convert QString to char*
        QByteArray baTemp = nameStr.toLocal8Bit();
        const char *fout = baTemp.data();
        //***************************************************
        // Writing scattered nii file
//        array_to_nii(nim_input, scatImageArr);
        nifti_set_filenames(nim_input,fout,1,1);
        nifti_image_write(nim_input);
        //***************************************************

        isSaved = imageObject->save(nameStr, 0, -1);
        if(isSaved) {
            QMessageBox::information(this, "Saved", "The image png is saved succesfully!");
        } else {
            QMessageBox::warning(this, "Saving Problem", "The png data is NOT saved");
        }
    }
}
