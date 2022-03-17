#include "foeedinpaintingdialog.h"
#include "ui_foeedinpaintingdialog.h"
#include "inpainting.h"
#include "supplementary_functions.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <cmath>
#include <QDoubleValidator>
#include <QIntValidator>

#include <QDebug>
#include "derivatives.h"

FOEEDInpaintingDialog::FOEEDInpaintingDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FOEEDInpaintingDialog)
{
    ui->setupUi(this);

    // The MSE tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_foeedInpt_mseTol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The Tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_foeedInpt_tol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The TimeStepSize lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_foeedInpt_timeStep->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The FED Inner Cycle lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_foeedInpt_FEDinnerCycle->setValidator(new QIntValidator(0, 999999999, this));
    // The FSI Inner Loop lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_foeedInpt_FSIInnerLoopSize->setValidator(new QIntValidator(0, 999999999, this));

    //Make the displayed images center of the label
    ui->label_foeedInpt_imgDisplay->setAlignment(Qt::AlignCenter);
    ui->label_foeedInpt_maskDisplay->setAlignment(Qt::AlignCenter);

    //Make FSI is chosen by default and hidden MRE tolerance &  line edits
    ui->radioButton_foeedInpt_FSISceheme->setChecked(true);
    ui->label_foeedInpt_mseTol->setVisible(false);
    ui->lineEdit_foeedInpt_mseTol->setVisible(false);
    ui->label_foeedInpt_FEDinnerCycle->setVisible(false);
    ui->lineEdit_foeedInpt_FEDinnerCycle->setVisible(false);

    //Set default value
    ui->lineEdit_foeedInpt_timeStep->setPlaceholderText("0.01");
    ui->lineEdit_foeedInpt_saveDatName->setPlaceholderText("inpaintedImage");
    ui->lineEdit_foeedInpt_FSIInnerLoopSize->setPlaceholderText("20");
    ui->lineEdit_foeedInpt_FEDinnerCycle->setPlaceholderText("20");
    ui->lineEdit_foeedInpt_tol->setPlaceholderText("0.0001");
    ui->lineEdit_foeedInpt_mseTol->setPlaceholderText("100");
}

FOEEDInpaintingDialog::~FOEEDInpaintingDialog()
{
    delete ui;
}

void FOEEDInpaintingDialog::on_radioButton_foeedInpt_ExplicitScheme_clicked()
{
    //Hide the Number of Inner Loop steps and MSE tolerance line edits and labels
    ui->lineEdit_foeedInpt_FSIInnerLoopSize->setVisible(false);
    ui->label_foeedInpt_FSIinnerLoopSize->setVisible(false);
    ui->lineEdit_foeedInpt_FEDinnerCycle->setVisible(false);
    ui->label_foeedInpt_FEDinnerCycle->setVisible(false);
    ui->lineEdit_foeedInpt_mseTol->setVisible(false);
    ui->label_foeedInpt_mseTol->setVisible(false);

    ui->lineEdit_foeedInpt_tol->setVisible(true);
    ui->label_foeedInpt_Tol->setVisible(true);
}

void FOEEDInpaintingDialog::on_radioButton_foeedInpt_FEDScheme_clicked()
{
    ui->lineEdit_foeedInpt_tol->setVisible(false);
    ui->label_foeedInpt_Tol->setVisible(false);
    ui->lineEdit_foeedInpt_FSIInnerLoopSize->setVisible(false);
    ui->label_foeedInpt_FSIinnerLoopSize->setVisible(false);

    ui->lineEdit_foeedInpt_mseTol->setVisible(true);
    ui->label_foeedInpt_mseTol->setVisible(true);
    ui->lineEdit_foeedInpt_FEDinnerCycle->setVisible(true);
    ui->label_foeedInpt_FEDinnerCycle->setVisible(true);
}

void FOEEDInpaintingDialog::on_radioButton_foeedInpt_FSISceheme_clicked()
{
    ui->lineEdit_foeedInpt_FEDinnerCycle->setVisible(false);
    ui->label_foeedInpt_FEDinnerCycle->setVisible(false);
    ui->lineEdit_foeedInpt_mseTol->setVisible(false);
    ui->label_foeedInpt_mseTol->setVisible(false);

    ui->lineEdit_foeedInpt_tol->setVisible(true);
    ui->label_foeedInpt_Tol->setVisible(true);
    ui->lineEdit_foeedInpt_FSIInnerLoopSize->setVisible(true);
    ui->label_foeedInpt_FSIinnerLoopSize->setVisible(true);
}

void FOEEDInpaintingDialog::on_pushButton_foeedInpt_origDataLoad_clicked()
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
    ui->label_foeedInpt_imgDisplay->setPixmap(image.scaled(ui->label_foeedInpt_imgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void FOEEDInpaintingDialog::on_pushButton_foeedInpt_loadMaskData_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

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
    QFile file("./img-outputs/masks/pixel_locations");
    if(!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Error", "File is NOT open.");
    }
    QString arr= file.readAll();
    randPxlStrArr = arr.split('\n');
    file.close();
    //***************************************************
    int imgWidth = nim_input_mask->dim[1];
    int imgHeight = nim_input_mask->dim[2];
 //    int imgDepth = nim_input_mask->dim[3];
    int imgTimeLen = nim_input_mask->dim[4];
 //   int dataDim = ui->comboBox_randMask->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_mask);
    // Create image and set to the label
    maskImageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*maskImageObject);
    ui->label_foeedInpt_maskDisplay->setPixmap(image.scaled(ui->label_foeedInpt_maskDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void FOEEDInpaintingDialog::on_pushButton_foeedInpt_run_clicked()
{
    int innerLoopSize, innerCycleSize;
    float timeStep, tol;
    bool zeros2mask = false;

    if(ui->checkBox_foeedIntp_zeros2mask->isChecked()) {
        zeros2mask = true;
    }

    if(ui->radioButton_foeedInpt_ExplicitScheme->isChecked()) {
        if(ui->lineEdit_foeedInpt_tol->text().isEmpty() || ui->lineEdit_foeedInpt_timeStep->text().isEmpty() || ui->label_foeedInpt_imgDisplay->pixmap()==0 // .isNull()
                || ui->label_foeedInpt_maskDisplay->pixmap() ==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter (time) step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_foeedInpt_timeStep->text()).toDouble();
            tol = (ui->lineEdit_foeedInpt_tol->text()).toDouble();

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

            imageArr = nii_to_array(nim_input);
            scatImageArr = nii_to_array(nim_input_mask);

            if(ui->checkBox_foeedInpt_monitor->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for Explicit Scheme has not been implemented yet :(");
                return;
            } else {
                if(imgDepth == 1) {
                    QMessageBox::information(this, "Monitoring", "Monitoring for Explicit Scheme has not been implemented yet :(");
                    return;
                } else {
                    QMessageBox::information(this, "Missing Functionality", "FOEED with Explicit Scheme has not been implemented yet :(");
                    return;
                }
            }

            array_to_nii(nim_input,scatImageArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_foeedInpt_imgDisplay->setPixmap(image.scaled(ui->label_foeedInpt_imgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    } else if(ui->radioButton_foeedInpt_FEDScheme->isChecked()) {
        qDebug() << "FED Scheme\n";

        if(ui->lineEdit_foeedInpt_timeStep->text().isEmpty() || ui->lineEdit_foeedInpt_mseTol->text().isEmpty() || ui->lineEdit_foeedInpt_FEDinnerCycle->text().isEmpty() || ui->label_foeedInpt_imgDisplay->pixmap() == 0 //.isNull()
                || ui->label_foeedInpt_maskDisplay->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_foeedInpt_timeStep->text()).toDouble();
            tol = (ui->lineEdit_foeedInpt_mseTol->text()).toDouble();
            innerCycleSize = (ui->lineEdit_foeedInpt_FEDinnerCycle->text()).toDouble();

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

            if(ui->checkBox_foeedInpt_monitor->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for FED Scheme has not been implemented yet :(");
                return;
            } else {
                QMessageBox::information(this, "Missing Functionality", "FOEED with FED Scheme has not been implemented yet :(");
                return;
            }

            array_to_nii(nim_input,scatImageArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_foeedInpt_imgDisplay->setPixmap(image.scaled(ui->label_foeedInpt_imgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    } else if(ui->radioButton_foeedInpt_FSISceheme->isChecked()) {
        qDebug() << "FSI Scheme\n";

        if(ui->lineEdit_foeedInpt_timeStep->text().isEmpty() || ui->lineEdit_foeedInpt_tol->text().isEmpty() || ui->lineEdit_foeedInpt_FSIInnerLoopSize->text().isEmpty() || ui->label_foeedInpt_imgDisplay->pixmap()==0 //.isNull()
                || ui->label_foeedInpt_maskDisplay->pixmap()==0)  //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_foeedInpt_timeStep->text()).toDouble();
            tol = (ui->lineEdit_foeedInpt_tol->text()).toDouble();
            innerLoopSize = (ui->lineEdit_foeedInpt_FSIInnerLoopSize->text()).toDouble();

            double *imageArr, *scatImageArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];
            float gridSpcX = nim_input->dx, gridSpcY = nim_input->dy, gridSpcZ = nim_input->dz, gridSpcT = nim_input->dt;

            // Decleartion and memory allocation
            int* randPxls = new int[sizeof(int) * (randPxlStrArr.size()-1)];

            for(int i=0; i<randPxlStrArr.size()-1; i++) {
                QString temp = randPxlStrArr[i];
                randPxls[i] = temp.toInt();
            }
            imageArr = nii_to_array(nim_input);
            scatImageArr = nii_to_array(nim_input_mask);

            if(ui->checkBox_foeedInpt_monitor->isChecked()) {
//                QMessageBox::information(this, "Monitoring", "Monitoring for FSI Scheme has not been implemented yet :(");
//                return;

//                qDebug() << "Voxel value: " << imageArr[56+imgHeight*112+imgHeight*imgWidth*38] << "Mask voxel value: " << scatImageArr[56+imgHeight*112+imgHeight*imgWidth*38];

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

//                //Initialization with average mask pixels
//                double vxlSum = 0;
//                randArrTraceIndex = 0;
//                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//                    if(i == randPxls[randArrTraceIndex]) {
//                        vxlSum = vxlSum + scatImageArr[i];
//                        randArrTraceIndex++;
//                    }
//                    continue;
//                }
//                double maskSize = randArrTraceIndex;
//                qDebug() << "Unknown voxel initialization: " << vxlSum/maskSize;
//                randArrTraceIndex = 0;
//                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//                    if(i == randPxls[randArrTraceIndex]) {
//                        randArrTraceIndex++;
//                        continue;
//                    }
//                    scatImageArr[i] = vxlSum/maskSize;
//                }
//                //End of initialization process

            //    mserror = mse(imageArr,scatImageArr,imgWidth*imgHeight);
            //    qDebug() << "MSE error: " << mserror << "\n";
            //    aaerror = aae(imageArr,scatImageArr,imgWidth*imgHeight);
            //    printf("Error: %lf\n", aaerror);
            //    qDebug() << "Here";
                l2normError = l2Norm(imageArr,scatImageArr,imgWidth*imgHeight*imgDepth);
                qDebug() << "l2 norm error: " << l2normError << "\n";

                N = innerLoopSize;
                while(l2normError > tol) {
                    for(int n=0; n<N; n++) {
                        array_to_nii(nim_input, scatImageArr);
                        unsigned char* niiImgArr = nii_to_ucharArray(nim_input);
                        // Create image and set to the label
                        imageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
                        image = QPixmap::fromImage(*imageObject);
                        ui->label_foeedInpt_imgDisplay->setPixmap(image.scaled(ui->label_foeedInpt_imgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//                            array_to_png(scatImageArr, maskImageObject, imgWidth, imgHeight);
//                            image = QPixmap::fromImage(*maskImageObject);
//                            ui->label_eedInpaintImg->setPixmap(image);
                        //***************************************************
                        //Convert QString to char*
//                        QByteArray baTemp = nameStr.toLocal8Bit();
//                        const char *fout = baTemp.data();
//                        const char *fout = "outData";
//                        //***************************************************
//                        // Writing scattered nii file
//                        array_to_nii(nim_input, scatImageArr);
//                        nifti_set_filenames(nim_input,fout,1,1);
//                        nifti_image_write(nim_input);
//                        //***************************************************
                        qApp->processEvents();

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

//                            double a1, a2, b1, b2, c1, c2;

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

//                                qDebug() << "Block1" << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum:" << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] << norm_i_square;
                            } else if(xy_norm_i_square == 0) {
                                v1[0] = 0;
                                v1[1] = 0;
                                v1[2] = dervZConv[i]/normi;

                                v2[0] = 0;
                                v2[1] = 1;
                                v2[2] = 0;

                                v3[0] = 1;
                                v3[1] = 0;
                                v3[2] = 0;

//                                qDebug() << "Block2" << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum:" << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] << norm_i_square;
                            } else {
                                v1[0] = dervXConv[i]/normi;
                                v1[1] = dervYConv[i]/normi;
                                v1[2] = dervZConv[i]/normi;

//                                double a = dervXConv[i]/normi;
//                                double b = dervYConv[i]/normi;
//                                double c = dervZConv[i]/normi;
//                                double a = v1[0];
//                                double b = v1[1];
//                                double c = v1[2];

                                v2[0] = dervYConv[i]/xy_norm_i;
                                v2[1] = -dervXConv[i]/xy_norm_i;
                                v2[2] = 0;

                                v3[0] = (dervXConv[i]*dervZConv[i])/(xy_norm_i*normi);
                                v3[1] = (dervYConv[i]*dervZConv[i])/(xy_norm_i*normi);
                                v3[2] = -xy_norm_i_square/(xy_norm_i*normi);

//                                double sv1 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
//                                qDebug() << "Block3" << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum:" << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] << dervXConv[i] << dervYConv[i] << dervZConv[i] << norm_i_square << normi << a*a+b*b+c*c;
//                                qDebug() << "Block3" << v1[0] << a << v1[1] << b << v1[2] << c << "Sum:" << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] << norm_i_square << a*a+b*b+c*c;
//                                return;
                            }

//                            //Define eigenvectors and entries for the diffusion tensor.
//                            double norm_i_square = dervXConv[i]*dervXConv[i] + dervYConv[i]*dervYConv[i] + dervZConv[i]*dervZConv[i], normi = sqrt(norm_i_square);
//                            double mutDiff_norm_i_square = (dervZConv[i]-dervYConv[i])*(dervZConv[i]-dervYConv[i]) + (dervXConv[i]-dervZConv[i])*(dervXConv[i]-dervZConv[i]) + (dervYConv[i]-dervXConv[i])*(dervYConv[i]-dervXConv[i]), mutDiff_norm_i = sqrt(mutDiff_norm_i_square);
// //                            double crossPr_norm_i_square = (dervZConv[i]-dervYConv)*(dervZConv[i]-dervYConv[i]) + (dervXConv[i]-dervZConv[i])*(dervXConv[i]-dervZConv[i]) + (dervYConv[i]-dervXConv[i])*(dervYConv[i]-dervXConv[i]), corssPr_norm_i = sqrt(crossPr_norm_i_square);
//                            double v1[2], v2[2], v3[3], m1, m2, m3, m4, m5, m6;
//                            double e1[3][3], e2[3][3], e3[3][3], e4[3][3], e5[3][3], e6[3][3];

// //                            dervXConv[i] = 1;
// //                            dervYConv[i] = 1;
// //                            dervZConv[i] = 1;

//                            if(norm_i_square == 0) {
//                                v1[0] = 1;
//                                v1[1] = 0;
//                                v1[2] = 0;

//                                v2[0] = 0;
//                                v2[1] = 1;
//                                v2[2] = 0;

//                                v3[0] = 0;
//                                v3[1] = 0;
//                                v3[2] = 1;

//                                qDebug() << "Block1" << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
//                            } else if(dervXConv[i] == dervYConv[i] && dervXConv[i] == dervZConv[i] && dervYConv[i] == dervZConv[i]) {
//                                v1[0] = 1/sqrt(3);
//                                v1[1] = 1/sqrt(3);
//                                v1[2] = 1/sqrt(3);

//                                v2[0] = 0;
//                                v2[1] = 0;
//                                v2[2] = 0;

//                                v3[0] = 0;
//                                v3[1] = 0;
//                                v3[2] = 0;

//                                qDebug() << "Block2" << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
// //                                qDebug() << dervXConv[i] << dervYConv[i] << dervZConv[i];
//                            } else {
// //                                v1[0] = dervXConv[i]/normi;
// //                                v1[1] = dervYConv[i]/normi;
// //                                v1[2] = dervZConv[i]/normi;
//                                v1[0] = dervXConv[i];
//                                v1[1] = dervYConv[i];
//                                v1[2] = dervZConv[i];

//                                v2[0] = dervZConv[i]-dervYConv[i];
//                                v2[1] = dervXConv[i]-dervZConv[i];
//                                v2[2] = dervYConv[i]-dervXConv[i];

//                                v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
//                                v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
//                                v3[2] = v1[0]*v2[1]-v1[1]*v2[0];

//                                //Normilization
//                                v1[0] = v1[0]/normi;
//                                v1[1] = v1[1]/normi;
//                                v1[2] = v1[2]/normi;

//                                v2[0] = v2[0]/mutDiff_norm_i;
//                                v2[1] = v2[1]/mutDiff_norm_i;
//                                v2[2] = v2[2]/mutDiff_norm_i;

//                                double v3Norm = sqrt(v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2]);
//                                v3[0] = v3[0]/v3Norm;
//                                v3[1] = v3[1]/v3Norm;
//                                v3[2] = v3[2]/v3Norm;

// //                                if(normi == 0)
// //                                    qDebug() << normi;
// //                                if(mutDiff_norm_i == 0)
// //                                    qDebug() << mutDiff_norm_i;
// //                                if(v3Norm == 0)
// //                                    qDebug() << v3Norm;
//                                qDebug() << "Block3" << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2] << "Vector: " << v1[0] << v1[1] << v1[2];
//                            }
//                            qDebug() << v1[0]*v1[0] << v1[1]*v1[1] << v1[2]*v1[2] << "Sum: " << v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];

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
//                            m4 = (m1+m2)/2;
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

                        delete [] niiImgArr;
                        niiImgArr = NULL;
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
            } else {
                if(imgDepth == 1) {
//                    eed_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight);
                    return;
                } else {
                    qDebug() << "3D FOEED_FSI\n";
                    qDebug() << gridSpcX << gridSpcY << gridSpcZ << gridSpcT;
                    foeed_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, zeros2mask);
//                    linfod_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
//                    int maskLen = randPxlStrArr.size()-1;
//                    mulThread_linfod_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, maskLen);
//                    mulThread_foeed_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
//                    RecurMulThread_linfod_3d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth);
                }
            }

            array_to_nii(nim_input,scatImageArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_foeedInpt_imgDisplay->setPixmap(image.scaled(ui->label_foeedInpt_imgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    }
}

void FOEEDInpaintingDialog::on_pushButton_foeedInpt_saveDatName_clicked()
{
    bool isSaved;

    if(ui->label_foeedInpt_imgDisplay->pixmap()==0 //.isNull()
            || ui->lineEdit_foeedInpt_saveDatName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and name is not given");
    } else {
        QDir dir("./img-outputs/FOEED_based_inpainting");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/FOEED_based_inpainting/";
        int percentage = randPxlStrArr.size()-1, imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

        percentage = (100*percentage)/(imgHeight*imgWidth*imgDepth);
        qDebug() << "Known data percentage: " << percentage;
        QString strPercentage = QString::number(percentage+1);
        qDebug() << "Known data percentage: " << strPercentage;

        if(ui->radioButton_foeedInpt_FSISceheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_foeedInpt_saveDatName->text() + "FSI" + "_innerLoopSize" + (ui->lineEdit_foeedInpt_FSIInnerLoopSize->text()) +  "_timeStepSize" + (ui->lineEdit_foeedInpt_timeStep->text()) + "_tol" + (ui->lineEdit_foeedInpt_tol->text()) + "_" + strPercentage + "percent_of_pixels";
        } else if(ui->radioButton_foeedInpt_ExplicitScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_foeedInpt_saveDatName->text() + "_timeStepSize" + (ui->lineEdit_foeedInpt_timeStep->text()) + "_tol" + (ui->lineEdit_foeedInpt_tol->text()) + "_" + strPercentage + "percent_of_pixels";
        } else if(ui->radioButton_foeedInpt_FEDScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_foeedInpt_saveDatName->text() + "FED" + "_innerCycleSize" + (ui->lineEdit_foeedInpt_FEDinnerCycle->text()) + "_timeStepSize" + (ui->lineEdit_foeedInpt_timeStep->text()) + "_tol" + (ui->lineEdit_foeedInpt_mseTol->text()) + "_" + strPercentage + "percent_of_pixels";
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
