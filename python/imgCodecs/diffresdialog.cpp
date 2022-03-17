#include "diffresdialog.h"
#include "ui_diffresdialog.h"
#include "supplementary_functions.h"
#include "nifti1_io.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>

#include <QDebug>

DiffResDialog::DiffResDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DiffResDialog)
{
    ui->setupUi(this);

    ui->lineEdit_diffRes_ImgName->setPlaceholderText("DiffResidual");
}

DiffResDialog::~DiffResDialog()
{
    delete ui;
}

void DiffResDialog::on_pushButton_diffRes_origImgLoad_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nimOrig = nifti_image_read(fin, 1);
    if(!nimOrig) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nimOrig->dim[1];
    int imgHeight = nimOrig->dim[2];
//    int imgDepth = nimOrig->dim[3];
    int imgTimeLen = nimOrig->dim[4];

    unsigned char* niiImgArr = nii_to_ucharArray(nimOrig);
    // Create image and set to the label
    origImageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageOrig = QPixmap::fromImage(*origImageObject);
    ui->label_diffRes_OrigImg->setPixmap(imageOrig.scaled(ui->label_diffRes_OrigImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_xorRes_OrigImg->setPixmap(imageMask);
}

void DiffResDialog::on_pushButton_diffRes_recImgLoad_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nimRec = nifti_image_read(fin, 1);
    if(!nimRec) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nimRec->dim[1];
    int imgHeight = nimRec->dim[2];
//    int imgDepth = nimRec->dim[3];
    int imgTimeLen = nimRec->dim[4];

    unsigned char* niiImgArr = nii_to_ucharArray(nimRec);
    // Create image and set to the label
    recImageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageRec = QPixmap::fromImage(*recImageObject);
    ui->label_diffRes_recImg->setPixmap(imageRec.scaled(ui->label_diffRes_recImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_xorRes_recImg->setPixmap(imageMask);
}

void DiffResDialog::on_pushButton_diffRes_run_clicked()
{
    if(ui->label_diffRes_OrigImg->pixmap()==0 || ui->label_diffRes_recImg->pixmap() ==0) // .isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter load all image data");
    } else {
        int imgWidth = nimOrig->dim[1], imgHeight = nimOrig->dim[2], imgDepth = nimOrig->dim[3];
        double *origImageArr, *recImageArr;

        //Memory allocation for residual array.
        short int* resArr = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];

        origImageArr = nii_to_array(nimOrig);
        recImageArr = nii_to_array(nimRec);

        //Post-processing: Clipping **********************************
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(recImageArr[i] < 0) {
                recImageArr[i] = 0;
            }
        }
        //End of post-processing *************************************

        double mse = 0;
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            resArr[i] = (short int)(origImageArr[i] - recImageArr[i]);
            mse = mse + resArr[i]*resArr[i];
        }
        mse = mse/(imgWidth*imgHeight*imgDepth);
        qDebug() << mse;

        si_array_to_nii(nimRec,resArr);


// //Combining EED and FOEED reconstructions ****************************************************************************************************************************************************
// //        QString dataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));
//        QString dataPath = "/home/ikram/Desktop/project/2package/build-imgCodecs-Desktop-Debug/img-outputs/randMask/20percent/foeedFSI.nii";
//        qDebug() << dataPath;
//        //***************************************************
//        //Convert QString to char*
//        QByteArray baTemp = dataPath.toLocal8Bit();
//        const char *fin = baTemp.data();
//        //***************************************************
//        // Read input dataset, including data
//        nifti_image* nimImg = nifti_image_read(fin, 1);
//        if(!nimImg) {
//            fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
//            return;
//        }
//        //***************************************************
//        double* recImageArr1 = nii_to_array(nimImg);

//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            double temp, temp1;
//            temp = abs(origImageArr[i] - recImageArr[i]);
//            temp1 = abs(origImageArr[i] - recImageArr1[i]);

// //            qDebug() << temp << temp1 << i << origImageArr[i] << recImageArr[i] << recImageArr1[i];

//            if(temp > temp1) {
//                resArr[i] = (short int)recImageArr1[i];
//            } else {
//                resArr[i] = (short int)recImageArr[i];
//            }
//        }

//        double mse = 0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            double diff = (origImageArr[i] - resArr[i]);
//            mse = mse + diff*diff;
//        }
//        mse = mse/(imgWidth*imgHeight*imgDepth);
//        qDebug() << mse;

//        si_array_to_nii(nimRec,resArr);
// //End of Combining EED and FOEED reconstructions ****************************************************************************************************************************************************


// //MSE based residual ****************************************************************************************************************************************************
//        //Post-processing: Clipping **********************************
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(recImageArr[i] < 0) {
//                recImageArr[i] = 0;
//            }
//        }
//        //End of post-processing *************************************
//        short int* errArr = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
//        double mse = 0, sign = 1, tempVal;

//        qDebug() << origImageArr[635483] << recImageArr[635483];

//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(origImageArr[i] < recImageArr[i]) {
//                sign = -1;
// //                qDebug() << "Here" << origImageArr[i] << recImageArr[i] << sign << (origImageArr[i] < recImageArr[i]) << i;
//            } else {
//                sign = 1;
//            }
// //            if(i == 635483) {
// //                qDebug() << origImageArr[635483] << recImageArr[635483] << sign << (origImageArr[635483] < recImageArr[635483]) << (385 < 341);
// //            }
//            tempVal = sign*(sqrt(abs(origImageArr[i] - recImageArr[i])));
// //            tempVal = round(tempVal);
//            resArr[i] = (short int)tempVal;
// //            errArr[i] = abs(origImageArr[i] - recImageArr[i]) - resArr[i]*resArr[i];
//            errArr[i] = abs(origImageArr[i] - recImageArr[i]) - resArr[i]*resArr[i] - abs(resArr[i]);

//            mse = mse + resArr[i]*resArr[i];
//        }
//        mse = mse/(imgWidth*imgHeight*imgDepth);
//        qDebug() << mse;

//        si_array_to_nii(nimRec,errArr);
// //End of MSE based residual ********************************************************************************************************************************************

// //Adding ****************************************************************************************************************************************************
// //Combining mask and residual
//for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//    resArr[i] = (origImageArr[i] + recImageArr[i]);
//}

//si_array_to_nii(nimRec,resArr);
// //End of Combining mask and residual
// //****************************************************************************************************************************************************


//        // Writing scattered nii file
//        si_array_to_nii(nimRec, resArr);
//        QString nameStr = "./img-outputs/residuals/";
//        nameStr = nameStr + ui->lineEdit_xorRes_ImgName->text() +  "_" + "DIFF_residual" + ".nii";
//        //***************************************************
//        //Convert QString to char*
//        QByteArray baTemp = nameStr.toLocal8Bit();
//        const char *fout = baTemp.data();
//        //***************************************************
//        nifti_set_filenames(nimRec,fout,1,1);
//        nifti_image_write(nimRec);
//        //***************************************************

        si_scaleToUchar(resArr, imgWidth*imgHeight*imgDepth);
        si_array_to_png(resArr, recImageObject, imgWidth, imgHeight);

        imageRec = QPixmap::fromImage(*recImageObject);
        ui->label_diffRes_recImg->setPixmap(imageRec.scaled(ui->label_diffRes_recImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//        ui->label_randMaks_ProcImg->setPixmap(imageMask);

        delete [] origImageArr;
        origImageArr = NULL;
        delete [] recImageArr;
        recImageArr = NULL;
        delete [] resArr;
        resArr = NULL;

        QMessageBox::information(this, "Saved", "The code has run succesfully!");
    }
}

void DiffResDialog::on_pushButton_diffRes_saveImg_clicked()
{
    bool isSaved;

    if(ui->label_diffRes_recImg->pixmap()==0 //.isNull()
            || ui->lineEdit_diffRes_ImgName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
    } else {
        QDir dir("./img-outputs/residuals/difference");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/residuals/difference/";
        nameStr = nameStr + ui->lineEdit_diffRes_ImgName->text() +  "_" + "DIFF_residual" + ".png";
        isSaved = recImageObject->save(nameStr, 0, -1);

        //***************************************************
        //Convert QString to char*
        QByteArray baTemp = nameStr.toLocal8Bit();
        const char *fout = baTemp.data();
        //***************************************************
        // Writing scattered nii file
//        array_to_nii(nim_input, scatImageArr);
        nifti_set_filenames(nimRec,fout,1,1);
        nifti_image_write(nimRec);
        //***************************************************
    }
    if(isSaved) {
        QMessageBox::information(this, "Saved", "The image is saved succesfully!");
    } else {
        QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
    }
}
