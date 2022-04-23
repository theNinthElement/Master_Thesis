#include "eed_4_4d_dialog.h"
#include "ui_eed_4_4d_dialog.h"
#include "inpainting.h"
#include "supplementary_functions.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QDoubleValidator>
#include <QIntValidator>


#include <QDebug>

EED_4_4D_Dialog::EED_4_4D_Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EED_4_4D_Dialog)
{
    ui->setupUi(this);

    // The MSE tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_eed44dInpaintMSETol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The Tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_eed44dInpaintTol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The TimeStepSize lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_eed44dInpaintTimeStep->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The FED Inner Cycle lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_eed44dInpaintInnerFEDCycSize->setValidator(new QIntValidator(0, 999999999, this));
    // The FSI Inner Loop lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_eed44dInpaintFSILoopSoze->setValidator(new QIntValidator(0, 999999999, this));

    //Make the displayed images center of the label
    ui->label_eed44dInpaintImg->setAlignment(Qt::AlignCenter);
    ui->label_eed44dInpaintMaskImg->setAlignment(Qt::AlignCenter);

    //Make FSI is chosen by default and hidden MRE tolerance &  line edits
    ui->radioButton_eed44dInpaintFSI->setChecked(true);
    ui->label_eed44dInpaintMSETol->setVisible(false);
    ui->lineEdit_eed44dInpaintMSETol->setVisible(false);
    ui->label_eed44dInpaintInnerFEDCycSize->setVisible(false);
    ui->lineEdit_eed44dInpaintInnerFEDCycSize->setVisible(false);

    ui->pushButton_brainBinMaskLoad->setVisible(false);

    //Set default value
    ui->lineEdit_eed44dInpaintTimeStep->setPlaceholderText("0.1");
    ui->lineEdit_eed44dInpaintFileName->setPlaceholderText("inpaintedImage");
    ui->lineEdit_eed44dInpaintFSILoopSoze->setPlaceholderText("20");
    ui->lineEdit_eed44dInpaintInnerFEDCycSize->setPlaceholderText("20");
    ui->lineEdit_eed44dInpaintTol->setPlaceholderText("0.00001");
    ui->lineEdit_eed44dInpaintMSETol->setPlaceholderText("100");
}

EED_4_4D_Dialog::~EED_4_4D_Dialog()
{
    delete ui;
}

void EED_4_4D_Dialog::on_radioButton_eed44dInpaintExplitScheme_clicked()
{
    //Hide the Number of Inner Loop steps and MSE tolerance line edits and labels
    ui->lineEdit_eed44dInpaintFSILoopSoze->setVisible(false);
    ui->label_eed44dInpaintFSILoopSize->setVisible(false);
    ui->lineEdit_eed44dInpaintInnerFEDCycSize->setVisible(false);
    ui->label_eed44dInpaintInnerFEDCycSize->setVisible(false);
    ui->lineEdit_eed44dInpaintMSETol->setVisible(false);
    ui->label_eed44dInpaintMSETol->setVisible(false);

    ui->lineEdit_eed44dInpaintTol->setVisible(true);
    ui->label_eed44dInpaintTol->setVisible(true);
}

void EED_4_4D_Dialog::on_radioButton_eed44dInpaintFED_clicked()
{
    ui->lineEdit_eed44dInpaintTol->setVisible(false);
    ui->label_eed44dInpaintTol->setVisible(false);
    ui->lineEdit_eed44dInpaintFSILoopSoze->setVisible(false);
    ui->label_eed44dInpaintFSILoopSize->setVisible(false);

    ui->lineEdit_eed44dInpaintMSETol->setVisible(true);
    ui->label_eed44dInpaintMSETol->setVisible(true);
    ui->lineEdit_eed44dInpaintInnerFEDCycSize->setVisible(true);
    ui->label_eed44dInpaintInnerFEDCycSize->setVisible(true);
}

void EED_4_4D_Dialog::on_radioButton_eed44dInpaintFSI_clicked()
{
    ui->lineEdit_eed44dInpaintInnerFEDCycSize->setVisible(false);
    ui->label_eed44dInpaintInnerFEDCycSize->setVisible(false);
    ui->lineEdit_eed44dInpaintMSETol->setVisible(false);
    ui->label_eed44dInpaintMSETol->setVisible(false);

    ui->lineEdit_eed44dInpaintTol->setVisible(true);
    ui->label_eed44dInpaintTol->setVisible(true);
    ui->lineEdit_eed44dInpaintFSILoopSoze->setVisible(true);
    ui->label_eed44dInpaintFSILoopSize->setVisible(true);
}

void EED_4_4D_Dialog::on_pushButton_eed44dInpaintUplImg_clicked()
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

    qDebug() << "Fourth Dimension Length: " << imgTimeLen;

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input);
    // Create image and set to the label
    imageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    image = QPixmap::fromImage(*imageObject);
    ui->label_eed44dInpaintImg->setPixmap(image.scaled(ui->label_eed44dInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void EED_4_4D_Dialog::on_pushButton_eed44dInpaintUplMask_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));
//    qDebug() << "passed 1";
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
//    //***************************************************
    //Writing pixels locations to randPxl array
    QFile file("./img-outputs/masks/pixel_locations");
    if(!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, "Error", "File is NOT open.");
    }
//    qDebug() << "passed 2";

//    //***************************************************
    int imgWidth = nim_input_mask->dim[1];
    int imgHeight = nim_input_mask->dim[2];
    int imgDepth = nim_input_mask->dim[3];
    int imgTimeLen = nim_input_mask->dim[4];

    randPxlStrArr4D = new QStringList[sizeof(QStringList)*imgTimeLen];

    if(imgTimeLen > 1)
    {
        QString arr= file.readAll();
        QString s;
        for(int t = 0; t< imgTimeLen; t++) {
//            randPxlStrArr4D[t] = new QStringList[sizeof(QStringList)];
            s = arr.split("####")[t];
//            qDebug() << s.length();
            randPxlStrArr4D[t] = s.split('\n');
//            qDebug() << "passed 3";
        }
        file.close();
    } else
    {
        QString arr= file.readAll();
        randPxlStrArr = arr.split('\n');
        file.close();
    }

//    qDebug() << "passed 3";

    qDebug() << "Fourth Dimension Length: " << imgTimeLen;

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_mask);
    // Create image and set to the label
    maskImageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);

    image = QPixmap::fromImage(*maskImageObject);
    ui->label_eed44dInpaintMaskImg->setPixmap(image.scaled(ui->label_eed44dInpaintMaskImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void EED_4_4D_Dialog::on_pushButton_eed44dInpaintRun_clicked()
{
    int innerLoopSize;
    float timeStep, tol;

    if(ui->radioButton_eed44dInpaintExplitScheme->isChecked()) {
        if(ui->lineEdit_eed44dInpaintTol->text().isEmpty() || ui->lineEdit_eed44dInpaintTimeStep->text().isEmpty() || ui->label_eed44dInpaintImg->pixmap()==0 //.isNull()
                || ui->label_eed44dInpaintMaskImg->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter (time) step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
//            QMessageBox::warning(this, "Missing Functionality", "Explicit Scheme is not added yet :(");
            qDebug() << "Explicit Scheme\n";

            timeStep = (ui->lineEdit_eed44dInpaintTimeStep->text()).toDouble();
            tol = (ui->lineEdit_eed44dInpaintTol->text()).toDouble();

            double **imageArr, **scatImageArr, **refImageArr, **binInpaintingMaskArr, *binArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3], imgTimeLen = nim_input->dim[4];
            float gridSpcX = nim_input->dx, gridSpcY = nim_input->dy, gridSpcZ = nim_input->dz, gridSpcT = nim_input->dt;

            // Memory allocation for inpaintedImageArr.
            double** inpaintedImageArr = new double*[sizeof(double)*imgTimeLen];
            for(int i = 0; i < imgTimeLen; ++i) {
                inpaintedImageArr[i] = new double[sizeof(double)*imgWidth*imgHeight*imgDepth];
            }

            imageArr = nii4d_to_array(nim_input);
//            binInpaintingMaskArr = nii4d_to_array_mask(nim_input_mask);
//            binArr = nii_to_array(nim_input_mask);
            scatImageArr = nii4d_to_array(nim_input_mask);
            refImageArr = nii4d_to_array(nim_input_refer);

            binInpaintingMaskArr = nii4d_to_array(nim_input_mask);

//            for(int i = 0; i < imgTimeLen; ++i) {
//                for(int j=0; j<imgWidth*imgHeight*imgDepth; j++) {
//                    int tempIndx = binArr[j];
//                    if(tempIndx == 1) {
//                        inpaintedImageArr[i][j] = imageArr[i][j];
//                    } else {
//                        inpaintedImageArr[i][j] = 0;
//                    }
//                }
//            }
//            qDebug() << "So far so good";


            if(ui->checkBox_eed44dInpaintMonitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for Explicit Scheme has not been implemented yet :(");
                return;
            } else {
                for(int i = 0; i < imgTimeLen; ++i) {
                    for(int j=0; j<imgWidth*imgHeight*imgDepth; j++) {
                        int tempIndx = (int)binInpaintingMaskArr[i][j];
                        if(tempIndx == 1) {
                            inpaintedImageArr[i][j] = imageArr[i][j];
                        } else {
                            inpaintedImageArr[i][j] = imageArr[i][j];
                        }
                    }
                }
                qDebug() << "4D EED Explicit Scheme\n";
                qDebug() << gridSpcX << gridSpcY << gridSpcZ << gridSpcT;
                qDebug() << imgWidth << imgHeight << imgDepth << imgTimeLen;
                eed_4d_signalDropoutImputation_ExplicitScheme(tol, timeStep, binInpaintingMaskArr, imageArr, inpaintedImageArr, refImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
            }

            array_to_nii4d(nim_input,inpaintedImageArr);
        }
    } else if(ui->radioButton_eed44dInpaintFED->isChecked()) {
        if(ui->lineEdit_eed44dInpaintTimeStep->text().isEmpty() || ui->lineEdit_eed44dInpaintMSETol->text().isEmpty() || ui->lineEdit_eed44dInpaintInnerFEDCycSize->text().isEmpty() || ui->label_eed44dInpaintImg->pixmap()==0 //.isNull()
                || ui->label_eed44dInpaintMaskImg->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            QMessageBox::warning(this, "Missing Functionality", "Explicit Scheme is not added yet :(");
        }
    } else if(ui->radioButton_eed44dInpaintFSI->isChecked()) {
        qDebug() << "FSI Scheme\n";

        if(ui->lineEdit_eed44dInpaintTimeStep->text().isEmpty() || ui->lineEdit_eed44dInpaintTol->text().isEmpty() || ui->lineEdit_eed44dInpaintFSILoopSoze->text().isEmpty() || ui->label_eed44dInpaintImg->pixmap()==0 //.isNull()
                || ui->label_eed44dInpaintMaskImg->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_eed44dInpaintTimeStep->text()).toDouble();
            tol = (ui->lineEdit_eed44dInpaintTol->text()).toDouble();
            innerLoopSize = (ui->lineEdit_eed44dInpaintFSILoopSoze->text()).toDouble();

            double **imageArr, **scatImageArr, **refImageArr, **binInpaintingMaskArr;
//            double *binBrainMaskArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3], imgTimeLen = nim_input->dim[4];
            float gridSpcX = nim_input->dx, gridSpcY = nim_input->dy, gridSpcZ = nim_input->dz, gridSpcT = nim_input->dt;

            // Memory allocation for inpaintedImageArr.
            double** inpaintedImageArr = new double*[sizeof(double)*imgTimeLen];
            for(int i = 0; i < imgTimeLen; ++i) {
                inpaintedImageArr[i] = new double[sizeof(double)*imgWidth*imgHeight*imgDepth];
            }


            // Decleartion and memory allocation
            int** randPxls = new int*[sizeof(int) * imgTimeLen];
            for (int i = 0; i < imgTimeLen; i++) {
                randPxls [i] = new int[sizeof(int) * (randPxlStrArr4D[i].size()-1)];
            }
            for (int t = 0; t < imgTimeLen; t++) {
                for(int i=0; i<randPxlStrArr4D[t].size()-1; i++) {
                    QString temp = randPxlStrArr4D[t][i];
                    randPxls[t][i] = temp.toInt();
               }
            }

            imageArr = nii4d_to_array(nim_input);
//            scatImageArr = nii4d_to_array_mask(nim_input_mask);
            scatImageArr = nii4d_to_array(nim_input_mask);
            refImageArr = nii4d_to_array(nim_input_refer);

            binInpaintingMaskArr = nii4d_to_array(nim_input_mask);
//            binBrainMaskArr = nii_to_array(nim_input_brain_mask);

            if(ui->checkBox_eed44dInpaintMonitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for FSI Scheme has not been implemented yet :(");
                return;
            } else {
//                    qDebug() << "Slope and Inter: " << nim_input->scl_slope << nim_input->scl_inter << nim_input->byteorder << nim_input->datatype;
//                    qDebug() << "S form code and Q form code: " << nim_input->sform_code << nim_input->qform_code;
//                    //*********************************************************************************************************************************************************************

                for(int i = 0; i < imgTimeLen; ++i) {
                    for(int j=0; j<imgWidth*imgHeight*imgDepth; j++) {
                 //     qDebug() << scatImageArr[i][j];
                        int tempIndx = (int)scatImageArr[i][j];
 //                        int tempIndxBin = (int)binBrainMaskArr[j];
 //                        if(tempIndx == 1 && tempIndxBin == 0) {
                        if(tempIndx == 1) {
                            inpaintedImageArr[i][j] = imageArr[i][j];
 //                            scatImageArr[i][j] = 0;

                 //         qDebug() << "Entered: ";
                 //         inpaintedImageArr[i][j] = 0;
                 //         continue;
                        } else {
                 //         qDebug() << "Else: " << tempIndx;
                            inpaintedImageArr[i][j] = imageArr[i][j];
                        }
                    }
                }
                qDebug() << "4D EED FSI\n";
                qDebug() << gridSpcX << gridSpcY << gridSpcZ << gridSpcT;
                qDebug() << imgWidth << imgHeight << imgDepth << imgTimeLen;
//                eed_with_qSpace_4d_signalDropoutImputation_FSI(tol, timeStep, innerLoopSize, binInpaintingMaskArr, imageArr, inpaintedImageArr, refImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
//                eed_4d_signalDropoutImputation_FSI(tol, timeStep, innerLoopSize, binInpaintingMaskArr, imageArr, inpaintedImageArr, refImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
//                average_4d_signalDropoutImputation(binInpaintingMaskArr, imageArr, inpaintedImageArr, refImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen);
//                foeed_4d_signalDropoutImputation_FSI(tol, timeStep, innerLoopSize, binInpaintingMaskArr, imageArr, inpaintedImageArr, refImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
//                eed_4d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, inpaintedImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
//                foeed_4d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, inpaintedImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
//                linear_4d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, inpaintedImageArr, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ);
//                st_eed_4d_inpainting_FSI(tol, timeStep, innerLoopSize, binInpaintingMaskArr, imageArr, randPxls, imgWidth, imgHeight, imgDepth, imgTimeLen, gridSpcX, gridSpcY, gridSpcZ, gridSpcT, true);
//                eed_3d_to_4d_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageArr, imageArr, refImageArr, randPxls, imgWidth, imgHeight, imgDepth, gridSpcX, gridSpcY, gridSpcZ, true);
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

            array_to_nii4d(nim_input,inpaintedImageArr);
//              array_to_nii4d(nim_input,binInpaintingMaskArr);

//            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth*imgTimeLen);
//            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
//            image = QPixmap::fromImage(*imageObject);
//            ui->label_eed44dInpaintImg->setPixmap(image.scaled(ui->label_eed44dInpaintImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    }
}

void EED_4_4D_Dialog::on_pushButton_eed44dInpaintSaveImg_clicked()
{
    bool isSaved;

    if(ui->label_eed44dInpaintImg->pixmap()==0 //.isNull()
            || ui->lineEdit_eed44dInpaintFileName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and name is not given");
    } else {
        QDir dir("./img-outputs/edge_enhancing_diffusion_based_inpainting");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/edge_enhancing_diffusion_based_inpainting/";
        int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

//        percentage = (100*percentage)/(imgHeight*imgWidth*imgDepth);
//        qDebug() << "Known data percentage: " << percentage;
//        QString strPercentage = QString::number(percentage+1);
//        qDebug() << "Known data percentage: " << strPercentage;

        if(ui->radioButton_eed44dInpaintFSI->isChecked()) {
            nameStr = nameStr + ui->lineEdit_eed44dInpaintFileName->text() + "FSI" + "_innerLoopSize" + (ui->lineEdit_eed44dInpaintFSILoopSoze->text()) +  "_timeStepSize" + (ui->lineEdit_eed44dInpaintTimeStep->text()) + "_tol" + (ui->lineEdit_eed44dInpaintTol->text()) + ".png";
        } else if(ui->radioButton_eed44dInpaintExplitScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_eed44dInpaintFileName->text() + "_timeStepSize" + (ui->lineEdit_eed44dInpaintTimeStep->text()) + "_tol" + (ui->lineEdit_eed44dInpaintTol->text()) + ".png";
        } else if(ui->radioButton_eed44dInpaintFED->isChecked()) {
            nameStr = nameStr + ui->lineEdit_eed44dInpaintFileName->text() + "FED" + "_innerCycleSize" + (ui->lineEdit_eed44dInpaintInnerFEDCycSize->text()) + "_timeStepSize" + (ui->lineEdit_eed44dInpaintTimeStep->text()) + "_tol" + (ui->lineEdit_eed44dInpaintMSETol->text()) + ".png";
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

void EED_4_4D_Dialog::on_pushButton_eed44dReferImgData_clicked()
{
    QString referDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = referDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input_refer = nifti_image_read(fin, 1);
    if(!nim_input_refer) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_input_refer->dim[1];
    int imgHeight = nim_input_refer->dim[2];
//    int imgDepth = nim_input_refer->dim[3];
    int imgTimeLen = nim_input_refer->dim[4];

    qDebug() << "Fourth Dimension Length: " << imgTimeLen;
    unsigned char* niiReferArr = nii_to_ucharArray(nim_input_refer);
}

void EED_4_4D_Dialog::on_radioButton_eed_brainBinMask_clicked()
{
    ui->pushButton_brainBinMaskLoad->setVisible(true);
}

void EED_4_4D_Dialog::on_radioButton_eed_brainBinMask_without_clicked()
{
    ui->pushButton_brainBinMaskLoad->setVisible(false);
}

void EED_4_4D_Dialog::on_pushButton_brainBinMaskLoad_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input_brain_mask = nifti_image_read(fin, 1);
    if(!nim_input_brain_mask) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }

//    int imgWidth = nim_input_brain_mask->dim[1];
//    int imgHeight = nim_input_brain_mask->dim[2];
//    int imgDepth = nim_input_brain_mask->dim[3];
    int imgTimeLen = nim_input_brain_mask->dim[4];

    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Binary Brain Mask should be 3D!");
        return;
    }
}
