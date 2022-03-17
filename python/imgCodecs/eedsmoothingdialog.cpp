#include "eedsmoothingdialog.h"
#include "ui_eedsmoothingdialog.h"
#include "supplementary_functions.h"
#include "smoothing.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QDoubleValidator>

#include <QDebug>
#include "derivatives.h"

EEDSmoothingDialog::EEDSmoothingDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::EEDSmoothingDialog)
{
    ui->setupUi(this);

    // The Iter Number lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_EEDSmoothing_iterNum->setValidator(new QIntValidator(0, 999999999, this));
    // The Number of Cycles lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_EEDSmoothing_FEDInnerCycles->setValidator(new QIntValidator(0, 999999999, this));
    // The Total Stopping Time lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_EEDSmoothing_TotTime->setValidator(new QIntValidator(0, 999999999, this));
    // The Step size lineedit will only accept doubles between 0 and 999999999 with up to decimals digits after the decimal point.
    ui->lineEdit_EEDSmoothing_timeStep->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The Contrast Parameter lineedit will only accept doubles between 0 and 999999999 with up to decimals digits after the decimal point.
    ui->lineEdit_EEDSmoothing_contPar->setValidator(new QDoubleValidator(0, 999999999, 10, this));

    //Make the displayed image center of the label
    ui->label_EEDSmoothing_ImgDisplay->setAlignment(Qt::AlignCenter);

    //Make FED is chosen by default and hidden stepSize & iteration Number line edits
    ui->radioButton_EEDSmoothing_FED->setChecked(true);
    ui->label__EEDSmoothing_timeStep->setVisible(false);
    ui->lineEdit_EEDSmoothing_timeStep->setVisible(false);
    ui->label__EEDSmoothing_iterNum->setVisible(false);
    ui->lineEdit_EEDSmoothing_iterNum->setVisible(false);

    //Set default value
    ui->lineEdit_EEDSmoothing_iterNum->setPlaceholderText("100");
    ui->lineEdit_EEDSmoothing_timeStep->setPlaceholderText("0.15");
    ui->lineEdit_EEDSmoothing_contPar->setPlaceholderText("3.5");
    ui->lineEdit__EEDSmoothing_saveImgName->setPlaceholderText("smoothedImage");
    ui->lineEdit_EEDSmoothing_FEDInnerCycles->setPlaceholderText("25");
    ui->lineEdit_EEDSmoothing_TotTime->setPlaceholderText("500");
}

EEDSmoothingDialog::~EEDSmoothingDialog()
{
    delete ui;
}

void EEDSmoothingDialog::on_radioButton_EEDSmoothing_explScheme_clicked()
{
    //Unhide the stepSize & iteration Number line edits
    ui->label__EEDSmoothing_timeStep->setVisible(true);
    ui->lineEdit_EEDSmoothing_timeStep->setVisible(true);
    ui->label__EEDSmoothing_iterNum->setVisible(true);
    ui->lineEdit_EEDSmoothing_iterNum->setVisible(true);

    //Hide the Number of Cycles and Total Stopping time line edits
    ui->label_EEDSmoothing_FEDInnerCycles->setVisible(false);
    ui->lineEdit_EEDSmoothing_FEDInnerCycles->setVisible(false);
    ui->label_EEDSmoothing_totTime->setVisible(false);
    ui->lineEdit_EEDSmoothing_TotTime->setVisible(false);
}

void EEDSmoothingDialog::on_radioButton_EEDSmoothing_FED_clicked()
{
    //Hiden the stepSize & iteration Number line edits
    ui->label__EEDSmoothing_timeStep->setVisible(false);
    ui->lineEdit_EEDSmoothing_timeStep->setVisible(false);
    ui->label__EEDSmoothing_iterNum->setVisible(false);
    ui->lineEdit_EEDSmoothing_iterNum->setVisible(false);

    //Unhide the Number of Cycles and Total Stopping time line edits
    ui->label_EEDSmoothing_FEDInnerCycles->setVisible(true);
    ui->lineEdit_EEDSmoothing_FEDInnerCycles->setVisible(true);
    ui->label_EEDSmoothing_totTime->setVisible(true);
    ui->lineEdit_EEDSmoothing_TotTime->setVisible(true);
}

void EEDSmoothingDialog::on_pushButton_EEDSmoothing_svaeImgName_clicked()
{
    bool isSaved;

    if(ui->label_EEDSmoothing_ImgDisplay->pixmap()==0  //.isNull()
            || ui->lineEdit__EEDSmoothing_saveImgName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and name is not given");
    } else {
        QDir dir("./img-outputs/edge_enhancing_diffusion_based_smoothing");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/edge_enhancing_diffusion_based_smoothing/";
//        int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

        if(ui->radioButton_EEDSmoothing_explScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit__EEDSmoothing_saveImgName->text() + "_timeStepSize" + (ui->lineEdit_EEDSmoothing_timeStep->text()) + "_stopTime" + (ui->lineEdit_EEDSmoothing_TotTime->text()) + "_contPar" + (ui->lineEdit_EEDSmoothing_contPar->text()) + ".png";
        } else if(ui->radioButton_EEDSmoothing_FED->isChecked()) {
            nameStr = nameStr + ui->lineEdit__EEDSmoothing_saveImgName->text() + "FED" + "_innerCycleSize" + (ui->lineEdit_EEDSmoothing_FEDInnerCycles->text()) + "_timeStepSize" + (ui->lineEdit_EEDSmoothing_timeStep->text()) + "_stopTime" + (ui->lineEdit_EEDSmoothing_TotTime->text()) + "_contPar" + (ui->lineEdit_EEDSmoothing_contPar->text()) + ".png";
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

void EEDSmoothingDialog::on_pushButton__EEDSmoothing_LoadImgData_clicked()
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
    ui->label_EEDSmoothing_ImgDisplay->setPixmap(image.scaled(ui->label_EEDSmoothing_ImgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void EEDSmoothingDialog::on_pushButton_EEDSmoothing_run_clicked()
{
    int innerCycleSize, iterNum;
    float timeStep, contPar;

    if(ui->radioButton_EEDSmoothing_explScheme->isChecked()) {
        if(ui->lineEdit_EEDSmoothing_iterNum->text().isEmpty() || ui->lineEdit_EEDSmoothing_timeStep->text().isEmpty() || ui->label_EEDSmoothing_ImgDisplay->pixmap()==0 //.isNull()
                || ui->lineEdit_EEDSmoothing_contPar->text().isEmpty()) {
            QMessageBox::warning(this, "Missing Input", "Please enter (time) step size and/or tol and/or image");
        } else {
            timeStep = (ui->lineEdit_EEDSmoothing_timeStep->text()).toDouble();
            contPar = (ui->lineEdit_EEDSmoothing_contPar->text()).toDouble();
            iterNum = (ui->lineEdit_EEDSmoothing_iterNum->text()).toDouble();

            double *imageArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

            imageArr = nii_to_array(nim_input);

            if(imgDepth == 1) {
//                eed_inpainting(tol, timeStep, scatImageArr, imageArr, randPxls, imgWidth, imgHeight);
                QMessageBox::warning(this, "Missing Function", "EED smothing with Explicit Scheme for 3D has not been implemented yet :(");
                return;
            } else {
                qDebug() << "3D EED Explicit Scheme\n";

                eed_3d_smoothing(timeStep, iterNum, contPar, imageArr, imgWidth, imgHeight, imgDepth);

//                QMessageBox::warning(this, "Missing Function", "EED smothing with Explicit Scheme for 2D has not been implemented yet :(");
//                return;
            }

            array_to_nii(nim_input,imageArr);

            scaleToUchar(imageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(imageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_EEDSmoothing_ImgDisplay->setPixmap(image.scaled(ui->label_EEDSmoothing_ImgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    } else if(ui->radioButton_EEDSmoothing_FED->isChecked()) {
        qDebug() << "FED Scheme\n";

        if(ui->lineEdit_EEDSmoothing_timeStep->text().isEmpty() || ui->lineEdit_EEDSmoothing_contPar->text().isEmpty() || ui->lineEdit_EEDSmoothing_FEDInnerCycles->text().isEmpty() || ui->label_EEDSmoothing_ImgDisplay->pixmap() == 0 //.isNull()
                || ui->label_EEDSmoothing_totTime->pixmap() ==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_EEDSmoothing_timeStep->text()).toDouble();
            contPar = (ui->lineEdit_EEDSmoothing_contPar->text()).toDouble();
            innerCycleSize = (ui->lineEdit_EEDSmoothing_FEDInnerCycles->text()).toDouble();

            double *imageArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];

            imageArr = nii_to_array(nim_input);

//            eed_inpaintingFED_steps(tol, timeStep, innerCycleSize, scatImageArr, imageArr, randPxls, imgWidth, imgHeight);
            QMessageBox::warning(this, "Missing Function", "EED smothing with FED Scheme has not been implemented yet :(");
            return;

            array_to_nii(nim_input,imageArr);

            scaleToUchar(imageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(imageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_EEDSmoothing_ImgDisplay->setPixmap(image.scaled(ui->label_EEDSmoothing_ImgDisplay->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageArr;
            imageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    }
}
