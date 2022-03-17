#include "primizedialog.h"
#include "ui_primizedialog.h"
#include "supplementary_functions.h"
#include "masks.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QTextStream>
#include <algorithm>

#include <QDebug>

PrimizeDialog::PrimizeDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::PrimizeDialog)
{
    ui->setupUi(this);

    //Make the displayed image center of the label
    ui->label_primImg_origImg->setAlignment(Qt::AlignCenter);

    //Set default value
    ui->lineEdit_primImgName->setPlaceholderText("maskImg");

    //Hide comboBox as of now
    ui->comboBox_primImg_ImgDim->setVisible(false);
}

PrimizeDialog::~PrimizeDialog()
{
    delete ui;
}

void PrimizeDialog::on_pushButton_primImg_imgUpl_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input_origImg = nifti_image_read(fin, 1);
    if(!nim_input_origImg) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_input_origImg->dim[1];
    int imgHeight = nim_input_origImg->dim[2];
    int imgDepth = nim_input_origImg->dim[3];
    int imgTimeLen = nim_input_origImg->dim[4];
//    int dataDim = ui->comboBox_regMask_ImgDim->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_origImg);
    // Create image and set to the label
    imagePrimImgObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imagePrimImg = QPixmap::fromImage(*imagePrimImgObject);
    ui->label_primImg_origImg->setPixmap(imagePrimImg.scaled(ui->label_primImg_origImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);
}

void PrimizeDialog::on_pushButton_primImgSave_clicked()
{
    bool isSaved;

    if(ui->label_primImg_origImg->pixmap()==0 //.isNull()
            || ui->lineEdit_primImgName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
    } else {
        QDir dir("./img-outputs/primize");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/primize/";
        nameStr = nameStr + ui->lineEdit_primImgName->text() +  "_" + "primized" + ".png";
        isSaved = imagePrimImgObject->save(nameStr, 0, -1);
    }
    if(isSaved) {
        QMessageBox::information(this, "Saved", "The image is saved succesfully!");
    } else {
        QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
    }
}

void PrimizeDialog::on_pushButton_primImg_run_clicked()
{
    int imgWidth, imgHeight, imgDepth, imgTimeLen;

    if(ui->label_primImg_origImg->pixmap()==0) //.isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter percentage and/or load an image");
    } else {
//        int regGridMaskArrLen = 0;                                        // It is needed for the reg.grd. selection by sampling size not by a ratio!
        double *imageArr;

        imgWidth = nim_input_origImg->dim[1];
        imgHeight = nim_input_origImg->dim[2];
        imgDepth = nim_input_origImg->dim[3];
        imgTimeLen = nim_input_origImg->dim[4];

        //Check if image data is a 4D data
        if(imgTimeLen != 1) {
            QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
            return;
        }

        imageArr = nii_to_array(nim_input_origImg);

        double maxVoxVal = *(max_element(imageArr, imageArr+imgWidth*imgHeight*imgDepth));
        qDebug() << "Max: " << maxVoxVal;
//        primize(imageArr, imgWidth, imgHeight, imgDepth, maxVoxVal);                                                    // To take a primized image
        primeCountImage(imageArr, imgWidth, imgHeight, imgDepth);                                                         // To take a prime counts image

        // Writing image to nifti file
        array_to_nii(nim_input_origImg, imageArr);
        QString nameStr = "./img-outputs/primize/";
        nameStr = nameStr + ui->lineEdit_primImgName->text() +  "_" + "primized" + ".nii";
        //***************************************************
        //Convert QString to char*
        QByteArray baTemp = nameStr.toLocal8Bit();
        const char *fout = baTemp.data();
        //***************************************************
        nifti_set_filenames(nim_input_origImg, fout, 1, 1);
        nifti_image_write(nim_input_origImg);

        scaleToUchar(imageArr, imgWidth*imgHeight*imgDepth);
        array_to_png(imageArr, imagePrimImgObject, imgWidth, imgHeight);

        imagePrimImg = QPixmap::fromImage(*imagePrimImgObject);
        ui->label_primImg_origImg->setPixmap(imagePrimImg.scaled(ui->label_primImg_origImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//        ui->label_randMaks_ProcImg->setPixmap(imageMask);

        delete [] imageArr;
        imageArr = NULL;

        QMessageBox::information(this, "Saved", "The code has run succesfully!");
    }
}
