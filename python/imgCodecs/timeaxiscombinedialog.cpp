#include "timeaxiscombinedialog.h"
#include "ui_timeaxiscombinedialog.h"
#include "supplementary_functions.h"
#include "nifti1_io.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>

#include <QDebug>

TimeAxisCombineDialog::TimeAxisCombineDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TimeAxisCombineDialog)
{
    ui->setupUi(this);

    ui->lineEdit_timeComb_saveImgName->setPlaceholderText("Time-axisCombined");
}

TimeAxisCombineDialog::~TimeAxisCombineDialog()
{
    delete ui;
}

void TimeAxisCombineDialog::on_pushButton_timeComb_UplData1_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nimData1 = nifti_image_read(fin, 1);
    if(!nimData1) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nimData1->dim[1];
    int imgHeight = nimData1->dim[2];
    int imgDepth = nimData1->dim[3];
    int imgTimeLen = nimData1->dim[4];

    qDebug() << imgWidth << imgHeight << imgDepth << imgTimeLen;

    unsigned char* niiImgArr = nii_to_ucharArray(nimData1);
    // Create image and set to the label
    imgDataObject1 = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageData1 = QPixmap::fromImage(*imgDataObject1);
    ui->label_timeComb_imgData1->setPixmap(imageData1.scaled(ui->label_timeComb_imgData1->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void TimeAxisCombineDialog::on_pushButton_timeComb_UplData2_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nimData2 = nifti_image_read(fin, 1);
    if(!nimData2) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nimData2->dim[1];
    int imgHeight = nimData2->dim[2];
    int imgDepth = nimData2->dim[3];
    int imgTimeLen = nimData2->dim[4];

    qDebug() << imgWidth << imgHeight << imgDepth << imgTimeLen;

    unsigned char* niiImgArr = nii_to_ucharArray(nimData2);
    // Create image and set to the label
    imgDataObject2 = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageData2 = QPixmap::fromImage(*imgDataObject2);
    ui->label_timeComb_imgData2->setPixmap(imageData2.scaled(ui->label_timeComb_imgData2->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void TimeAxisCombineDialog::on_pushButton_timeComb_saveImg_clicked()
{
    bool isSaved;

    if(ui->label_timeComb_imgData1->pixmap() ==0 //.isNull()
            || ui->label_timeComb_imgData2->pixmap() ==0 //.isNull()
            || ui->lineEdit_timeComb_saveImgName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
    } else {
        QDir dir("./img-outputs/time_axis_combine");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/time_axis_combine/";
        nameStr = nameStr + ui->lineEdit_timeComb_saveImgName->text() +  "_" + "timeAxisComb" + ".png";
        isSaved = imgDataObject1->save(nameStr, 0, -1);

        //***************************************************
        //Convert QString to char*
        QByteArray baTemp = nameStr.toLocal8Bit();
        const char *fout = baTemp.data();
        //***************************************************
        // Writing scattered nii file
//        array_to_nii(nim_input, scatImageArr);
        nifti_set_filenames(nimOutData,fout,1,1);
        nifti_image_write(nimOutData);
        //***************************************************
    }
    if(isSaved) {
        QMessageBox::information(this, "Saved", "The image is saved succesfully!");
    } else {
        QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
    }
}

void TimeAxisCombineDialog::on_pushButton_timeComb_combine_clicked()
{
    if(ui->label_timeComb_imgData1->pixmap()==0 //.isNull()
            || ui->label_timeComb_imgData2->pixmap()==0) //.isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter load all image data");
    } else {
        int imgWidth1 = nimData1->dim[1], imgHeight1 = nimData1->dim[2], imgDepth1 = nimData1->dim[3], imgTimeLen1 = nimData1->dim[4];
        int imgWidth2 = nimData2->dim[1], imgHeight2 = nimData2->dim[2], imgDepth2 = nimData2->dim[3], imgTimeLen2 = nimData2->dim[4];

        if(imgWidth1 != imgWidth2 || imgHeight1 != imgHeight2 || imgDepth1 != imgDepth2) {
            QMessageBox::warning(this, "Invalid File Sizes", "Spatial Dimentions should be equal!");
            return;
        }

//        qDebug() << imgTimeLen1 << imgTimeLen2;

        double *imgArr1, *imgArr2;

//        //Memory allocation for residual array.
        double* combArr = new double[sizeof(double) * imgWidth1*imgHeight1*imgDepth1*(imgTimeLen1+imgTimeLen2)];
//        short int* resArr1 = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];

//        qDebug() << imgWidth1*imgHeight1*imgDepth1*imgTimeLen1 << imgWidth1*imgHeight1*imgDepth1*(imgTimeLen1+imgTimeLen2);

        imgArr1 = nii_to_array(nimData1);
        imgArr2 = nii_to_array(nimData2);

//        nimData2->dim[4] = 2;
//        nifti_update_dims_from_array(nimData2);                                                              // changing according sizes nt etc.
//        qDebug() << nimData2->dim[4];

        qDebug() << "Hello1";

        for(int i=0; i<imgWidth1*imgHeight1*imgDepth1*imgTimeLen1; i++) {
            combArr[i] = imgArr1[i];
        }
        qDebug() << "Hello2";
        for(int i=imgWidth1*imgHeight1*imgDepth1*imgTimeLen1; i<imgWidth1*imgHeight1*imgDepth1*(imgTimeLen1+imgTimeLen2); i++) {
            combArr[i] = imgArr2[i-imgWidth1*imgHeight1*imgDepth1*imgTimeLen1];
        }
        qDebug() << "Hello3";

//        nimData1->dim[4] = imgTimeLen1+imgTimeLen2;
        nimOutData->dim[4] = 2;
        nimOutData->dim[0] = 4;
        nifti_update_dims_from_array(nimOutData);   // changing according sizes nt etc.
        qDebug() << nimOutData->dim[4];
        qDebug() << nimOutData->dim[0];
        qDebug() << "Hello4";
        array_to_nii(nimOutData,combArr);
        qDebug() << "Hello5";

//        si_array_to_nii(nimData1,combArr);
//        si_scaleToUchar(imgArr1, imgWidth1*imgHeight1*imgDepth1);
//        si_array_to_png(resArr, recImageObject, imgWidth, imgHeight);

//        imageRec = QPixmap::fromImage(*recImageObject);
//        ui->label_xorRes_recImg->setPixmap(imageRec.scaled(ui->label_xorRes_recImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//        ui->label_randMaks_ProcImg->setPixmap(imageMask);

        delete [] imgArr1;
        imgArr1 = NULL;
        delete [] imgArr2;
        imgArr2 = NULL;
        delete [] combArr;
        combArr = NULL;

        QMessageBox::information(this, "Finished", "The code has run succesfully!");
    }
}

void TimeAxisCombineDialog::on_pushButton_outData_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nimOutData = nifti_image_read(fin, 1);
    if(!nimOutData) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nimOutData->dim[1];
    int imgHeight = nimOutData->dim[2];
    int imgDepth = nimOutData->dim[3];
    int imgTimeLen = nimOutData->dim[4];

    qDebug() << imgWidth << imgHeight << imgDepth << imgTimeLen;
}
