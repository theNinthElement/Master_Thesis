#include "randommask.h"
#include "ui_randommask.h"
#include "supplementary_functions.h"
#include "masks.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QTextStream>

#include <QDebug>

RandomMask::RandomMask(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RandomMask)
{
    ui->setupUi(this);

    // The Percentage lineedit will accept double between 1 and 100
    ui->lineEdit_randMask_percentage->setValidator(new QDoubleValidator(0, 100, 10, this));

    //Make the displayed image center of the label
    ui->label_randMaks_ProcImg->setAlignment(Qt::AlignCenter);
    //Make the displayed image center of the label
    ui->label_randMask_origImg->setAlignment(Qt::AlignCenter);

    //Set default value
    ui->lineEdit_randMask_percentage->setPlaceholderText("5");
    ui->lineEdit_randMask_fileName->setPlaceholderText("maskImg");

    //Make sliders invisible by default
    ui->horizontalSlider_radnMaskProcImg->setVisible(false);
    ui->horizontalSlider_randMaskOrigImg->setVisible(false);

    //Hide comboBox as of now
    ui->comboBox_randMask->setVisible(false);
}

RandomMask::~RandomMask()
{
    delete ui;
}

void RandomMask::on_pushButton_randMask_imgUpl_clicked()
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
    int imgDepth = nim_input->dim[3];
    int imgTimeLen = nim_input->dim[4];
//    int dataDim = ui->comboBox_randMask->currentIndex();

    qDebug() << "Read successfully: " << imgWidth << imgHeight << imgDepth << imgTimeLen;

    if(imgDepth == 1) {
        ui->horizontalSlider_radnMaskProcImg->setVisible(false);
        ui->horizontalSlider_randMaskOrigImg->setVisible(false);
    }

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input);

    qDebug() << "Hiiii";

    // Create image and set to the label
    imageMaskObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Aha you have entered 4D data");
//        return;
    }
    imageMask = QPixmap::fromImage(*imageMaskObject);
    ui->label_randMask_origImg->setPixmap(imageMask.scaled(ui->label_randMask_origImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);

    if(imgDepth != 1) {
        ui->horizontalSlider_randMaskOrigImg->setVisible(true);
    }
}

void RandomMask::on_pushButton_randMask_run_clicked()
{
    int imgWidth, imgHeight, imgDepth, imgTimeLen;
    double percent;

    if(ui->lineEdit_randMask_percentage->text().isEmpty() || ui->label_randMask_origImg->pixmap()==0) //.isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter percentage and/or load an image");
    } else {

        percent = (ui->lineEdit_randMask_percentage->text()).toDouble();
        imgWidth = nim_input->dim[1];
        imgHeight = nim_input->dim[2];
        imgDepth = nim_input->dim[3];
        imgTimeLen = nim_input->dim[4];

        qDebug() << percent << imgWidth << imgHeight << imgDepth << imgTimeLen;

        //Check if image data is a 4D data
        if(imgTimeLen != 1) {
            qDebug() << "Here we are about to start";

            double* *imageArr = new double*[sizeof(double) * imgTimeLen];
            double** randPxls = new double*[sizeof(double) * imgTimeLen];
            // Decleration and Memory allocation
            double** scatImageArr = new double*[sizeof(double) * imgTimeLen];
            for (int i =0;i<imgTimeLen;i++)
            {
                scatImageArr[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
                imageArr[i] = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];
                randPxls[i] =  new double[sizeof(double) * (int)((imgWidth*imgHeight*imgDepth*percent)/100)];
            }

            qDebug() << "declared a scatImageArray";


            //Initialize scarImageArr to the black(=0) image
            for (int v =0; v <imgTimeLen; v++) {
                for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                    scatImageArr[v][i] = 0;
                }
                qDebug() << "scattered data set for Volume " << v;
                randPxls[v] = randDiffPixels(percent, imgWidth*imgHeight*imgDepth);
                qDebug() << "random pixels generated for volume  " << v;
            }

            qDebug() << "random pixels generated";

            QDir dir("./img-outputs/masks");
            if (!dir.exists())
                dir.mkpath(".");
            QFile file("./img-outputs/masks/pixel_locations");
            if(!file.open(QFile::WriteOnly | QFile::Text)) {
                QMessageBox::warning(this, "Error", "File is NOT open");
            }

            imageArr = nii4d_to_array(nim_input);
            qDebug() << "Nii2Array successfully";

            QByteArray temp;
            for (int v = 0; v < imgTimeLen; v++) {
                int randArrSize = (int)(percent*imgHeight*imgWidth*imgDepth)/100;

                for(int i=0; i<randArrSize; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
                    int ind = randPxls[v][i];
                    scatImageArr[v][ind] = imageArr[v][ind];

                    char buf[9];
                    ::sprintf(buf, "%d", ind);
                    temp.append(buf);
                    temp.append("\n");
                }
                temp.append("####");
                temp.append("\n");
            }
            file.write(temp);
            file.flush();
            file.close();

            qDebug() << "Pixels written for 4d data";

            // Writing scattered nii file
            array_to_nii4d(nim_input, scatImageArr);
            QString nameStr = "./img-outputs/masks/";
            nameStr = nameStr + ui->lineEdit_randMask_fileName->text() +  "_" + (ui->lineEdit_randMask_percentage->text()) + "randMask" + ".nii";
            //***************************************************
            //Convert QString to char*
            QByteArray baTemp = nameStr.toLocal8Bit();
            const char *fout = baTemp.data();
            //***************************************************
            nifti_set_filenames(nim_input,fout,1,1);
            nifti_image_write(nim_input);
            //***************************************************

//            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
//            array_to_png(scatImageArr, imageMaskObject, imgWidth, imgHeight);

            imageMask = QPixmap::fromImage(*imageMaskObject);
            ui->label_randMaks_ProcImg->setPixmap(imageMask.scaled(ui->label_randMaks_ProcImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
    //        ui->label_randMaks_ProcImg->setPixmap(imageMask);

            if(imgDepth != 1) {
                ui->horizontalSlider_radnMaskProcImg->setVisible(true);
            }

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;
            delete [] randPxls;
            randPxls = NULL;

            QMessageBox::information(this, "Saved", "The code has run succesfully!");


        }else {
            double* imageArr, *randPxls;

            // Decleration and Memory allocation
            double* scatImageArr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

            //Initialize scarImageArr to the black(=0) image
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                scatImageArr[i] = 0;
            }
            randPxls = randDiffPixels(percent, imgWidth*imgHeight*imgDepth);

            QDir dir("./img-outputs/masks");
            if (!dir.exists())
                dir.mkpath(".");
            QFile file("./img-outputs/masks/pixel_locations");
            if(!file.open(QFile::WriteOnly | QFile::Text)) {
                QMessageBox::warning(this, "Error", "File is NOT open");
            }

            imageArr = nii_to_array(nim_input);
            qDebug() << "Nii2Array successfully";

            QByteArray temp;
            int randArrSize = (int)(percent*imgHeight*imgWidth*imgDepth)/100;

            for(int i=0; i<randArrSize; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
                int ind = randPxls[i];
                scatImageArr[ind] = imageArr[ind];

                char buf[9];
                ::sprintf(buf, "%d", ind);
                temp.append(buf);
                temp.append("\n");
            }
            file.write(temp);
            file.flush();
            file.close();

    //        qDebug() << "So far successfully";

            // Writing scattered nii file
            array_to_nii(nim_input, scatImageArr);
            QString nameStr = "./img-outputs/masks/";
            nameStr = nameStr + ui->lineEdit_randMask_fileName->text() +  "_" + (ui->lineEdit_randMask_percentage->text()) + "randMask" + ".nii";
            //***************************************************
            //Convert QString to char*
            QByteArray baTemp = nameStr.toLocal8Bit();
            const char *fout = baTemp.data();
            //***************************************************
            nifti_set_filenames(nim_input,fout,1,1);
            nifti_image_write(nim_input);
            //***************************************************

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageMaskObject, imgWidth, imgHeight);

            imageMask = QPixmap::fromImage(*imageMaskObject);
            ui->label_randMaks_ProcImg->setPixmap(imageMask.scaled(ui->label_randMaks_ProcImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
    //        ui->label_randMaks_ProcImg->setPixmap(imageMask);

            if(imgDepth != 1) {
                ui->horizontalSlider_radnMaskProcImg->setVisible(true);
            }

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;
            delete [] randPxls;
            randPxls = NULL;

            QMessageBox::information(this, "Saved", "The code has run succesfully!");
        }
    }

    qDebug() << "So far successfully";
}

void RandomMask::on_pushButton_randMask_fileNameSave_clicked()
{
    bool isSaved;

    if(ui->label_randMaks_ProcImg->pixmap()==0 //.isNull()
            || ui->lineEdit_randMask_fileName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
    } else {
        QDir dir("./img-outputs/masks");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/masks/";
        nameStr = nameStr + ui->lineEdit_randMask_fileName->text() +  "_" + (ui->lineEdit_randMask_percentage->text()) + "percentage_of_random_pixels" + ".png";
        isSaved = imageMaskObject->save(nameStr, 0, -1);
    }
    if(isSaved) {
        QMessageBox::information(this, "Saved", "The image is saved succesfully!");
    } else {
        QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
    }
}
