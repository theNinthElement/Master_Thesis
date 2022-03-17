#include "slicewiseeeddialog.h"
#include "ui_slicewiseeeddialog.h"
#include "inpainting.h"
#include "supplementary_functions.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>

#include <QDebug>
#include "derivatives.h"

SlicewiseEEDDialog::SlicewiseEEDDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::SlicewiseEEDDialog)
{
    ui->setupUi(this);

    // The MSE tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_mse_tol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The Tolerance lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_2deed_tol->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The TimeStepSize lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_2deed_timeStep->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The FED Inner Cycle lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_innerCycleSize->setValidator(new QIntValidator(0, 999999999, this));
    // The FSI Inner Loop lineedit will only accept integers between 0 and 999999999
    ui->lineEdit_innerLoopSIze->setValidator(new QIntValidator(0, 999999999, this));

    //Make the displayed images center of the label
    ui->label_imgShow->setAlignment(Qt::AlignCenter);
    ui->label_maskShow->setAlignment(Qt::AlignCenter);

    //Make FSI is chosen by default and hidden MRE tolerance &  line edits
    ui->radioButton_FSIScheme->setChecked(true);
    ui->label_mse_tol->setVisible(false);
    ui->lineEdit_mse_tol->setVisible(false);
    ui->label_innerCycleSize->setVisible(false);
    ui->lineEdit_innerCycleSize->setVisible(false);

    //Set default value
    ui->lineEdit_2deed_timeStep->setPlaceholderText("0.25");
    ui->lineEdit_saveImgName->setPlaceholderText("inpaintedImage");
    ui->lineEdit_innerLoopSIze->setPlaceholderText("20");
    ui->lineEdit_innerCycleSize->setPlaceholderText("20");
    ui->lineEdit_2deed_tol->setPlaceholderText("0.00001");
    ui->lineEdit_mse_tol->setPlaceholderText("100");
}

SlicewiseEEDDialog::~SlicewiseEEDDialog()
{
    delete ui;
}

void SlicewiseEEDDialog::on_radioButton_2deed_ExScheme_clicked()
{
    //Hide the Number of Inner Loop steps and MSE tolerance line edits and labels
    ui->lineEdit_innerLoopSIze->setVisible(false);
    ui->label_innerLoopSIze->setVisible(false);
    ui->lineEdit_innerCycleSize->setVisible(false);
    ui->label_innerCycleSize->setVisible(false);
    ui->lineEdit_mse_tol->setVisible(false);
    ui->label_mse_tol->setVisible(false);

    ui->lineEdit_2deed_tol->setVisible(true);
    ui->label_2deed_tol->setVisible(true);
}

void SlicewiseEEDDialog::on_radioButton_FEDScheme_clicked()
{
    ui->lineEdit_2deed_tol->setVisible(false);
    ui->label_2deed_tol->setVisible(false);
    ui->lineEdit_innerLoopSIze->setVisible(false);
    ui->label_innerLoopSIze->setVisible(false);

    ui->lineEdit_mse_tol->setVisible(true);
    ui->label_mse_tol->setVisible(true);
    ui->lineEdit_innerCycleSize->setVisible(true);
    ui->label_innerCycleSize->setVisible(true);
}

void SlicewiseEEDDialog::on_radioButton_FSIScheme_clicked()
{
    ui->lineEdit_innerCycleSize->setVisible(false);
    ui->label_innerCycleSize->setVisible(false);
    ui->lineEdit_mse_tol->setVisible(false);
    ui->label_mse_tol->setVisible(false);

    ui->lineEdit_2deed_tol->setVisible(true);
    ui->label_2deed_tol->setVisible(true);
    ui->lineEdit_innerLoopSIze->setVisible(true);
    ui->label_innerLoopSIze->setVisible(true);
}

void SlicewiseEEDDialog::on_pushButton_loadImg_clicked()
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

    ui->spinBox_2dEED->setRange(1,imgDepth);

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input);
    // Create image and set to the label
    imageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*imageObject);
    ui->label_imgShow->setPixmap(image.scaled(ui->label_imgShow->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void SlicewiseEEDDialog::on_pushButton_mask_clicked()
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
//    int dataDim = ui->comboBox_randMask->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_mask);
    // Create image and set to the label
    maskImageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*maskImageObject);
    ui->label_maskShow->setPixmap(image.scaled(ui->label_maskShow->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void SlicewiseEEDDialog::on_pushButton_run2dEED_clicked()
{
    int innerLoopSize, innerCycleSize;
    float timeStep, tol;

    if(ui->radioButton_2deed_ExScheme->isChecked()) {
        if(ui->lineEdit_2deed_tol->text().isEmpty() || ui->lineEdit_2deed_timeStep->text().isEmpty() || ui->label_imgShow->pixmap()==0 //.isNull()
                || ui->label_maskShow->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter (time) step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_2deed_timeStep->text()).toDouble();
            tol = (ui->lineEdit_2deed_tol->text()).toDouble();

            double *imageArr, *scatImageArr;
            double **imageSliceArr, **scatImageSliceArr;
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

//            imageArr = nii_to_array(nim_input);
            scatImageArr = nii_to_array(nim_input_mask);

            imageSliceArr = nii_to_arrayof_2d_arrays(nim_input);
            scatImageSliceArr = nii_to_arrayof_2d_arrays(nim_input_mask);


//            qDebug() << imageSliceArr[40][5351] << imageSliceArr[40][5353] << imageSliceArr[40][5354] << imageSliceArr[40][5356] << imageSliceArr[40][5451] << imageSliceArr[40][5561];

            if(ui->checkBox_monitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for Explicit Scheme has not been implemented yet :(");
                return;
            } else {
                int leng = randPxlStrArr.size()-1;
                int** randSlicePxls = linear2twoDim(randPxls,leng,imgWidth,imgHeight,imgDepth);

                int sliceNum = 49;
                eed_inpainting(tol, timeStep, scatImageSliceArr[sliceNum], imageSliceArr[sliceNum], randSlicePxls[sliceNum], imgWidth, imgHeight);

//                double* mseArr = new double[sizeof(double) * imgDepth];
//                double* aaeArr = new double[sizeof(double) * imgDepth];
//                for(int i=81; i<84; i++) {
//                    eed_inpainting(tol, timeStep, scatImageSliceArr[i], imageSliceArr[i], randSlicePxls[i], imgWidth, imgHeight);

//                    mseArr[i] = mse(imageSliceArr[i],scatImageSliceArr[i],imgWidth*imgHeight);
//                    aaeArr[i] = aae(imageSliceArr[i],scatImageSliceArr[i],imgWidth*imgHeight);
//                }
//                for(int i=81; i<84; i++) {
//                    qDebug() << i << mseArr[i] << aaeArr[i];
//                }
            }

            array_to_nii(nim_input,scatImageSliceArr[0]);
            arrayof_2d_arrays_to_nii(nim_input,scatImageSliceArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_imgShow->setPixmap(image.scaled(ui->label_imgShow->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageSliceArr;
            imageArr = NULL;
            delete [] scatImageSliceArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    } else if(ui->radioButton_FEDScheme->isChecked()) {
        qDebug() << "FED Scheme\n";
    } else if(ui->radioButton_FSIScheme->isChecked()) {
        qDebug() << "FSI Scheme\n";

        if(ui->lineEdit_2deed_timeStep->text().isEmpty() || ui->lineEdit_2deed_tol->text().isEmpty() || ui->lineEdit_innerLoopSIze->text().isEmpty() || ui->label_imgShow->pixmap()==0 //.isNull()
                || ui->label_maskShow->pixmap()==0) //.isNull())
        {
            QMessageBox::warning(this, "Missing Input", "Please enter step size and/or tol and/or image and/or mask image and/or mask txt file");
        } else {
            timeStep = (ui->lineEdit_2deed_timeStep->text()).toDouble();
            tol = (ui->lineEdit_2deed_tol->text()).toDouble();
            innerLoopSize = (ui->lineEdit_innerLoopSIze->text()).toDouble();

            double *imageArr, *scatImageArr;
            double **imageSliceArr, **scatImageSliceArr;
            int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];
            float gridSpcX = nim_input->dx, gridSpcY = nim_input->dy, gridSpcZ = nim_input->dz, gridSpcT = nim_input->dt;

            // Decleartion and memory allocation
            int* randPxls = new int[sizeof(int) * (randPxlStrArr.size()-1)];

            for(int i=0; i<randPxlStrArr.size()-1; i++) {
                QString temp = randPxlStrArr[i];
                randPxls[i] = temp.toInt();
            }
            imageArr = nii_to_array(nim_input);
//            scatImageArr = nii_to_array(nim_input_mask);

            imageSliceArr = nii_to_arrayof_2d_arrays(nim_input);
            scatImageSliceArr = nii_to_arrayof_2d_arrays(nim_input_mask);

            if(ui->checkBox_monitoring->isChecked()) {
                QMessageBox::information(this, "Monitoring", "Monitoring for FSI Scheme has not been implemented yet :(");
                return;
            } else {
                int leng = randPxlStrArr.size()-1;
                int** randSlicePxls = linear2twoDim(randPxls,leng,imgWidth,imgHeight,imgDepth);

//                int sliceNum = 80;
//                eed_inpainting_FSI(tol, timeStep, scatImageSliceArr[sliceNum], imageSliceArr[sliceNum], randSlicePxls[sliceNum], imgWidth, imgHeight);

                double* mseArr = new double[sizeof(double) * imgDepth];
                double* aaeArr = new double[sizeof(double) * imgDepth];
                for(int i=0; i<26; i++) {
                    eed_inpainting_FSI(tol, timeStep, innerLoopSize, scatImageSliceArr[i], imageSliceArr[i], randSlicePxls[i], imgWidth, imgHeight);

                    mseArr[i] = mse(imageSliceArr[i],scatImageSliceArr[i],imgWidth*imgHeight);
                    aaeArr[i] = aae(imageSliceArr[i],scatImageSliceArr[i],imgWidth*imgHeight);
                }
                for(int i=0; i<26; i++) {
                    qDebug() << i << mseArr[i] << aaeArr[i];
                }
            }

            array_to_nii(nim_input,scatImageSliceArr[0]);
            arrayof_2d_arrays_to_nii(nim_input,scatImageSliceArr);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageObject, imgWidth, imgHeight);
            image = QPixmap::fromImage(*imageObject);
            ui->label_imgShow->setPixmap(image.scaled(ui->label_imgShow->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//            ui->label_eedInpaintImg->setPixmap(image);

            delete [] imageSliceArr;
            imageArr = NULL;
            delete [] scatImageSliceArr;
            scatImageArr = NULL;

            QMessageBox::information(this, "Finished running", "The code has run succesfully!");
        }
    }
}

void SlicewiseEEDDialog::on_pushButton_saveImgName_clicked()
{
    bool isSaved;

    if(ui->label_imgShow->pixmap()==0 //.isNull()
            || ui->lineEdit_saveImgName->text().isEmpty()) {
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

        if(ui->radioButton_FSIScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_saveImgName->text() + "FSI" + "_innerLoopSize" + (ui->lineEdit_innerLoopSIze->text()) +  "_timeStepSize" + (ui->lineEdit_2deed_timeStep->text()) + "_tol" + (ui->lineEdit_2deed_tol->text()) + "_" + strPercentage + "percent_of_pixels" + ".png";
        } else if(ui->radioButton_2deed_ExScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_saveImgName->text() + "_timeStepSize" + (ui->lineEdit_2deed_timeStep->text()) + "_tol" + (ui->lineEdit_2deed_tol->text()) + "_" + strPercentage + "percent_of_pixels" + ".png";
        } else if(ui->radioButton_FEDScheme->isChecked()) {
            nameStr = nameStr + ui->lineEdit_saveImgName->text() + "FED" + "_innerCycleSize" + (ui->lineEdit_innerCycleSize->text()) + "_timeStepSize" + (ui->lineEdit_2deed_timeStep->text()) + "_tol" + (ui->lineEdit_mse_tol->text()) + "_" + strPercentage + "percent_of_pixels" + ".png";
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
