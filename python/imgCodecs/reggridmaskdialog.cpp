#include "reggridmaskdialog.h"
#include "ui_reggridmaskdialog.h"
#include "supplementary_functions.h"
#include "masks.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QTextStream>

#include <QDebug>

RegGridMaskDialog::RegGridMaskDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::RegGridMaskDialog)
{
    ui->setupUi(this);

    // The Percentage lineedit will only accept integers between 1 and 100
    ui->lineEdit_regMask_percentage->setValidator(new QIntValidator(1, 100, this));

    //Make the displayed image center of the label
    ui->label_regGridMaskProcImg->setAlignment(Qt::AlignCenter);
    //Make the displayed image center of the label
    ui->label_regMaskOrigImg->setAlignment(Qt::AlignCenter);

    //Set default value
    ui->lineEdit_regMask_percentage->setPlaceholderText("5");
    ui->lineEdit_regMaskImgName->setPlaceholderText("maskImg");

    //Make sliders invisible by default
    ui->horizontalSlider_regGridMaskProcImg->setVisible(false);
    ui->horizontalSlider_regGridOrigImg->setVisible(false);

    //Hide comboBox as of now
    ui->comboBox_regMask_ImgDim->setVisible(false);
}

RegGridMaskDialog::~RegGridMaskDialog()
{
    delete ui;
}

void RegGridMaskDialog::on_pushButton_regMask_imgUpl_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_grid_mask = nifti_image_read(fin, 1);                                                               // copy image to mask image
    nim_input_reg = nifti_image_read(fin, 1);                                          // Read for using later to get smaller size mask
    if(!nim_input_reg) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_input_reg->dim[1];
    int imgHeight = nim_input_reg->dim[2];
    int imgDepth = nim_input_reg->dim[3];
    int imgTimeLen = nim_input_reg->dim[4];
//    int dataDim = ui->comboBox_regMask_ImgDim->currentIndex();

//    float gridSpcX = nim_input_reg->dx, gridSpcY = nim_input_reg->dy, gridSpcZ = nim_input_reg->dz, gridSpcT = nim_input_reg->dt;

//    qDebug() << "anisotropic data gridSpcX : "<< gridSpcX << "   gridSpcY : "<< gridSpcY << "   gridSpcZ : "<< gridSpcZ << "   gridSpcT : "<< gridSpcT;

    if(imgDepth == 1) {
        ui->horizontalSlider_regGridMaskProcImg->setVisible(false);
        ui->horizontalSlider_regGridOrigImg->setVisible(false);
    }

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_reg);
    // Create image and set to the label
    imageRegMaskObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
//    if(imgTimeLen != 1) {
//        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
//        return;
//    }
    imageRegMask = QPixmap::fromImage(*imageRegMaskObject);
    ui->label_regMaskOrigImg->setPixmap(imageRegMask.scaled(ui->label_regMaskOrigImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);

    if(imgDepth != 1) {
        ui->horizontalSlider_regGridOrigImg->setVisible(true);
    }
}

void RegGridMaskDialog::on_pushButton_regMask_run_clicked()
{
    int imgWidth, imgHeight, imgDepth, imgTimeLen, ratio;

    if(ui->lineEdit_regMask_percentage->text().isEmpty() || ui->label_regMaskOrigImg->pixmap()==0) //.isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter percentage and/or load an image");
    } else {
        int regGridMaskArrLen = 0;                                        // It is needed for the reg.grd. selection by sampling size not by a ratio!
        double *imageArr, *regGridPxls;

        ratio = (ui->lineEdit_regMask_percentage->text()).toInt();
        imgWidth = nim_input_reg->dim[1];
        imgHeight = nim_input_reg->dim[2];
        imgDepth = nim_input_reg->dim[3];
        imgTimeLen = nim_input_reg->dim[4];

        //Check if image data is a 4D data
        if(imgTimeLen != 1) {
//            double** imageArr = nii4d_to_array(nim_grid_mask);

            // Decleration and Memory allocation
            double** scatImageArr = new double*[sizeof(double)*imgTimeLen];
            for(int i = 0; i < imgTimeLen; ++i) {
                scatImageArr[i] = new double[sizeof(double)*imgWidth*imgHeight*imgDepth];
            }

            //Initialize mask image with zeros  *************************************************************************************
            for(int i = 0; i < imgTimeLen; ++i) {
                for(int j=0; j<imgWidth*imgHeight*imgDepth; j++) {
                    scatImageArr[i][j] = 0;
                }
            }
            double ** binRegGrid;
            binRegGrid = regGridPixelsVolum_4D(ratio, imgWidth, imgHeight, imgDepth, imgTimeLen);                  // 4D image Volumetric regular grid by sampling size
            array_to_nii4d(nim_input_reg,binRegGrid);

            // Writing scattered nii
//            array_to_nii(nim_input_reg, scatImageArr);
            QString nameStr = "./img-outputs/masks/";
            nameStr = nameStr + ui->lineEdit_regMaskImgName->text() +  "_" + (ui->lineEdit_regMask_percentage->text()) + "regGridMask_4D" + ".nii";
            //***************************************************
            //Convert QString to char*
            QByteArray baTemp = nameStr.toLocal8Bit();
            const char *fout = baTemp.data();
            //***************************************************
            nifti_set_filenames(nim_input_reg, fout, 1, 1);
            nifti_image_write(nim_input_reg);

//            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
//            array_to_png(scatImageArr, imageRegMaskObject, imgWidth, imgHeight);

//            imageRegMask = QPixmap::fromImage(*imageRegMaskObject);
//            ui->label_regGridMaskProcImg->setPixmap(imageRegMask.scaled(ui->label_regGridMaskProcImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//    //        ui->label_randMaks_ProcImg->setPixmap(imageMask);

            if(imgDepth != 1) {
                ui->horizontalSlider_regGridMaskProcImg->setVisible(true);
            }

            delete [] scatImageArr;
            scatImageArr = NULL;
            delete [] binRegGrid;
            binRegGrid = NULL;

            QMessageBox::information(this, "Saved", "The code has run succesfully!");

        } else {
            // Decleration and Memory allocation
            double* scatImageArr = new double[sizeof(double) * imgWidth*imgHeight*imgDepth];

            //Initialize scarImageArr to the black(=0) image
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                scatImageArr[i] = 0;
            }
    //        regGridPxls = regGridPixels(ratio, imgWidth*imgHeight*imgDepth);                                              // Regular grid by ratio
    //        regGridPxls = regGridPixelsSlicewise(ratio, imgWidth, imgHeight, imgDepth, regGridMaskArrLen);                  // Slicewise regular grid by sampling size
            regGridPxls = regGridPixelsVolum(ratio, imgWidth, imgHeight, imgDepth, regGridMaskArrLen);                  // Volumetric regular grid by sampling size

            QDir dir("./img-outputs/masks");
            if (!dir.exists())
                dir.mkpath(".");
            QFile file("./img-outputs/masks/pixel_locations");
            if(!file.open(QFile::WriteOnly | QFile::Text)) {
                QMessageBox::warning(this, "Error", "File is NOT open");
            }

            imageArr = nii_to_array(nim_input_reg);

    //        int len = (int)((imgWidth*imgHeight*imgDepth)/ratio);
    //        int divisor = 1;
    //        for (int i=2; i<len; i++) {
    //            if (len%i==0) {
    //                divisor = i;
    //                break;
    //            }
    //        }
    //        qDebug() << divisor;
    //        qDebug() << nim_grid_mask->dim[1] << nim_grid_mask->dim[2] << nim_grid_mask->dim[3] << nim_grid_mask->dim[4];
    //        int dim1 = nim_grid_mask->dim[1]/(pow(ratio,1.0/3));
    //        int dim2 = nim_grid_mask->dim[2]/(pow(ratio,1.0/3));
    //        int dim3 = nim_grid_mask->dim[3]/(pow(ratio,1.0/3));
    //        int len = (imgWidth*imgHeight*imgDepth)/ratio;
    //        qDebug() << dim1 << dim2 << dim3 << len;
    //        dim3 = len/(dim1*dim2)+1;
    //        qDebug() << dim3;
    //        nim_grid_mask->dim[1] = dim1;
    //        nim_grid_mask->dim[2] = dim2;
    //        nim_grid_mask->dim[3] = dim3;
    //        nim_grid_mask->dim[1] = divisor;
    //        nim_grid_mask->dim[2] = len/divisor;
    //        nim_grid_mask->dim[3] = 1;
    //        nim_grid_mask->dim[0] = 2;
    //        nifti_update_dims_from_array(nim_grid_mask);   // changing according sizes nx, ny, nz etc.
    //        qDebug() << nim_grid_mask->dim[0] << nim_grid_mask->dim[1] << nim_grid_mask->dim[2] << nim_grid_mask->dim[3] << nim_grid_mask->dim[4];
    //        double* maskImgArr = new double[sizeof(double) * (dim1*dim2*dim3)];



            QByteArray temp;
    //        for(int i=0; i<((imgHeight*imgWidth*imgDepth)/ratio); i++) {                                 // Used reg.grid.arr. by ratio length in the for loop
            for(int i=0; i<regGridMaskArrLen; i++) {                                                       // Used reg.grid.arr. by sampling size length in the for loop
                int ind = regGridPxls[i];
                scatImageArr[ind] = imageArr[ind];
    //            maskImgArr[i] = imageArr[ind];

                char buf[10];
                ::sprintf(buf, "%d", ind);
                temp.append(buf);
                temp.append("\n");
            }
            file.write(temp);
            file.flush();
            file.close();


    //        qDebug() << nim_input_reg->scl_slope << nim_input_reg->scl_inter;
    //        for(int i=0; i<((imgHeight*imgWidth*imgDepth)/ratio); i++) {
    //            imageArr[i] = nim_input_reg->scl_slope*imageArr[i] + nim_input_reg->scl_inter;
    //        }



    //        for(int i=((imgHeight*imgWidth*imgDepth)/ratio); i<dim1*dim2*dim3; i++) {                                 // Filling mask image with zeros
    //            maskImgArr[i] = 0;
    //        }
    //        // Writing mask nii
    //        array_to_nii(nim_grid_mask, maskImgArr);
    //        QString nameStrMask = "./img-outputs/masks/reg-grid_for_tests";
    //        nameStrMask = nameStrMask + ui->lineEdit_regMaskImgName->text() +  "_" + (ui->lineEdit_regMask_percentage->text()) + "regGridMask" + ".nii";
    //        //***************************************************
    //        //Convert QString to char*
    //        QByteArray baTempMask = nameStrMask.toLocal8Bit();
    //        const char *foutMask = baTempMask.data();
    //        //***************************************************
    //        nifti_set_filenames(nim_grid_mask, foutMask, 1, 1);
    //        nifti_image_write(nim_grid_mask);


            // Writing scattered nii
            array_to_nii(nim_input_reg, scatImageArr);
            QString nameStr = "./img-outputs/masks/";
            nameStr = nameStr + ui->lineEdit_regMaskImgName->text() +  "_" + (ui->lineEdit_regMask_percentage->text()) + "regGridMask" + ".nii";
            //***************************************************
            //Convert QString to char*
            QByteArray baTemp = nameStr.toLocal8Bit();
            const char *fout = baTemp.data();
            //***************************************************
            nifti_set_filenames(nim_input_reg, fout, 1, 1);
            nifti_image_write(nim_input_reg);

            scaleToUchar(scatImageArr, imgWidth*imgHeight*imgDepth);
            array_to_png(scatImageArr, imageRegMaskObject, imgWidth, imgHeight);

            imageRegMask = QPixmap::fromImage(*imageRegMaskObject);
            ui->label_regGridMaskProcImg->setPixmap(imageRegMask.scaled(ui->label_regGridMaskProcImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
    //        ui->label_randMaks_ProcImg->setPixmap(imageMask);

            if(imgDepth != 1) {
                ui->horizontalSlider_regGridMaskProcImg->setVisible(true);
            }

            delete [] imageArr;
            imageArr = NULL;
            delete [] scatImageArr;
            scatImageArr = NULL;
            delete [] regGridPxls;
            regGridPxls = NULL;

            QMessageBox::information(this, "Saved", "The code has run succesfully!");
        }
    }
}

void RegGridMaskDialog::on_pushButton_regMaskImgSave_clicked()
{
    bool isSaved;
    int imgTimeLen = nim_input_reg->dim[4];

    if(imgTimeLen != 1) {
        if(ui->lineEdit_regMaskImgName->text().isEmpty()) {
            QMessageBox::warning(this, "Missing Input", "No file name is given");
        } else {
            QDir dir("./img-outputs/masks");
            if (!dir.exists())
                dir.mkpath(".");

            QString nameStr = "./img-outputs/masks/";
            nameStr = nameStr + ui->lineEdit_regMaskImgName->text() +  "_" + (ui->lineEdit_regMask_percentage->text()) + "ratio_of_all_pixels" + ".png";
            isSaved = imageRegMaskObject->save(nameStr, 0, -1);
        }
        if(isSaved) {
            QMessageBox::information(this, "Saved", "The image is saved succesfully!");
        } else {
            QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
        }
    } else {
        if(ui->label_regGridMaskProcImg->pixmap()==0 || ui->lineEdit_regMaskImgName->text().isEmpty()) {
            QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
        } else {
            QDir dir("./img-outputs/masks");
            if (!dir.exists())
                dir.mkpath(".");

            QString nameStr = "./img-outputs/masks/";
            nameStr = nameStr + ui->lineEdit_regMaskImgName->text() +  "_" + (ui->lineEdit_regMask_percentage->text()) + "ratio_of_all_pixels" + ".png";
            isSaved = imageRegMaskObject->save(nameStr, 0, -1);
        }
        if(isSaved) {
            QMessageBox::information(this, "Saved", "The image is saved succesfully!");
        } else {
            QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
        }
    }

}
