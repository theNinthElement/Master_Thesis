#include "directneighbordilationdialog.h"
#include "ui_directneighbordilationdialog.h"
#include "supplementary_functions.h"
#include "masks.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <QTextStream>

#include <QDebug>

DirectNeighborDilationDialog::DirectNeighborDilationDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::DirectNeighborDilationDialog)
{
    ui->setupUi(this);

    //Make the displayed image center of the label
    ui->label_directDilationOrigImg->setAlignment(Qt::AlignCenter);
    //Make the displayed image center of the label
    ui->label_directDilationMaskImg->setAlignment(Qt::AlignCenter);
    //Make the displayed image center of the label
    ui->label_directDilationRecImg->setAlignment(Qt::AlignCenter);

    //Set default value
    ui->lineEdit_directDilationImgName->setPlaceholderText("updatedMaskImg");

    //Hide comboBox as of now
    ui->comboBox_directDilation_ImgDim->setVisible(false);
}

DirectNeighborDilationDialog::~DirectNeighborDilationDialog()
{
    delete ui;
}

void DirectNeighborDilationDialog::on_pushButton_directDilation_imgUpl_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_input_orig = nifti_image_read(fin, 1);
    if(!nim_input_orig) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_input_orig->dim[1];
    int imgHeight = nim_input_orig->dim[2];
//    int imgDepth = nim_input_orig->dim[3];
    int imgTimeLen = nim_input_orig->dim[4];
//    int dataDim = ui->comboBox_regMask_ImgDim->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input_orig);
    // Create image and set to the label
    imageDirDilOrigImgObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageDirDilOrigImg = QPixmap::fromImage(*imageDirDilOrigImgObject);
    ui->label_directDilationOrigImg->setPixmap(imageDirDilOrigImg.scaled(ui->label_directDilationOrigImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);
}

void DirectNeighborDilationDialog::on_pushButton_directDilation_mskUpl_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_grid_mask = nifti_image_read(fin, 1);
    if(!nim_grid_mask) {
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
        QMessageBox::warning(this, "Error", "File is NOT open.");
    }
    QString arr= file.readAll();
    randPxlStrArr = arr.split('\n');
    file.close();
    //***************************************************
    int imgWidth = nim_grid_mask->dim[1];
    int imgHeight = nim_grid_mask->dim[2];
//    int imgDepth = nim_grid_mask->dim[3];
    int imgTimeLen = nim_grid_mask->dim[4];
//    int dataDim = ui->comboBox_regMask_ImgDim->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_grid_mask);
    // Create image and set to the label
    imageDirDilMaskObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageDirDilMask = QPixmap::fromImage(*imageDirDilMaskObject);
    ui->label_directDilationMaskImg->setPixmap(imageDirDilMask.scaled(ui->label_directDilationMaskImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);
}

void DirectNeighborDilationDialog::on_pushButton_directDilation_recUpl_clicked()
{
    QString imageDataPath = QFileDialog::getOpenFileName(this,tr("Open File"),"",tr("NIfTI (*.nii)" ));

    //***************************************************
    //Convert QString to char*
    QByteArray baTemp = imageDataPath.toLocal8Bit();
    const char *fin = baTemp.data();
    //***************************************************
    // Read input dataset, including data
    nim_grid_recons = nifti_image_read(fin, 1);
    if(!nim_grid_recons) {
        fprintf(stderr,"Failed to read NIfTI image data from '%s'\n", fin);
        return;
    }
    //***************************************************
    int imgWidth = nim_grid_recons->dim[1];
    int imgHeight = nim_grid_recons->dim[2];
//    int imgDepth = nim_grid_recons->dim[3];
    int imgTimeLen = nim_grid_recons->dim[4];
//    int dataDim = ui->comboBox_regMask_ImgDim->currentIndex();

    unsigned char* niiImgArr = nii_to_ucharArray(nim_grid_recons);
    // Create image and set to the label
    imageDirDilRecImgObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    imageDirDilRecImg = QPixmap::fromImage(*imageDirDilRecImgObject);
    ui->label_directDilationRecImg->setPixmap(imageDirDilRecImg.scaled(ui->label_directDilationRecImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_randMask_origImg->setPixmap(imageMask);
}

void DirectNeighborDilationDialog::on_pushButton_directDilationImgSave_clicked()
{
    bool isSaved;

    if(ui->label_directDilationOrigImg->pixmap()==0 //.isNull()
            || ui->label_directDilationMaskImg->pixmap()==0 //.isNull()
            || ui->lineEdit_directDilationImgName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
    } else {
        QDir dir("./img-outputs/masks");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/masks/";
        nameStr = nameStr + ui->lineEdit_directDilationImgName->text() +  "_" + "dilated" + ".png";
        isSaved = imageDirDilMaskObject->save(nameStr, 0, -1);
    }
    if(isSaved) {
        QMessageBox::information(this, "Saved", "The image is saved succesfully!");
    } else {
        QMessageBox::warning(this, "Saving Problem", "The image is NOT saved");
    }
}

void DirectNeighborDilationDialog::on_pushButton_directDilation_run_clicked()
{
    int imgWidth, imgHeight, imgDepth, imgTimeLen;
    bool excludImg = false;

    if(ui->label_directDilationOrigImg->pixmap()==0 //.isNull()
            || ui->label_directDilationMaskImg->pixmap()==0) //.isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter percentage and/or load an image");
    } else {
        int dilatedMaskLen = 0;
        double *imageArr, *recImageArr, *mskImageArr, *dilatedPxls;

        imgWidth = nim_grid_mask->dim[1];
        imgHeight = nim_grid_mask->dim[2];
        imgDepth = nim_grid_mask->dim[3];
        imgTimeLen = nim_grid_mask->dim[4];

        //Check if image data is a 4D data
        if(imgTimeLen != 1) {
            QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
            return;
        }
        // Decleartion and memory allocation
        int* regGridPxls = new int[sizeof(int) * (randPxlStrArr.size()-1)];
        for(int i=0; i<randPxlStrArr.size()-1; i++) {
            QString temp = randPxlStrArr[i];
            regGridPxls[i] = temp.toInt();
        }
        imageArr = nii_to_array(nim_input_orig);
        recImageArr = nii_to_array(nim_grid_recons);
        mskImageArr = nii_to_array(nim_grid_mask);

//        dilatedPxls = dilatedMask2D(imgWidth, imgHeight, imgDepth, regGridPxls, dilatedMaskLen, excludImg);
        dilatedPxls = dilatedMask3D(imgWidth, imgHeight, imgDepth, regGridPxls, dilatedMaskLen, excludImg);
//        dilatedPxls = dilatedMask3D_partial(imgWidth, imgHeight, imgDepth, regGridPxls, dilatedMaskLen, excludImg, 1);
        qDebug() << "Dilated Mask Size: " << dilatedMaskLen;

        QDir dir("./img-outputs/masks");
        if (!dir.exists())
            dir.mkpath(".");
        QFile file("./img-outputs/masks/pixel_locations_dilated");
        if(!file.open(QFile::WriteOnly | QFile::Text)) {
            QMessageBox::warning(this, "Error", "File is NOT open");
        }
        QByteArray temp;

        if(excludImg == false) {
            //Initialize maskImageArr to the black(=0) image
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                mskImageArr[i] = 0;
            }
        } else {
            //Initialize maskImageArr to the reconstructed image
            for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
                mskImageArr[i] = recImageArr[i];
            }
        }

        for(int i=0; i<dilatedMaskLen; i++) {
            int ind = (int)dilatedPxls[i];                              // For get a mask with direct neighbors
            mskImageArr[ind] = imageArr[ind];                              // For lossless mask pixel values

            char buf[9];
            ::sprintf(buf, "%d", ind);
            temp.append(buf);
            temp.append("\n");
        }
        file.write(temp);
        file.flush();
        file.close();

        // Writing mask nii
        array_to_nii(nim_grid_mask, mskImageArr);
        QString nameStrMask = "./img-outputs/masks/_";
        nameStrMask = nameStrMask + ui->lineEdit_directDilationImgName->text() + "regGridDilated" + ".nii";
        //***************************************************
        //Convert QString to char*
        QByteArray baTempMask = nameStrMask.toLocal8Bit();
        const char *foutMask = baTempMask.data();
        //***************************************************
        nifti_set_filenames(nim_grid_mask, foutMask, 1, 1);
        nifti_image_write(nim_grid_mask);

        scaleToUchar(mskImageArr, imgWidth*imgHeight*imgDepth);

        imageDirDilMask = QPixmap::fromImage(*imageDirDilMaskObject);
        ui->label_directDilationMaskImg->setPixmap(imageDirDilMask.scaled(ui->label_directDilationMaskImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
//        ui->label_randMaks_ProcImg->setPixmap(imageMask);

        delete [] imageArr;
        imageArr = NULL;
        delete [] mskImageArr;
        mskImageArr = NULL;
        delete [] recImageArr;
        recImageArr = NULL;
        delete [] regGridPxls;
        regGridPxls = NULL;

//                // For Test: Decleration and Memory allocation
//                double* testImg = new double[sizeof(double) * 4*4*4];
//                double* testRecImg = new double[sizeof(double) * 4*4*4];
//                testImg[0] = 10;
//                testImg[1] = 8;
//                testImg[2] = 6;
//                testImg[3] = 3;
//                testImg[4] = 7;
//                testImg[5] = 9;
//                testImg[6] = 6;
//                testImg[7] = 7;
//                testImg[8] = 1;
//                testImg[9] = 0;
//                testImg[10] = 0;
//                testImg[11] = 0;
//                testImg[12] = 0;
//                testImg[13] = 0;
//                testImg[14] = 2;
//                testImg[15] = 1;
//                testImg[16] = 2;
//                testImg[17] = 3;
//                testImg[18] = 4;
//                testImg[19] = 5;
//                testImg[20] = 2;
//                testImg[21] = 3;
//                testImg[22] = 4;
//                testImg[23] = 5;
//                testImg[24] = 5;
//                testImg[25] = 1;
//                testImg[26] = 2;
//                testImg[27] = 1;
//                testImg[28] = 2;
//                testImg[29] = 6;
//                testImg[30] = 7;
//                testImg[31] = 8;
//                testImg[32] = 9;
//                testImg[33] = 10;
//                testImg[34] = 11;
//                testImg[35] = 14;
//                testImg[36] = 18;
//                testImg[37] = 16;
//                testImg[38] = 13;
//                testImg[39] = 17;
//                testImg[40] = 19;
//                testImg[41] = 16;
//                testImg[42] = 17;
//                testImg[43] = 11;
//                testImg[44] = 10;
//                testImg[45] = 10;
//                testImg[46] = 10;
//                testImg[47] = 10;
//                testImg[48] = 20;
//                testImg[49] = 22;
//                testImg[50] = 21;
//                testImg[51] = 12;
//                testImg[52] = 13;
//                testImg[53] = 14;
//                testImg[54] = 15;
//                testImg[55] = 12;
//                testImg[56] = 13;
//                testImg[57] = 14;
//                testImg[58] = 15;
//                testImg[59] = 15;
//                testImg[60] = 11;
//                testImg[61] = 12;
//                testImg[62] = 11;
//                testImg[63] = 12;
//                testImg[64] = 10;
//                testImg[65] = 8;
//                testImg[66] = 6;
//                testImg[67] = 3;
//                testImg[68] = 7;
//                testImg[69] = 9;
//                testImg[70] = 6;
//                testImg[71] = 7;
//                testImg[72] = 1;
//                testImg[73] = 0;
//                testImg[74] = 0;
//                testImg[75] = 0;
//                testImg[76] = 0;
//                testImg[77] = 0;
//                testImg[78] = 2;
//                testImg[79] = 1;
//                testImg[80] = 2;
//                testImg[81] = 3;
//                testImg[82] = 4;
//                testImg[83] = 5;
//                testImg[84] = 2;
//                testImg[85] = 3;
//                testImg[86] = 4;
//                testImg[87] = 5;
//                testImg[88] = 5;
//                testImg[89] = 1;
//                testImg[90] = 2;
//                testImg[91] = 1;
//                testImg[92] = 2;
//                testImg[93] = 6;
//                testImg[94] = 7;
//                testImg[95] = 8;
//                testImg[96] = 9;
//                testImg[97] = 10;
//                testImg[98] = 11;
//                testImg[99] = 14;
//                testImg[100] = 18;
//                testImg[101] = 16;
//                testImg[102] = 13;
//                testImg[103] = 17;
//                testImg[104] = 19;
//                testImg[105] = 16;
//                testImg[106] = 17;
//                testImg[107] = 11;
//                testImg[108] = 10;
//                testImg[109] = 10;
//                testImg[110] = 10;
//                testImg[111] = 10;
//                testImg[112] = 20;
//                testImg[113] = 22;
//                testImg[114] = 21;
//                testImg[115] = 12;
//                testImg[116] = 13;
//                testImg[117] = 14;
//                testImg[118] = 15;
//                testImg[119] = 12;
//                testImg[120] = 13;
//                testImg[121] = 14;
//                testImg[122] = 15;
//                testImg[123] = 15;
//                testImg[124] = 11;
//                testImg[125] = 12;
//                testImg[126] = 11;
//                testImg[127] = 12;
//                testImg[128] = 18;
//                testImg[129] = 16;
//                testImg[130] = 13;
//                testImg[131] = 17;
//                testImg[132] = 19;
//                testImg[133] = 16;
//                testImg[134] = 17;
//                testImg[135] = 11;
//                testImg[136] = 10;
//                testImg[137] = 10;
//                testImg[138] = 10;
//                testImg[139] = 10;
//                testImg[140] = 20;
//                testImg[141] = 22;
//                testImg[142] = 21;
//                testImg[143] = 12;
//                testImg[144] = 13;
//                testImg[145] = 14;
//                testImg[146] = 15;
//                testImg[147] = 12;
//                testImg[148] = 13;
//                testImg[149] = 14;
//                testImg[150] = 15;
//                testImg[151] = 15;
//                testImg[152] = 11;
//                testImg[153] = 12;
//                testImg[154] = 11;
//                testImg[155] = 12;
//                testImg[156] = 14;
//                testImg[157] = 15;
//                testImg[158] = 15;
//                testImg[159] = 11;
//                testImg[160] = 12;
//                testImg[161] = 11;
//                testImg[162] = 12;
//                testImg[163] = 18;
//                testImg[164] = 16;
//                testImg[165] = 13;
//                testImg[166] = 17;
//                testImg[167] = 19;
//                testImg[168] = 16;
//                testImg[169] = 17;
//                testImg[170] = 11;
//                testImg[171] = 10;
//                testImg[172] = 10;
//                testImg[173] = 10;
//                testImg[174] = 10;

//                testRecImg[0] = 10;
//                testRecImg[1] = 4;
//                testRecImg[2] = 6;
//                testRecImg[3] = 3;
//                testRecImg[4] = 6;
//                testRecImg[5] = 9;
//                testRecImg[6] = 6;
//                testRecImg[7] = 9;
//                testRecImg[8] = 2;
//                testRecImg[9] = 0;
//                testRecImg[10] = 1;
//                testRecImg[11] = 2;
//                testRecImg[12] = 0;
//                testRecImg[13] = 1;
//                testRecImg[14] = 1;
//                testRecImg[15] = 1;
//                testRecImg[16] = 3;
//                testRecImg[17] = 4;
//                testRecImg[18] = 4;
//                testRecImg[19] = 6;
//                testRecImg[20] = 1;
//                testRecImg[21] = 3;
//                testRecImg[22] = 0;
//                testRecImg[23] = 0;
//                testRecImg[24] = 5;
//                testRecImg[25] = 0;
//                testRecImg[26] = 0;
//                testRecImg[27] = 1;
//                testRecImg[28] = 1;
//                testRecImg[29] = 8;
//                testRecImg[30] = 7;
//                testRecImg[31] = 0;
//                testRecImg[32] = 10;
//                testRecImg[33] = 10;
//                testRecImg[34] = 7;
//                testRecImg[35] = 7;
//                testRecImg[36] = 18;
//                testRecImg[37] = 9;
//                testRecImg[38] = 17;
//                testRecImg[39] = 17;
//                testRecImg[40] = 17;
//                testRecImg[41] = 16;
//                testRecImg[42] = 17;
//                testRecImg[43] = 12;
//                testRecImg[44] = 11;
//                testRecImg[45] = 10;
//                testRecImg[46] = 12;
//                testRecImg[47] = 13;
//                testRecImg[48] = 20;
//                testRecImg[49] = 23;
//                testRecImg[50] = 22;
//                testRecImg[51] = 12;
//                testRecImg[52] = 10;
//                testRecImg[53] = 13;
//                testRecImg[54] = 15;
//                testRecImg[55] = 15;
//                testRecImg[56] = 13;
//                testRecImg[57] = 14;
//                testRecImg[58] = 15;
//                testRecImg[59] = 15;
//                testRecImg[60] = 11;
//                testRecImg[61] = 11;
//                testRecImg[62] = 10;
//                testRecImg[63] = 12;
//                testRecImg[64] = 10;
//                testRecImg[65] = 8;
//                testRecImg[66] = 6;
//                testRecImg[67] = 3;
//                testRecImg[68] = 7;
//                testRecImg[69] = 9;
//                testRecImg[70] = 6;
//                testRecImg[71] = 7;
//                testRecImg[72] = 1;
//                testRecImg[73] = 0;
//                testRecImg[74] = 0;
//                testRecImg[75] = 0;
//                testRecImg[76] = 0;
//                testRecImg[77] = 0;
//                testRecImg[78] = 2;
//                testRecImg[79] = 1;
//                testRecImg[80] = 2;
//                testRecImg[81] = 3;
//                testRecImg[82] = 4;
//                testRecImg[83] = 5;
//                testRecImg[84] = 2;
//                testRecImg[85] = 3;
//                testRecImg[86] = 4;
//                testRecImg[87] = 5;
//                testRecImg[88] = 5;
//                testRecImg[89] = 1;
//                testRecImg[90] = 2;
//                testRecImg[91] = 1;
//                testRecImg[92] = 2;
//                testRecImg[93] = 6;
//                testRecImg[94] = 7;
//                testRecImg[95] = 8;
//                testRecImg[96] = 9;
//                testRecImg[97] = 10;
//                testRecImg[98] = 11;
//                testRecImg[99] = 14;
//                testRecImg[100] = 18;
//                testRecImg[101] = 16;
//                testRecImg[102] = 13;
//                testRecImg[103] = 17;
//                testRecImg[104] = 19;
//                testRecImg[105] = 16;
//                testRecImg[106] = 17;
//                testRecImg[107] = 11;
//                testRecImg[108] = 10;
//                testRecImg[109] = 10;
//                testRecImg[110] = 10;
//                testRecImg[111] = 10;
//                testRecImg[112] = 20;
//                testRecImg[113] = 22;
//                testRecImg[114] = 21;
//                testRecImg[115] = 12;
//                testRecImg[116] = 13;
//                testRecImg[117] = 14;
//                testRecImg[118] = 15;
//                testRecImg[119] = 12;
//                testRecImg[120] = 13;
//                testRecImg[121] = 14;
//                testRecImg[122] = 15;
//                testRecImg[123] = 15;
//                testRecImg[124] = 11;
//                testRecImg[125] = 12;
//                testRecImg[126] = 11;
//                testRecImg[127] = 12;
//                testRecImg[128] = 18;
//                testRecImg[129] = 16;
//                testRecImg[130] = 13;
//                testRecImg[131] = 17;
//                testRecImg[132] = 19;
//                testRecImg[133] = 16;
//                testRecImg[134] = 17;
//                testRecImg[135] = 11;
//                testRecImg[136] = 10;
//                testRecImg[137] = 10;
//                testRecImg[138] = 10;
//                testRecImg[139] = 10;
//                testRecImg[140] = 20;
//                testRecImg[141] = 22;
//                testRecImg[142] = 21;
//                testRecImg[143] = 12;
//                testRecImg[144] = 13;
//                testRecImg[145] = 14;
//                testRecImg[146] = 15;
//                testRecImg[147] = 12;
//                testRecImg[148] = 13;
//                testRecImg[149] = 14;
//                testRecImg[150] = 15;
//                testRecImg[151] = 15;
//                testRecImg[152] = 11;
//                testRecImg[153] = 12;
//                testRecImg[154] = 11;
//                testRecImg[155] = 12;
//                testRecImg[156] = 14;
//                testRecImg[157] = 15;
//                testRecImg[158] = 15;
//                testRecImg[159] = 11;
//                testRecImg[160] = 12;
//                testRecImg[161] = 11;
//                testRecImg[162] = 12;
//                testRecImg[163] = 18;
//                testRecImg[164] = 16;
//                testRecImg[165] = 13;
//                testRecImg[166] = 17;
//                testRecImg[167] = 19;
//                testRecImg[168] = 16;
//                testRecImg[169] = 17;
//                testRecImg[170] = 11;
//                testRecImg[171] = 10;
//                testRecImg[172] = 10;
//                testRecImg[173] = 10;
//                testRecImg[174] = 10;
//                int* rpxl1 = new int[sizeof(int) * 8];
//                rpxl1[0] = 0;
//                rpxl1[1] = 6;
//                rpxl1[2] = 17;
//                rpxl1[3] = 33;
//                rpxl1[4] = 90;
//                rpxl1[5] = 115;
//                rpxl1[6] = 141;
//                rpxl1[7] = 168;
//                int dilMsk = 0;
//                double* rpxl2 = dilatedMask3D_partial(7, 5, 5, rpxl1, dilMsk, false, 1);
//                qDebug() << dilMsk;
//                for(int i=0; i<50; i++) {
//                    qDebug() << "Dilat: " << i << rpxl2[i] << testImg[(int)rpxl2[i]];
//                }

        QMessageBox::information(this, "Saved", "The code has run succesfully!");
    }
}
