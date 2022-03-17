#include "xorresdialog.h"
#include "ui_xorresdialog.h"
#include "supplementary_functions.h"
#include "nifti1_io.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <vector>
#include <algorithm>

#include <QDebug>

XORResDialog::XORResDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::XORResDialog)
{
    ui->setupUi(this);

    ui->lineEdit_xorRes_ImgName->setPlaceholderText("XORresidual");
}

XORResDialog::~XORResDialog()
{
    delete ui;
}

void XORResDialog::on_pushButton_xorRes_origImgLoad_clicked()
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
    ui->label_xorRes_OrigImg->setPixmap(imageOrig.scaled(ui->label_xorRes_OrigImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_xorRes_OrigImg->setPixmap(imageMask);
}

void XORResDialog::on_pushButton_xorRes_recImgLoad_clicked()
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
    ui->label_xorRes_recImg->setPixmap(imageRec.scaled(ui->label_xorRes_recImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
//    ui->label_xorRes_recImg->setPixmap(imageMask);
}

void XORResDialog::on_pushButton_xorRes_run_clicked()
{
    if(ui->label_xorRes_OrigImg->pixmap()==0 //.isNull()
            || ui->label_xorRes_recImg->pixmap()==0) //.isNull())
    {
        QMessageBox::warning(this, "Missing Input", "Please enter load all image data");
    } else {
        int imgWidth = nimOrig->dim[1], imgHeight = nimOrig->dim[2], imgDepth = nimOrig->dim[3];
        double *origImageArr, *recImageArr;

        //Memory allocation for residual array.
        short int* resArr = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
        short int* resArr1 = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];

        origImageArr = nii_to_array(nimOrig);
        recImageArr = nii_to_array(nimRec);

        //Post-processing: Clipping **********************************
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            if(recImageArr[i] < 0) {
                recImageArr[i] = 0;
            }
        }
        //End of post-processing *************************************

        double sum = 0;
        //XOR opearation
        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
            resArr[i] = (short int)origImageArr[i] ^ (short int)recImageArr[i];
 //            sum = sum + resArr[i]*resArr[i];
        }
 //        qDebug() << sum/(imgWidth*imgHeight*imgDepth);

        si_array_to_nii(nimRec,resArr);


// // //Half XOR residual *************************************************************************************************************************************
//        //Post-processing: Clipping **********************************
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(recImageArr[i] < 0) {
//                recImageArr[i] = 0;
//            }
//        }
//        //End of post-processing *************************************

//        double sum = 0;
//        //XOR opearation
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            int sign = 1;
//            resArr[i] = (short int)origImageArr[i] ^ (short int)recImageArr[i];
//            if(resArr[i]%2 != 0) {
//                sign = -1;
//                resArr[i] = resArr[i]-1;
//                resArr[i] = resArr[i]/2;
//            } else {
//                sign = 1;
//                resArr[i] = resArr[i]/2;
//            }
//            resArr[i] = sign*resArr[i];
//            sum = sum + resArr[i]*resArr[i];
//        }
//        qDebug() << "Half XOR: " << sum/(imgWidth*imgHeight*imgDepth);

//        si_array_to_nii(nimRec,resArr);
// // //End of Half XOR residual *************************************************************************************************************************************


// // //Partion into 3 array k times *************************************************************************************************************************************
//        short int* resArr2 = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
//        // Array to store all prime less than and equal to 10^6
//        vector <short int> primes;
//        SieveOfEratosthenes(primes,4096);
//        //        for (vector<short int>::iterator it = primes.begin()+1 ; it != primes.end(); ++it) {
//        //            qDebug() << ' ' << *it;
//        //            qDebug() << '\n';
//        //        }
//        int sign = 1;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//       //             if(recImageArr[i] < 0) {
//       //                 sign = -1;
//       //             } else {
//       //                 sign = 1;
//       //             }
//            short int tempVal = (short int)recImageArr[i];
//            short int* a;
//            QString bin_divTracer;

//            tempVal = tempVal*sign;

// //1 division *********************************************************************************************************
//            if(tempVal<=7) {
//                resArr2[i] = 0;
//                resArr1[i] = sign*tempVal;
//                resArr[i] = 0;
//                continue;
//            }

// //            qDebug() << tempVal;

//            // Chech if number is even
//            if (tempVal%2 != 0) {
//                bin_divTracer = bin_divTracer + "1";
//                tempVal = tempVal-1;
//            } else {
//                bin_divTracer = bin_divTracer + "0";
//            }
//            tempVal = tempVal/2;
//            if (tempVal%2 != 0) {
//                bin_divTracer = bin_divTracer + "1";
//                tempVal = tempVal-1;
//            } else {
//                bin_divTracer = bin_divTracer + "0";
//            }
//            bool ok;
//            int divTracer = bin_divTracer.toInt(&ok,2);
// //            qDebug() << tempVal << bin_divTracer << divTracer;
// //End of 1 division *********************************************************************************************************

// //2 division ****************************************************************************************************************
//            if(tempVal<=15) {
//                resArr2[i] = 0;
//                resArr1[i] = sign*tempVal;
//                resArr[i] = 0;
//                continue;
//            }

// //            qDebug() << tempVal;

//            // Chech if number is even or less than 3
//            if (tempVal%2 != 0) {
//                bin_divTracer = bin_divTracer + "1";
//                tempVal = tempVal-1;
//            } else {
//                bin_divTracer = bin_divTracer + "0";
//            }
//            tempVal = tempVal/2;
//            if (tempVal%2 != 0) {
//                bin_divTracer = bin_divTracer + "1";
//                tempVal = tempVal-1;
//            } else {
//                bin_divTracer = bin_divTracer + "0";
//            }
//            tempVal = tempVal/2;
//            if (tempVal%2 != 0) {
//                bin_divTracer = bin_divTracer + "1";
//                tempVal = tempVal-1;
//            } else {
//                bin_divTracer = bin_divTracer + "0";
//            }
//            bool ok;
//            int divTracer = bin_divTracer.toInt(&ok,2);
// //            qDebug() << tempVal << bin_divTracer << divTracer;
// //End 2 division *********************************************************************************************************

// //3 division ****************************************************************************************************************
//        if(tempVal<=31) {
//            resArr2[i] = tempVal % 8;
//            resArr1[i] = sign*((int)tempVal/8);
//            resArr[i] = 0;

// //            qDebug() << tempVal << resArr[i] << resArr1[i] << resArr2[i];
//            continue;
//        }
// //        qDebug() << tempVal;

//        // Chech if number is even
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        bool ok;
//        int divTracer = bin_divTracer.toInt(&ok,2);
// //        qDebug() << tempVal << bin_divTracer << divTracer;
// //End 3 division *********************************************************************************************************

// //4 division ****************************************************************************************************************
//        if(tempVal<=63) {
//            resArr2[i] = tempVal % 8;
//            resArr1[i] = sign*((int)tempVal/8);
//            resArr[i] = 0;

//            //            qDebug() << tempVal << resArr[i] << resArr1[i] << resArr2[i];
//            continue;
//        }
//            //        qDebug() << tempVal;

//        // Chech if number is even
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        tempVal = tempVal/2;
//        if (tempVal%2 != 0) {
//            bin_divTracer = bin_divTracer + "1";
//            tempVal = tempVal-1;
//        } else {
//            bin_divTracer = bin_divTracer + "0";
//        }
//        bool ok;
//        int divTracer = bin_divTracer.toInt(&ok,2);
//        //        qDebug() << tempVal << bin_divTracer << divTracer;
// //End 4 division *********************************************************************************************************

//            a = findPrimes(primes, tempVal);
// //        //            qDebug() << a[0] << a[1];
//            resArr2[i] = divTracer;
//            resArr1[i] = sign*countPrimes(a[0]);
//            resArr[i] = countPrimes(a[1]);
//        }

//        QDir dir("./img-outputs/");
//        if (!dir.exists())
//            dir.mkpath(".");
//        QFile file0("./img-outputs/res0");
//        if(!file0.open(QFile::WriteOnly | QFile::Text)) {
//            QMessageBox::warning(this, "Error", "File is NOT open");
//        }
//        QByteArray temp0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
//            int ind = resArr[i];
//            char buf[7];
//            ::sprintf(buf, "%d", ind);
//            temp0.append(buf);
//            temp0.append("\n");
//        }
//        file0.write(temp0);
//        file0.flush();
//        file0.close();

//        QFile file1("./img-outputs/res1");
//        if(!file1.open(QFile::WriteOnly | QFile::Text)) {
//            QMessageBox::warning(this, "Error", "File is NOT open");
//        }
//        QByteArray temp1;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
//            int ind = resArr1[i];
//            char buf[7];
//            ::sprintf(buf, "%d", ind);
//            temp1.append(buf);
//            temp1.append("\n");
//        }
//        file1.write(temp1);
//        file1.flush();
//        file1.close();

//        QFile file2("./img-outputs/res2");
//        if(!file2.open(QFile::WriteOnly | QFile::Text)) {
//            QMessageBox::warning(this, "Error", "File is NOT open");
//        }
//        QByteArray temp2;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
//            int ind = resArr2[i];
//            char buf[7];
//            ::sprintf(buf, "%d", ind);
//            temp2.append(buf);
//            temp2.append("\n");
//        }
//        file2.write(temp2);
//        file2.flush();
//        file2.close();

//        si_array_to_nii(nimRec,resArr2);
// //End Partion into 3 array k times *************************************************************************************************************************************************************************************


// // //Partion into 3 array *************************************************************************************************************************************
//        short int* resArr2 = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
//        // Array to store all prime less than and equal to 10^6
//        vector <short int> primes;
//        SieveOfEratosthenes(primes,4096);
// //        for (vector<short int>::iterator it = primes.begin()+1 ; it != primes.end(); ++it) {
// //            qDebug() << ' ' << *it;
// //            qDebug() << '\n';
// //        }
//        int sign = 1;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
// //             if(recImageArr[i] < 0) {
// //                 sign = -1;
// //             } else {
// //                 sign = 1;
// //             }
//            short int tempVal = (short int)recImageArr[i];
//            short int* a;

//            tempVal = tempVal*sign;

//            if(tempVal<=1) {
//                resArr2[i] = sign*tempVal;
//                resArr1[i] = 0;
//                resArr[i] = 0;
//                continue;
//            } else if(tempVal == 2) {
//                resArr2[i] = 0;
//                resArr1[i] = sign*tempVal;
//                resArr[i] = 0;
//                continue;
//            } else if(tempVal == 3) {
//                resArr2[i] = 0;
//                resArr1[i] = sign*tempVal;
//                resArr[i] = 0;
//                continue;
//            }

//            // Chech if number is even or less than 3
//            if (tempVal%2 != 0) {
//                resArr2[i] = 1;
//                tempVal = tempVal-1;
//            } else {
//                resArr2[i] = 0;
//            }

// //            qDebug() << tempVal;

//            a = findPrimes(primes, tempVal);
// //            qDebug() << a[0] << a[1];
//            resArr1[i] = sign*countPrimes(a[0]);
//            resArr[i] = countPrimes(a[1]);
//        }

//        QDir dir("./img-outputs/");
//        if (!dir.exists())
//            dir.mkpath(".");
//        QFile file0("./img-outputs/res0");
//        if(!file0.open(QFile::WriteOnly | QFile::Text)) {
//            QMessageBox::warning(this, "Error", "File is NOT open");
//        }
//        QByteArray temp0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
//            int ind = resArr[i];
//            char buf[7];
//            ::sprintf(buf, "%d", ind);
//            temp0.append(buf);
//            temp0.append("\n");
//        }
//        file0.write(temp0);
//        file0.flush();
//        file0.close();

//        QFile file1("./img-outputs/res1");
//        if(!file1.open(QFile::WriteOnly | QFile::Text)) {
//            QMessageBox::warning(this, "Error", "File is NOT open");
//        }
//        QByteArray temp1;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
//            int ind = resArr1[i];
//            char buf[7];
//            ::sprintf(buf, "%d", ind);
//            temp1.append(buf);
//            temp1.append("\n");
//        }
//        file1.write(temp1);
//        file1.flush();
//        file1.close();

//        QFile file2("./img-outputs/res2");
//        if(!file2.open(QFile::WriteOnly | QFile::Text)) {
//            QMessageBox::warning(this, "Error", "File is NOT open");
//        }
//        QByteArray temp2;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {                                 // Randomly choosen X percent of pixels are used to create scattered image
//            int ind = resArr2[i];
//            char buf[7];
//            ::sprintf(buf, "%d", ind);
//            temp2.append(buf);
//            temp2.append("\n");
//        }
//        file2.write(temp2);
//        file2.flush();
//        file2.close();

//        si_array_to_nii(nimRec,resArr2);
// //  *************************************************************************************************************************************************************************************

// //Prime Map conversion of the combination of residual and mask *************************************************************************************************************************************
//        int maxElement = 0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(maxElement < abs(recImageArr[i])) {
//                maxElement = abs(recImageArr[i]);
//            }
//        }
//        vector<short int> primeArr;
//        SieveOfEratosthenes(primeArr,maxElement);

//        //'Number of prime less than given intiger' map for prime voxel values of combined residual
//        int count = 0, sign = 1;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(recImageArr[i] < 0) {
//                sign = -1;
//            } else {
//                sign = 1;
//            }
//            short int tempVal = (short int)recImageArr[i];
//            if(std::find(primeArr.begin(), primeArr.end(), abs(tempVal)) != primeArr.end()) {
//                /* v contains x */
//                resArr[i] = countPrimes(abs(tempVal))*sign;
//                count++;
//            } else {
//                /* v does not contain x */
//                resArr[i] = tempVal;
//            }
//        }
//        qDebug() << "Number of primes: " << count;
//        qDebug() << recImageArr[956833] << resArr[956833];

//        si_array_to_nii(nimRec,resArr);
// //*************************************************************************************************************************************************************************************************



//        //2-power Residual Part*************************************************************************************************************************************
//        //Post-processing: Clipping **********************************
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(recImageArr[i] < 0) {
//                recImageArr[i] = 0;
//            }
//        }
//        //End of post-processing *************************************
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            short int temValR = (short int)recImageArr[i];
//            short int temValO = (short int)origImageArr[i];
//            short int powTemValR = (int)log2(temValR);
//            short int powTemValO = (int)log2(temValO);
//            int zeroR = 1, zeroO = 1;

//            if(temValR == 0) {
//                zeroR = 0;
//            }
//            if(temValO == 0) {
//                zeroO = 0;
//            }
// //            resArr[i] = powTemValO - powTemValR;
//            resArr[i] = powTemValO ^ powTemValR;
// //            resArr1[i] = origImageArr[i] - zeroO*pow(2,powTemValO) - recImageArr[i] + zeroR*pow(2,powTemValR);
//            resArr1[i] = (short int)(origImageArr[i] - zeroO*pow(2,powTemValO)) ^ (short int)(recImageArr[i] - zeroR*pow(2,powTemValR));
//        }

//        qDebug() << recImageArr[101] << resArr[101] << origImageArr[101] << resArr1[101];

//        si_array_to_nii(nimRec,resArr);

//        // Writing scattered nii file
//        si_array_to_nii(nimOrig,resArr1);
//        QString nameStr = "./img-outputs/residuals/xor/";
//        nameStr = nameStr + ui->lineEdit_xorRes_ImgName->text() +  "_" + "XOR_residual_residual" + ".nii";
//        //***************************************************
//        //Convert QString to char*
//        QByteArray baTemp = nameStr.toLocal8Bit();
//        const char *fout = baTemp.data();
//        //***************************************************
//        nifti_set_filenames(nimOrig,fout,1,1);
//        nifti_image_write(nimOrig);
//        //***************************************************
//        //*********************************************************************************************************************************************************


//        //Prime Map Residual Part*************************************************************************************************************************************
//        //'Number of prime less than given intiger' map
//        int maxElement = 0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            short int temVal = (short int)recImageArr[i];
//            resArr[i] = countPrimes(temVal);
// //            resArr[i] = (short int)temVal;
//            if(temVal > maxElement) {                                                  // Trace maximum for getting array of primes less than or equal to maximum
//                maxElement = temVal;
//            }
//        }
//        qDebug() << "Maximum: " << maxElement;

//        vector<short int> primeArr;
//        SieveOfEratosthenes(primeArr,maxElement);
// //        qDebug() << primeArr[0];
//        qDebug() << "Vector size: " << primeArr.size();
// //        for(int i=0; i<primeArr.size(); i++) {
// //            qDebug() << primeArr[i];
// //        }
//        //Memory allocation for residual array.
//        short int* nthPrime = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {

// //            //This part fonly for primeMap smoothing test*******
// //            if(origImageArr[i] < 0) {
// //                origImageArr[i] = 0;
// //            }
// //            short int tem = origImageArr[i];
// //            nthPrime[i] = primeArr[tem];
// //            qDebug() << tem << nthPrime[i];
//            //**************************************************

//             short int tem = resArr[i];
//             nthPrime[i] = primeArr[tem];
// //             qDebug() << "Nth Prime: " << nthPrime[i] << primeArr[tem] << tem << recImageArr[i];
//        }

//        qDebug() << countPrimes(19);
//        int a = resArr[101];
//        qDebug() << recImageArr[101] << resArr[101] << nthPrime[101] << primeArr[a];
//        qDebug() << "Nth Prime: " << nthPrime[101] << primeArr[a] << a << recImageArr[101];
//        //XOR opearation
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
// //             resArr1[i] = nthPrime[i] ^ (short int)recImageArr[i];
//              resArr1[i] =  (short int)recImageArr[i] - nthPrime[i];

// //                        resArr1[i] =  nthPrime[i];

// //            qDebug() << resArr1[i];
//        }

//        si_array_to_nii(nimRec,resArr);

//        // Writing scattered nii file
//        si_array_to_nii(nimOrig,resArr1);
//        QString nameStr = "./img-outputs/residuals/xor/";
//        nameStr = nameStr + ui->lineEdit_xorRes_ImgName->text() +  "_" + "XOR_residual_residual" + ".nii";
//        //***************************************************
//        //Convert QString to char*
//        QByteArray baTemp = nameStr.toLocal8Bit();
//        const char *fout = baTemp.data();
//        //***************************************************
//        nifti_set_filenames(nimOrig,fout,1,1);
//        nifti_image_write(nimOrig);
//        //***************************************************
//        //*********************************************************************************************************************************************************


//        //Near Lossless Compression Part*************************************************************************************************************************************
//        //'Number of prime less than given intiger' map
//        int maxElement = 0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            short int temVal = (short int)recImageArr[i];
//            resArr[i] = countPrimes(temVal);
//            if(temVal > maxElement) {                                                  // Trace maximum for getting array of primes less than or equal to maximum
//                maxElement = temVal;
//            }
//        }

//        qDebug() << "Maximum: " << maxElement;

//        vector<short int> primeArr;
//        SieveOfEratosthenes(primeArr,maxElement);
// //        qDebug() << primeArr[0];
//        qDebug() << "Vector size: " << primeArr.size();
// //        for(int i=0; i<primeArr.size(); i++) {
// //            qDebug() << primeArr[i];
// //        }

//        //Post-processing: Clipping **********************************
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(origImageArr[i] < 0) {
//                origImageArr[i] = 0;
//            }
//        }
//        //Memory allocation for residual array.
//        short int* nthPrime = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            short int tem = resArr[i];
//            nthPrime[i] = primeArr[tem];
// //            qDebug() << "Nth Prime: " << nthPrime[i] << primeArr[tem] << tem << recImageArr[i];
//        }

//        qDebug() << countPrimes(19);
//        int a = resArr[101];
//        qDebug() << recImageArr[101] << resArr[101] << nthPrime[101] << primeArr[a];
//        qDebug() << "Nth Prime: " << nthPrime[101] << primeArr[a] << a << recImageArr[101];
// //        double mseEr = 0;
//        //XOR opearation
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            resArr1[i] =  nthPrime[i] ^ (short int)origImageArr[i];
// //            mseEr = mseEr + (origImageArr[i]-resArr1[i])*(origImageArr[i]-resArr1[i]);
//        }
// //        mseEr = mseEr/(imgWidth*imgHeight*imgDepth);
// //        qDebug() << "MSE: " << mseEr;

//        int index = 46+61*imgWidth+61*imgWidth*imgHeight;
//        qDebug() << origImageArr[index] << recImageArr[index] << resArr[index] << resArr1[index];

//        si_array_to_nii(nimRec,resArr);

//        // Writing scattered nii file
//        si_array_to_nii(nimOrig,resArr1);
//        QString nameStr = "./img-outputs/residuals/xor/";
//        nameStr = nameStr + ui->lineEdit_xorRes_ImgName->text() +  "_" + "near_lossless" + ".nii";
//        //***************************************************
//        //Convert QString to char*
//        QByteArray baTemp = nameStr.toLocal8Bit();
//        const char *fout = baTemp.data();
//        //***************************************************
//        nifti_set_filenames(nimOrig,fout,1,1);
//        nifti_image_write(nimOrig);
//        //***************************************************
//        //*********************************************************************************************************************************************************


//        //Near Lossless Compression Part with Difference*************************************************************************************************************************************
//        //'Number of prime less than given intiger' map
//        int maxElement = 0;
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            short int temVal = (short int)recImageArr[i];
//            if(recImageArr[i] < 0) {
//                temVal = (-1)*temVal;
//                resArr[i] = (-1)*countPrimes(temVal);
//            } else {
//                resArr[i] = countPrimes(temVal);
//            }
//            if(temVal > maxElement) {                                                  // Trace maximum for getting array of primes less than or equal to maximum
//                maxElement = temVal;
//            }
//        }

//        qDebug() << "Maximum: " << maxElement;

//        vector<short int> primeArr;
//        SieveOfEratosthenes(primeArr,maxElement);
// //        qDebug() << primeArr[0];
//        qDebug() << "Vector size: " << primeArr.size();
// //        for(int i=0; i<primeArr.size(); i++) {
// //            qDebug() << primeArr[i];
// //        }

//        //Post-processing: Clipping **********************************
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            if(origImageArr[i] < 0) {
//                origImageArr[i] = 0;
//            }
//        }
//        //Memory allocation for residual array.
//        short int* nthPrime = new short int[sizeof(short int) * imgWidth*imgHeight*imgDepth];
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            short int tem = resArr[i];
//            if(resArr[i] < 0) {
//                tem = (-1)*tem;
//                nthPrime[i] = (-1)*primeArr[tem];
//            } else {
//                nthPrime[i] = primeArr[tem];
//            }
// //            qDebug() << "Nth Prime: " << nthPrime[i] << primeArr[tem] << tem << recImageArr[i];
//        }

// //        qDebug() << countPrimes(19);
// //        int a = resArr[101];
// //        qDebug() << recImageArr[101] << resArr[101] << nthPrime[101] << primeArr[a];
// //        qDebug() << "Nth Prime: " << nthPrime[101] << primeArr[a] << a << recImageArr[101];
// //        double mseEr = 0;
//        //Difference opearation
//        for(int i=0; i<imgWidth*imgHeight*imgDepth; i++) {
//            resArr1[i] =   (short int)origImageArr[i] + nthPrime[i];
// //            mseEr = mseEr + (origImageArr[i]-resArr1[i])*(origImageArr[i]-resArr1[i]);

// //             resArr1[i] =   (short int)recImageArr[i] - nthPrime[i];
//        }
// //        mseEr = mseEr/(imgWidth*imgHeight*imgDepth);
// //        qDebug() << "MSE: " << mseEr;

//        int index = 46+61*imgWidth+61*imgWidth*imgHeight;
//        qDebug() << origImageArr[index] << recImageArr[index] << resArr[index] << resArr1[index];

//        si_array_to_nii(nimRec,resArr);

//        // Writing scattered nii file
//        si_array_to_nii(nimOrig,resArr1);
//        QString nameStr = "./img-outputs/residuals/xor/";
//        nameStr = nameStr + ui->lineEdit_xorRes_ImgName->text() +  "_" + "near_lossless" + ".nii";
//        //***************************************************
//        //Convert QString to char*
//        QByteArray baTemp = nameStr.toLocal8Bit();
//        const char *fout = baTemp.data();
//        //***************************************************
//        nifti_set_filenames(nimOrig,fout,1,1);
//        nifti_image_write(nimOrig);
//        //***************************************************
//        //*********************************************************************************************************************************************************



        si_scaleToUchar(resArr, imgWidth*imgHeight*imgDepth);
        si_array_to_png(resArr, recImageObject, imgWidth, imgHeight);

        imageRec = QPixmap::fromImage(*recImageObject);
        ui->label_xorRes_recImg->setPixmap(imageRec.scaled(ui->label_xorRes_recImg->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size);
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

void XORResDialog::on_pushButton_xorres_saveImg_clicked()
{
    bool isSaved;

    if(ui->label_xorRes_recImg->pixmap()==0 //.isNull()
            || ui->lineEdit_xorRes_ImgName->text().isEmpty()) {
        QMessageBox::warning(this, "Missing Input", "No image to save or/and file name is not given");
    } else {
        QDir dir("./img-outputs/residuals/xor");
        if (!dir.exists())
            dir.mkpath(".");

        QString nameStr = "./img-outputs/residuals/xor/";
        nameStr = nameStr + ui->lineEdit_xorRes_ImgName->text() +  "_" + "XOR_residual" + ".png";
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
