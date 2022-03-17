#include "chain_coding.h"
#include "ui_chain_coding.h"
#include <QMessageBox>
#include <QPixmap>
#include <QFileDialog>
#include <QImage>
#include <QFile>
#include <cmath>
#include "supplementary_functions.h"

#include <QDebug>
#include <opencv2/opencv.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <QDoubleValidator>

using namespace cv;

Chain_coding::Chain_coding(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Chain_coding)
{
    ui->setupUi(this);

    // The lower threshold lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_chainCode_cannyLowerThres->setValidator(new QDoubleValidator(0, 999999999, 10, this));
    // The upper threshold lineedit will only accept doubles between 0 and 999999999
    ui->lineEdit_chainCode_upperThres->setValidator(new QDoubleValidator(0, 999999999, 10, this));

    //Make the displayed images center of the label
    ui->label_chainCode_data->setAlignment(Qt::AlignCenter);

    //Set default value
    ui->lineEdit_chainCode_cannyLowerThres->setPlaceholderText("0.04");
    ui->lineEdit_chainCode_upperThres->setPlaceholderText("0.1");
}

Chain_coding::~Chain_coding()
{
    delete ui;
}

void Chain_coding::on_pushButton_chainCode_loadData_clicked()
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

    unsigned char* niiImgArr = nii_to_ucharArray(nim_input);
    // Create image and set to the label
    imageObject = new QImage(niiImgArr, imgWidth, imgHeight, QImage::Format_Indexed8);
    if(imgTimeLen != 1) {
        QMessageBox::warning(this, "Image Dimension Problem", "Unfortunately, 4D data does not fit yet!");
        return;
    }
    image = QPixmap::fromImage(*imageObject);
    ui->label_chainCode_data->setPixmap(image.scaled(ui->label_chainCode_data->size(), Qt::KeepAspectRatio, Qt::SmoothTransformation));  //Fit image size to label size
}

void Chain_coding::on_pushButton_chainCode_run_clicked()
{
    qDebug() << "Hello";
    double*** imageArr;
    double* edges;
    int imgWidth = nim_input->dim[1], imgHeight = nim_input->dim[2], imgDepth = nim_input->dim[3];
    float gridSpcX = nim_input->dx, gridSpcY = nim_input->dy, gridSpcZ = nim_input->dz, gridSpcT = nim_input->dt;

    imageArr = nii_to_3d_array(nim_input);
//    for(int i=0; i<imgHeight; i++) {
//        for(int j=0; j<imgWidth; j++) {
//            qDebug() << i << j << imageArr[40][i][j];
//        }
//    }

    Mat sliceImg(imgWidth,imgHeight,CV_64F);
    Mat dst, detected_edges, sliceImg_gray;
    const char* window_name = "Edge Map";
    int lowThreshold = 0;
    const int max_lowThreshold = 100;
    const int ratio = 3;
    const int kernel_size = 3;

    std::memcpy(sliceImg.data, imageArr[40], imgWidth*imgHeight*sizeof(double));
    dst.create(sliceImg.size(), sliceImg.type());
//    cv::cvtColor(sliceImg, sliceImg_gray, cv::COLOR_BGR2GRAY);
//    cv::namedWindow( window_name, cv::WINDOW_AUTOSIZE );
//    cv::createTrackbar( "Min Threshold:", window_name, &lowThreshold, max_lowThreshold, CannyThreshold );

    blur(sliceImg, detected_edges, Size(3,3));
//    cv::Canny(detected_edges, detected_edges, lowThreshold, lowThreshold*ratio, kernel_size);
    Canny(detected_edges, detected_edges, lowThreshold, lowThreshold*ratio, kernel_size);
    dst = Scalar::all(0);
    sliceImg.copyTo(dst, detected_edges);
//    cv::imshow(window_name, dst);

    namedWindow( "Gray image", WINDOW_AUTOSIZE );
    imshow( "Gray image", detected_edges);
    imwrite("detected_edges.jpg", sliceImg);
    waitKey(0);
    qDebug() << "End";
}
