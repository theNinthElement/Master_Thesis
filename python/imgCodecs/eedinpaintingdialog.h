#ifndef EEDINPAINTINGDIALOG_H
#define EEDINPAINTINGDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class EEDInpaintingDialog;
}

class EEDInpaintingDialog : public QDialog
{
    Q_OBJECT

public:
    explicit EEDInpaintingDialog(QWidget *parent = 0);
    ~EEDInpaintingDialog();

private Q_SLOTS:
    void on_radioButton_eedInpaintExplitScheme_clicked();

    void on_radioButton_eedInpaintFED_clicked();

    void on_radioButton_eedInpaintFSI_clicked();

    void on_pushButton_eedInpaintUplImg_clicked();

    void on_pushButton_eedInpaintUplMask_clicked();

    void on_pushButton_eedInpaintRun_clicked();

    void on_pushButton_eedInpaintSaveImg_clicked();

private:
    Ui::EEDInpaintingDialog *ui;

    QPixmap image;
    QImage  *imageObject;
    QImage  *maskImageObject;
    QStringList randPxlStrArr, eigVals, eigVecs;
    nifti_image *nim_input, *nim_input_mask, *nim_dti;
//    QVector<double> eigVecs;
};

#endif // EEDINPAINTINGDIALOG_H
