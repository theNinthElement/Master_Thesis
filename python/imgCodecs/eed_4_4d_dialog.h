#ifndef EED_4_4D_DIALOG_H
#define EED_4_4D_DIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class EED_4_4D_Dialog;
}

class EED_4_4D_Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit EED_4_4D_Dialog(QWidget *parent = 0);
    ~EED_4_4D_Dialog();

private Q_SLOTS:
    void on_radioButton_eed44dInpaintExplitScheme_clicked();

    void on_radioButton_eed44dInpaintFED_clicked();

    void on_radioButton_eed44dInpaintFSI_clicked();

    void on_pushButton_eed44dInpaintUplImg_clicked();

    void on_pushButton_eed44dInpaintUplMask_clicked();

    void on_pushButton_eed44dInpaintRun_clicked();

    void on_pushButton_eed44dInpaintSaveImg_clicked();

    void on_pushButton_eed44dReferImgData_clicked();

//    void on_checkBox_brainBinMask_clicked();

    void on_radioButton_eed_brainBinMask_clicked();

    void on_radioButton_eed_brainBinMask_without_clicked();

    void on_pushButton_brainBinMaskLoad_clicked();

private:
    Ui::EED_4_4D_Dialog *ui;

    QPixmap image;
    QImage  *imageObject;
    QImage  *maskImageObject;
    QStringList randPxlStrArr;
    nifti_image *nim_input, *nim_input_mask, *nim_input_refer, *nim_input_brain_mask;
};

#endif // EED_4_4D_DIALOG_H
