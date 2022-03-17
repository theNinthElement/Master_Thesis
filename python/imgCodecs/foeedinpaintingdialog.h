#ifndef FOEEDINPAINTINGDIALOG_H
#define FOEEDINPAINTINGDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class FOEEDInpaintingDialog;
}

class FOEEDInpaintingDialog : public QDialog
{
    Q_OBJECT

public:
    explicit FOEEDInpaintingDialog(QWidget *parent = 0);
    ~FOEEDInpaintingDialog();

private Q_SLOTS:
    void on_radioButton_foeedInpt_ExplicitScheme_clicked();

    void on_radioButton_foeedInpt_FEDScheme_clicked();

    void on_radioButton_foeedInpt_FSISceheme_clicked();

    void on_pushButton_foeedInpt_origDataLoad_clicked();

    void on_pushButton_foeedInpt_loadMaskData_clicked();

    void on_pushButton_foeedInpt_run_clicked();

    void on_pushButton_foeedInpt_saveDatName_clicked();

private:
    Ui::FOEEDInpaintingDialog *ui;

    QPixmap image;
    QImage  *imageObject;
    QImage  *maskImageObject;
    QStringList randPxlStrArr;
    nifti_image *nim_input, *nim_input_mask;
};

#endif // FOEEDINPAINTINGDIALOG_H
