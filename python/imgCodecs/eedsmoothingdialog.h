#ifndef EEDSMOOTHINGDIALOG_H
#define EEDSMOOTHINGDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class EEDSmoothingDialog;
}

class EEDSmoothingDialog : public QDialog
{
    Q_OBJECT

public:
    explicit EEDSmoothingDialog(QWidget *parent = 0);
    ~EEDSmoothingDialog();

private Q_SLOTS:
    void on_radioButton_EEDSmoothing_explScheme_clicked();

    void on_radioButton_EEDSmoothing_FED_clicked();

    void on_pushButton_EEDSmoothing_svaeImgName_clicked();

    void on_pushButton__EEDSmoothing_LoadImgData_clicked();

    void on_pushButton_EEDSmoothing_run_clicked();

private:
    Ui::EEDSmoothingDialog *ui;
    QPixmap image;
    QImage  *imageObject;
    nifti_image *nim_input;
};

#endif // EEDSMOOTHINGDIALOG_H
