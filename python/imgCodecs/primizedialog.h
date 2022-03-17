#ifndef PRIMIZEDIALOG_H
#define PRIMIZEDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class PrimizeDialog;
}

class PrimizeDialog : public QDialog
{
    Q_OBJECT

public:
    explicit PrimizeDialog(QWidget *parent = 0);
    ~PrimizeDialog();

private Q_SLOTS:
    void on_pushButton_primImg_imgUpl_clicked();

    void on_pushButton_primImgSave_clicked();

    void on_pushButton_primImg_run_clicked();

private:
    Ui::PrimizeDialog *ui;
    QPixmap imagePrimImg;
    QImage  *imagePrimImgObject;
    nifti_image *nim_input_origImg;
};

#endif // PRIMIZEDIALOG_H
