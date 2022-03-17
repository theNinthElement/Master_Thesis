#ifndef XORRESDIALOG_H
#define XORRESDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class XORResDialog;
}

class XORResDialog : public QDialog
{
    Q_OBJECT

public:
    explicit XORResDialog(QWidget *parent = 0);
    ~XORResDialog();

private Q_SLOTS:
    void on_pushButton_xorRes_origImgLoad_clicked();

    void on_pushButton_xorRes_recImgLoad_clicked();

    void on_pushButton_xorRes_run_clicked();

    void on_pushButton_xorres_saveImg_clicked();

private:
    Ui::XORResDialog *ui;
    QPixmap imageOrig,imageRec;
    QImage  *origImageObject, *recImageObject;
    nifti_image *nimOrig, *nimRec;
};

#endif // XORRESDIALOG_H
