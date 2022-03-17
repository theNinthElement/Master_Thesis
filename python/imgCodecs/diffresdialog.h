#ifndef DIFFRESDIALOG_H
#define DIFFRESDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class DiffResDialog;
}

class DiffResDialog : public QDialog
{
    Q_OBJECT

public:
    explicit DiffResDialog(QWidget *parent = 0);
    ~DiffResDialog();

private Q_SLOTS:
    void on_pushButton_diffRes_origImgLoad_clicked();

    void on_pushButton_diffRes_recImgLoad_clicked();

    void on_pushButton_diffRes_run_clicked();

    void on_pushButton_diffRes_saveImg_clicked();

private:
    Ui::DiffResDialog *ui;
    QPixmap imageOrig,imageRec;
    QImage  *origImageObject, *recImageObject;
    nifti_image *nimOrig, *nimRec;
};

#endif // DIFFRESDIALOG_H
