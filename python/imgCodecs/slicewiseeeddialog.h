#ifndef SLICEWISEEEDDIALOG_H
#define SLICEWISEEEDDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class SlicewiseEEDDialog;
}

class SlicewiseEEDDialog : public QDialog
{
    Q_OBJECT

public:
    explicit SlicewiseEEDDialog(QWidget *parent = 0);
    ~SlicewiseEEDDialog();

private Q_SLOTS:
    void on_radioButton_2deed_ExScheme_clicked();

    void on_radioButton_FEDScheme_clicked();

    void on_radioButton_FSIScheme_clicked();

    void on_pushButton_loadImg_clicked();

    void on_pushButton_mask_clicked();

    void on_pushButton_run2dEED_clicked();

    void on_pushButton_saveImgName_clicked();

private:
    Ui::SlicewiseEEDDialog *ui;

    QPixmap image;
    QImage  *imageObject;
    QImage  *maskImageObject;
    QStringList randPxlStrArr;
    nifti_image *nim_input, *nim_input_mask;
};

#endif // SLICEWISEEEDDIALOG_H
