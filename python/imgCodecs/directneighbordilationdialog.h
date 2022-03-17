#ifndef DIRECTNEIGHBORDILATIONDIALOG_H
#define DIRECTNEIGHBORDILATIONDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class DirectNeighborDilationDialog;
}

class DirectNeighborDilationDialog : public QDialog
{
    Q_OBJECT

public:
    explicit DirectNeighborDilationDialog(QWidget *parent = 0);
    ~DirectNeighborDilationDialog();

private Q_SLOTS:
    void on_pushButton_directDilation_imgUpl_clicked();

    void on_pushButton_directDilation_mskUpl_clicked();

    void on_pushButton_directDilation_recUpl_clicked();

    void on_pushButton_directDilationImgSave_clicked();

    void on_pushButton_directDilation_run_clicked();

private:
    Ui::DirectNeighborDilationDialog *ui;

    QPixmap imageDirDilMask, imageDirDilOrigImg, imageDirDilRecImg;
    QImage  *imageDirDilMaskObject, *imageDirDilOrigImgObject, *imageDirDilRecImgObject;
    QStringList randPxlStrArr;
    nifti_image *nim_input_orig, *nim_grid_mask, *nim_grid_recons;
};

#endif // DIRECTNEIGHBORDILATIONDIALOG_H
