#ifndef REGGRIDMASKDIALOG_H
#define REGGRIDMASKDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class RegGridMaskDialog;
}

class RegGridMaskDialog : public QDialog
{
    Q_OBJECT

public:
    explicit RegGridMaskDialog(QWidget *parent = 0);
    ~RegGridMaskDialog();

private Q_SLOTS:
    void on_pushButton_regMask_imgUpl_clicked();

    void on_pushButton_regMask_run_clicked();

    void on_pushButton_regMaskImgSave_clicked();

private:
    Ui::RegGridMaskDialog *ui;
    QPixmap imageRegMask;
    QImage  *imageRegMaskObject;
    nifti_image *nim_input_reg, *nim_grid_mask;
};

#endif // REGGRIDMASKDIALOG_H
