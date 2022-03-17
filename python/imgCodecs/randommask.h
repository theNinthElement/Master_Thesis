#ifndef RANDOMMASK_H
#define RANDOMMASK_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class RandomMask;
}

class RandomMask : public QDialog
{
    Q_OBJECT

public:
    explicit RandomMask(QWidget *parent = 0);
    ~RandomMask();

private Q_SLOTS:
    void on_pushButton_randMask_imgUpl_clicked();

    void on_pushButton_randMask_run_clicked();

    void on_pushButton_randMask_fileNameSave_clicked();

private:
    Ui::RandomMask *ui;
    QPixmap imageMask;
    QImage  *imageMaskObject;
    nifti_image *nim_input;
};

#endif // RANDOMMASK_H
