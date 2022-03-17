#ifndef TIMEAXISCOMBINEDIALOG_H
#define TIMEAXISCOMBINEDIALOG_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class TimeAxisCombineDialog;
}

class TimeAxisCombineDialog : public QDialog
{
    Q_OBJECT

public:
    explicit TimeAxisCombineDialog(QWidget *parent = 0);
    ~TimeAxisCombineDialog();

private Q_SLOTS:
    void on_pushButton_timeComb_UplData1_clicked();

    void on_pushButton_timeComb_UplData2_clicked();

    void on_pushButton_timeComb_saveImg_clicked();

    void on_pushButton_timeComb_combine_clicked();

    void on_pushButton_outData_clicked();

private:
    Ui::TimeAxisCombineDialog *ui;
    QPixmap imageData1,imageData2;
    QImage  *imgDataObject1, *imgDataObject2;
    nifti_image *nimData1, *nimData2, *nimOutData;
};

#endif // TIMEAXISCOMBINEDIALOG_H
