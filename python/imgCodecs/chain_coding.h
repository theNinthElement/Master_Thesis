#ifndef CHAIN_CODING_H
#define CHAIN_CODING_H

#include <QDialog>
#include "nifti1_io.h"

namespace Ui {
class Chain_coding;
}

class Chain_coding : public QDialog
{
    Q_OBJECT

public:
    explicit Chain_coding(QWidget *parent = 0);
    ~Chain_coding();

private Q_SLOTS:
    void on_pushButton_chainCode_loadData_clicked();

    void on_pushButton_chainCode_run_clicked();

private:
    Ui::Chain_coding *ui;
    QPixmap image;
    QImage  *imageObject;
    nifti_image *nim_input;
};

#endif // CHAIN_CODING_H
