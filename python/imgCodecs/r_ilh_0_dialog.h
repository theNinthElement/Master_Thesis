#ifndef R_ILH_0_DIALOG_H
#define R_ILH_0_DIALOG_H

#include <QDialog>

namespace Ui {
class R_ILH_0_Dialog;
}

class R_ILH_0_Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit R_ILH_0_Dialog(QWidget *parent = 0);
    ~R_ILH_0_Dialog();

private:
    Ui::R_ILH_0_Dialog *ui;
};

#endif // R_ILH_0_DIALOG_H
