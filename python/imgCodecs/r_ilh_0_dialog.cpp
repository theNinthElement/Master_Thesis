#include "r_ilh_0_dialog.h"
#include "ui_r_ilh_0_dialog.h"

R_ILH_0_Dialog::R_ILH_0_Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::R_ILH_0_Dialog)
{
    ui->setupUi(this);
}

R_ILH_0_Dialog::~R_ILH_0_Dialog()
{
    delete ui;
}
