#include "imgcodecs.h"
#include "ui_imgcodecs.h"

ImgCodecs::ImgCodecs(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ImgCodecs)
{
    ui->setupUi(this);
}

ImgCodecs::~ImgCodecs()
{
    delete ui;
}

void ImgCodecs::on_pushButton_RandMask_clicked()
{
//**********************************************************
    //Display Random Known Data Selection Window
    randMaskDialog = new RandomMask(this);
    randMaskDialog->show();
//**********************************************************
}

void ImgCodecs::on_pushButton_RegGridMask_clicked()
{
//**********************************************************
    //Display Regular Known Data Selection Window
    regGridMaskDialog = new RegGridMaskDialog(this);
    regGridMaskDialog->show();
//**********************************************************
}

void ImgCodecs::on_pushButton_clicked()
{
//**********************************************************
    //Display Regular Known Data Selection Window
    eedInpaintingDialog = new EEDInpaintingDialog(this);
    eedInpaintingDialog->show();
//**********************************************************
}

void ImgCodecs::on_pushButton_2_clicked()
{
//**********************************************************
    //Display Regular Known Data Selection Window
    foeedInpaintingDialog = new FOEEDInpaintingDialog(this);
    foeedInpaintingDialog->show();
//**********************************************************
}

void ImgCodecs::on_pushButton_xorRes_clicked()
{
//**********************************************************
    //Display XOR type Residual Window
    XORresDialog = new XORResDialog(this);
    XORresDialog->show();
//**********************************************************
}

void ImgCodecs::on_pushButton_diffRes_clicked()
{
//**********************************************************
    //Display Difference type Residual Window
    DiffresDialog = new DiffResDialog(this);
    DiffresDialog->show();
//**********************************************************
}

void ImgCodecs::on_pushButton__EEDSmoothing_clicked()
{
    eedsmoothingDialog = new EEDSmoothingDialog(this);
    eedsmoothingDialog->show();
}

void ImgCodecs::on_pushButton_ZaxisCombine_clicked()
{
    timeaxiscombinedialog = new TimeAxisCombineDialog(this);
    timeaxiscombinedialog->show();
}

void ImgCodecs::on_pushButton_chainCode_clicked()
{
    chain_coding = new Chain_coding(this);
    chain_coding->show();
}

void ImgCodecs::on_pushButton_2dEED_inpt_clicked()
{
    slicewiseEEDdialog = new SlicewiseEEDDialog(this);
    slicewiseEEDdialog->show();
}

void ImgCodecs::on_pushButton_eed_4_4D_clicked()
{
    eed_4d_inpaintingDiaglog = new EED_4_4D_Dialog(this);
    eed_4d_inpaintingDiaglog->show();
}

void ImgCodecs::on_pushButton_4_clicked()
{
    directNeighborDilationDialog = new DirectNeighborDilationDialog(this);
    directNeighborDilationDialog->show();
}

void ImgCodecs::on_pushButton_5_clicked()
{
    primizeDialog = new PrimizeDialog(this);
    primizeDialog->show();
}

void ImgCodecs::on_pushButton_codecPDE_LinHom_clicked()
{
    r_ilh_0_dialog = new R_ILH_0_Dialog(this);
    r_ilh_0_dialog->show();
}
