#ifndef IMGCODECS_H
#define IMGCODECS_H

#include <QMainWindow>
#include "randommask.h"
#include "reggridmaskdialog.h"
#include "eedinpaintingdialog.h"
#include "eedsmoothingdialog.h"
#include "foeedinpaintingdialog.h"
#include "xorresdialog.h"
#include "diffresdialog.h"
#include "timeaxiscombinedialog.h"
#include "chain_coding.h"
#include "slicewiseeeddialog.h"
#include "eed_4_4d_dialog.h"
#include "directneighbordilationdialog.h"
#include "primizedialog.h"
#include "r_ilh_0_dialog.h"

namespace Ui {
class ImgCodecs;
}

class ImgCodecs : public QMainWindow
{
    Q_OBJECT

public:
    explicit ImgCodecs(QWidget *parent = 0);
    ~ImgCodecs();

private Q_SLOTS:
    void on_pushButton_RandMask_clicked();

    void on_pushButton_RegGridMask_clicked();

    void on_pushButton_clicked();

    void on_pushButton_xorRes_clicked();

    void on_pushButton_diffRes_clicked();

    void on_pushButton_2_clicked();

    void on_pushButton__EEDSmoothing_clicked();

    void on_pushButton_ZaxisCombine_clicked();

    void on_pushButton_chainCode_clicked();

    void on_pushButton_2dEED_inpt_clicked();

    void on_pushButton_eed_4_4D_clicked();

    void on_pushButton_4_clicked();

    void on_pushButton_5_clicked();

    void on_pushButton_codecPDE_LinHom_clicked();

private:
    Ui::ImgCodecs *ui;
    RandomMask* randMaskDialog;
    RegGridMaskDialog* regGridMaskDialog;
    EEDInpaintingDialog* eedInpaintingDialog;
    FOEEDInpaintingDialog* foeedInpaintingDialog;
    XORResDialog* XORresDialog;
    DiffResDialog* DiffresDialog;
    EEDSmoothingDialog* eedsmoothingDialog;
    TimeAxisCombineDialog* timeaxiscombinedialog;
    Chain_coding* chain_coding;
    SlicewiseEEDDialog* slicewiseEEDdialog;
    EED_4_4D_Dialog* eed_4d_inpaintingDiaglog;
    DirectNeighborDilationDialog* directNeighborDilationDialog;
    PrimizeDialog* primizeDialog;
    R_ILH_0_Dialog* r_ilh_0_dialog;
};

#endif // IMGCODECS_H
