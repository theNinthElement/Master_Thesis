/********************************************************************************
** Form generated from reading UI file 'eed_4_4d_dialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EED_4_4D_DIALOG_H
#define UI_EED_4_4D_DIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_EED_4_4D_Dialog
{
public:
    QLabel *label_eed44dInpaintImg;
    QGroupBox *groupBox_eed44dInpaint;
    QHBoxLayout *horizontalLayout_7;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_eed44dInpaintTol;
    QLineEdit *lineEdit_eed44dInpaintTol;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_eed44dInpaintTimeStep;
    QLineEdit *lineEdit_eed44dInpaintTimeStep;
    QRadioButton *radioButton_eed44dInpaintExplitScheme;
    QRadioButton *radioButton_eed44dInpaintFED;
    QRadioButton *radioButton_eed44dInpaintFSI;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_eed44dInpaintMSETol;
    QLineEdit *lineEdit_eed44dInpaintMSETol;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_eed44dInpaintInnerFEDCycSize;
    QLineEdit *lineEdit_eed44dInpaintInnerFEDCycSize;
    QHBoxLayout *horizontalLayout;
    QLabel *label_eed44dInpaintFSILoopSize;
    QLineEdit *lineEdit_eed44dInpaintFSILoopSoze;
    QPushButton *pushButton_eed44dReferImgData;
    QPushButton *pushButton_eed44dInpaintUplImg;
    QPushButton *pushButton_eed44dInpaintUplMask;
    QPushButton *pushButton_brainBinMaskLoad;
    QFrame *line;
    QPushButton *pushButton_eed44dInpaintRun;
    QCheckBox *checkBox_eed44dInpaintMonitoring;
    QRadioButton *radioButton_eed_brainBinMask;
    QRadioButton *radioButton_eed_brainBinMask_without;
    QLabel *label_eed44dInpaintMaskImg;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_eed44dInpaintFileName;
    QLineEdit *lineEdit_eed44dInpaintFileName;
    QPushButton *pushButton_eed44dInpaintSaveImg;

    void setupUi(QDialog *EED_4_4D_Dialog)
    {
        if (EED_4_4D_Dialog->objectName().isEmpty())
            EED_4_4D_Dialog->setObjectName(QStringLiteral("EED_4_4D_Dialog"));
        EED_4_4D_Dialog->resize(1454, 1115);
        label_eed44dInpaintImg = new QLabel(EED_4_4D_Dialog);
        label_eed44dInpaintImg->setObjectName(QStringLiteral("label_eed44dInpaintImg"));
        label_eed44dInpaintImg->setGeometry(QRect(50, 40, 611, 601));
        groupBox_eed44dInpaint = new QGroupBox(EED_4_4D_Dialog);
        groupBox_eed44dInpaint->setObjectName(QStringLiteral("groupBox_eed44dInpaint"));
        groupBox_eed44dInpaint->setGeometry(QRect(1010, 50, 400, 411));
        horizontalLayout_7 = new QHBoxLayout(groupBox_eed44dInpaint);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        label_eed44dInpaintTol = new QLabel(groupBox_eed44dInpaint);
        label_eed44dInpaintTol->setObjectName(QStringLiteral("label_eed44dInpaintTol"));

        horizontalLayout_6->addWidget(label_eed44dInpaintTol);

        lineEdit_eed44dInpaintTol = new QLineEdit(groupBox_eed44dInpaint);
        lineEdit_eed44dInpaintTol->setObjectName(QStringLiteral("lineEdit_eed44dInpaintTol"));

        horizontalLayout_6->addWidget(lineEdit_eed44dInpaintTol);


        verticalLayout->addLayout(horizontalLayout_6);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_eed44dInpaintTimeStep = new QLabel(groupBox_eed44dInpaint);
        label_eed44dInpaintTimeStep->setObjectName(QStringLiteral("label_eed44dInpaintTimeStep"));

        horizontalLayout_2->addWidget(label_eed44dInpaintTimeStep);

        lineEdit_eed44dInpaintTimeStep = new QLineEdit(groupBox_eed44dInpaint);
        lineEdit_eed44dInpaintTimeStep->setObjectName(QStringLiteral("lineEdit_eed44dInpaintTimeStep"));

        horizontalLayout_2->addWidget(lineEdit_eed44dInpaintTimeStep);


        verticalLayout->addLayout(horizontalLayout_2);

        radioButton_eed44dInpaintExplitScheme = new QRadioButton(groupBox_eed44dInpaint);
        radioButton_eed44dInpaintExplitScheme->setObjectName(QStringLiteral("radioButton_eed44dInpaintExplitScheme"));

        verticalLayout->addWidget(radioButton_eed44dInpaintExplitScheme);

        radioButton_eed44dInpaintFED = new QRadioButton(groupBox_eed44dInpaint);
        radioButton_eed44dInpaintFED->setObjectName(QStringLiteral("radioButton_eed44dInpaintFED"));

        verticalLayout->addWidget(radioButton_eed44dInpaintFED);

        radioButton_eed44dInpaintFSI = new QRadioButton(groupBox_eed44dInpaint);
        radioButton_eed44dInpaintFSI->setObjectName(QStringLiteral("radioButton_eed44dInpaintFSI"));

        verticalLayout->addWidget(radioButton_eed44dInpaintFSI);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_eed44dInpaintMSETol = new QLabel(groupBox_eed44dInpaint);
        label_eed44dInpaintMSETol->setObjectName(QStringLiteral("label_eed44dInpaintMSETol"));

        horizontalLayout_3->addWidget(label_eed44dInpaintMSETol);

        lineEdit_eed44dInpaintMSETol = new QLineEdit(groupBox_eed44dInpaint);
        lineEdit_eed44dInpaintMSETol->setObjectName(QStringLiteral("lineEdit_eed44dInpaintMSETol"));

        horizontalLayout_3->addWidget(lineEdit_eed44dInpaintMSETol);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_eed44dInpaintInnerFEDCycSize = new QLabel(groupBox_eed44dInpaint);
        label_eed44dInpaintInnerFEDCycSize->setObjectName(QStringLiteral("label_eed44dInpaintInnerFEDCycSize"));

        horizontalLayout_5->addWidget(label_eed44dInpaintInnerFEDCycSize);

        lineEdit_eed44dInpaintInnerFEDCycSize = new QLineEdit(groupBox_eed44dInpaint);
        lineEdit_eed44dInpaintInnerFEDCycSize->setObjectName(QStringLiteral("lineEdit_eed44dInpaintInnerFEDCycSize"));

        horizontalLayout_5->addWidget(lineEdit_eed44dInpaintInnerFEDCycSize);


        verticalLayout->addLayout(horizontalLayout_5);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_eed44dInpaintFSILoopSize = new QLabel(groupBox_eed44dInpaint);
        label_eed44dInpaintFSILoopSize->setObjectName(QStringLiteral("label_eed44dInpaintFSILoopSize"));

        horizontalLayout->addWidget(label_eed44dInpaintFSILoopSize);

        lineEdit_eed44dInpaintFSILoopSoze = new QLineEdit(groupBox_eed44dInpaint);
        lineEdit_eed44dInpaintFSILoopSoze->setObjectName(QStringLiteral("lineEdit_eed44dInpaintFSILoopSoze"));

        horizontalLayout->addWidget(lineEdit_eed44dInpaintFSILoopSoze);


        verticalLayout->addLayout(horizontalLayout);

        pushButton_eed44dReferImgData = new QPushButton(groupBox_eed44dInpaint);
        pushButton_eed44dReferImgData->setObjectName(QStringLiteral("pushButton_eed44dReferImgData"));

        verticalLayout->addWidget(pushButton_eed44dReferImgData);

        pushButton_eed44dInpaintUplImg = new QPushButton(groupBox_eed44dInpaint);
        pushButton_eed44dInpaintUplImg->setObjectName(QStringLiteral("pushButton_eed44dInpaintUplImg"));

        verticalLayout->addWidget(pushButton_eed44dInpaintUplImg);

        pushButton_eed44dInpaintUplMask = new QPushButton(groupBox_eed44dInpaint);
        pushButton_eed44dInpaintUplMask->setObjectName(QStringLiteral("pushButton_eed44dInpaintUplMask"));

        verticalLayout->addWidget(pushButton_eed44dInpaintUplMask);

        pushButton_brainBinMaskLoad = new QPushButton(groupBox_eed44dInpaint);
        pushButton_brainBinMaskLoad->setObjectName(QStringLiteral("pushButton_brainBinMaskLoad"));

        verticalLayout->addWidget(pushButton_brainBinMaskLoad);

        line = new QFrame(groupBox_eed44dInpaint);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_eed44dInpaintRun = new QPushButton(groupBox_eed44dInpaint);
        pushButton_eed44dInpaintRun->setObjectName(QStringLiteral("pushButton_eed44dInpaintRun"));

        verticalLayout->addWidget(pushButton_eed44dInpaintRun);

        checkBox_eed44dInpaintMonitoring = new QCheckBox(groupBox_eed44dInpaint);
        checkBox_eed44dInpaintMonitoring->setObjectName(QStringLiteral("checkBox_eed44dInpaintMonitoring"));

        verticalLayout->addWidget(checkBox_eed44dInpaintMonitoring);

        radioButton_eed_brainBinMask = new QRadioButton(groupBox_eed44dInpaint);
        radioButton_eed_brainBinMask->setObjectName(QStringLiteral("radioButton_eed_brainBinMask"));

        verticalLayout->addWidget(radioButton_eed_brainBinMask);

        radioButton_eed_brainBinMask_without = new QRadioButton(groupBox_eed44dInpaint);
        radioButton_eed_brainBinMask_without->setObjectName(QStringLiteral("radioButton_eed_brainBinMask_without"));

        verticalLayout->addWidget(radioButton_eed_brainBinMask_without);


        horizontalLayout_7->addLayout(verticalLayout);

        label_eed44dInpaintMaskImg = new QLabel(EED_4_4D_Dialog);
        label_eed44dInpaintMaskImg->setObjectName(QStringLiteral("label_eed44dInpaintMaskImg"));
        label_eed44dInpaintMaskImg->setGeometry(QRect(960, 570, 450, 450));
        layoutWidget = new QWidget(EED_4_4D_Dialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(110, 960, 471, 31));
        horizontalLayout_4 = new QHBoxLayout(layoutWidget);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        label_eed44dInpaintFileName = new QLabel(layoutWidget);
        label_eed44dInpaintFileName->setObjectName(QStringLiteral("label_eed44dInpaintFileName"));

        horizontalLayout_4->addWidget(label_eed44dInpaintFileName);

        lineEdit_eed44dInpaintFileName = new QLineEdit(layoutWidget);
        lineEdit_eed44dInpaintFileName->setObjectName(QStringLiteral("lineEdit_eed44dInpaintFileName"));

        horizontalLayout_4->addWidget(lineEdit_eed44dInpaintFileName);

        pushButton_eed44dInpaintSaveImg = new QPushButton(layoutWidget);
        pushButton_eed44dInpaintSaveImg->setObjectName(QStringLiteral("pushButton_eed44dInpaintSaveImg"));

        horizontalLayout_4->addWidget(pushButton_eed44dInpaintSaveImg);


        retranslateUi(EED_4_4D_Dialog);

        QMetaObject::connectSlotsByName(EED_4_4D_Dialog);
    } // setupUi

    void retranslateUi(QDialog *EED_4_4D_Dialog)
    {
        EED_4_4D_Dialog->setWindowTitle(QApplication::translate("EED_4_4D_Dialog", "Edge-Enhancing Diffusion based inpainting for 4D data", Q_NULLPTR));
        label_eed44dInpaintImg->setText(QString());
        groupBox_eed44dInpaint->setTitle(QApplication::translate("EED_4_4D_Dialog", "Inputs", Q_NULLPTR));
        label_eed44dInpaintTol->setText(QApplication::translate("EED_4_4D_Dialog", "Tolerance", Q_NULLPTR));
        label_eed44dInpaintTimeStep->setText(QApplication::translate("EED_4_4D_Dialog", "Time Step Size", Q_NULLPTR));
        radioButton_eed44dInpaintExplitScheme->setText(QApplication::translate("EED_4_4D_Dialog", "Explicit Scheme", Q_NULLPTR));
        radioButton_eed44dInpaintFED->setText(QApplication::translate("EED_4_4D_Dialog", "Fast Explicit Diffusion (FED) Scheme", Q_NULLPTR));
        radioButton_eed44dInpaintFSI->setText(QApplication::translate("EED_4_4D_Dialog", "Fast Semi_iterative Scheme (FSI)", Q_NULLPTR));
        label_eed44dInpaintMSETol->setText(QApplication::translate("EED_4_4D_Dialog", "MSE Tolerance", Q_NULLPTR));
        label_eed44dInpaintInnerFEDCycSize->setText(QApplication::translate("EED_4_4D_Dialog", "Inner Cycle Size", Q_NULLPTR));
        label_eed44dInpaintFSILoopSize->setText(QApplication::translate("EED_4_4D_Dialog", "Inner Loop Size", Q_NULLPTR));
        pushButton_eed44dReferImgData->setText(QApplication::translate("EED_4_4D_Dialog", "Load Reference Image Data", Q_NULLPTR));
        pushButton_eed44dInpaintUplImg->setText(QApplication::translate("EED_4_4D_Dialog", "Load Image Data", Q_NULLPTR));
        pushButton_eed44dInpaintUplMask->setText(QApplication::translate("EED_4_4D_Dialog", "Load Mask Data", Q_NULLPTR));
        pushButton_brainBinMaskLoad->setText(QApplication::translate("EED_4_4D_Dialog", "Load Brain Binary Mask", Q_NULLPTR));
        pushButton_eed44dInpaintRun->setText(QApplication::translate("EED_4_4D_Dialog", "Run", Q_NULLPTR));
        checkBox_eed44dInpaintMonitoring->setText(QApplication::translate("EED_4_4D_Dialog", "Monitoring", Q_NULLPTR));
        radioButton_eed_brainBinMask->setText(QApplication::translate("EED_4_4D_Dialog", "with Brain Binary Mask", Q_NULLPTR));
        radioButton_eed_brainBinMask_without->setText(QApplication::translate("EED_4_4D_Dialog", "without Brain Binary Mask", Q_NULLPTR));
        label_eed44dInpaintMaskImg->setText(QString());
        label_eed44dInpaintFileName->setText(QApplication::translate("EED_4_4D_Dialog", "File Name", Q_NULLPTR));
        pushButton_eed44dInpaintSaveImg->setText(QApplication::translate("EED_4_4D_Dialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class EED_4_4D_Dialog: public Ui_EED_4_4D_Dialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EED_4_4D_DIALOG_H
