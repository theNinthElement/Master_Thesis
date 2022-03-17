/********************************************************************************
** Form generated from reading UI file 'foeedinpaintingdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_FOEEDINPAINTINGDIALOG_H
#define UI_FOEEDINPAINTINGDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QDialog>
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

class Ui_FOEEDInpaintingDialog
{
public:
    QLabel *label_foeedInpt_imgDisplay;
    QLabel *label_foeedInpt_maskDisplay;
    QGroupBox *groupBox_foeedInpt_Inputs;
    QHBoxLayout *horizontalLayout_7;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_foeedInpt_Tol;
    QLineEdit *lineEdit_foeedInpt_tol;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_foeedInpt_timeStep;
    QLineEdit *lineEdit_foeedInpt_timeStep;
    QRadioButton *radioButton_foeedInpt_ExplicitScheme;
    QRadioButton *radioButton_foeedInpt_FEDScheme;
    QRadioButton *radioButton_foeedInpt_FSISceheme;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_foeedInpt_mseTol;
    QLineEdit *lineEdit_foeedInpt_mseTol;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_foeedInpt_FSIinnerLoopSize;
    QLineEdit *lineEdit_foeedInpt_FSIInnerLoopSize;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_foeedInpt_FEDinnerCycle;
    QLineEdit *lineEdit_foeedInpt_FEDinnerCycle;
    QPushButton *pushButton_foeedInpt_origDataLoad;
    QPushButton *pushButton_foeedInpt_loadMaskData;
    QPushButton *pushButton_foeedInpt_run;
    QCheckBox *checkBox_foeedInpt_monitor;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_foeedInpt_saveDatName;
    QLineEdit *lineEdit_foeedInpt_saveDatName;
    QPushButton *pushButton_foeedInpt_saveDatName;
    QCheckBox *checkBox_foeedIntp_zeros2mask;

    void setupUi(QDialog *FOEEDInpaintingDialog)
    {
        if (FOEEDInpaintingDialog->objectName().isEmpty())
            FOEEDInpaintingDialog->setObjectName(QStringLiteral("FOEEDInpaintingDialog"));
        FOEEDInpaintingDialog->resize(1245, 1093);
        label_foeedInpt_imgDisplay = new QLabel(FOEEDInpaintingDialog);
        label_foeedInpt_imgDisplay->setObjectName(QStringLiteral("label_foeedInpt_imgDisplay"));
        label_foeedInpt_imgDisplay->setGeometry(QRect(90, 60, 531, 551));
        label_foeedInpt_maskDisplay = new QLabel(FOEEDInpaintingDialog);
        label_foeedInpt_maskDisplay->setObjectName(QStringLiteral("label_foeedInpt_maskDisplay"));
        label_foeedInpt_maskDisplay->setGeometry(QRect(710, 600, 471, 431));
        groupBox_foeedInpt_Inputs = new QGroupBox(FOEEDInpaintingDialog);
        groupBox_foeedInpt_Inputs->setObjectName(QStringLiteral("groupBox_foeedInpt_Inputs"));
        groupBox_foeedInpt_Inputs->setGeometry(QRect(740, 90, 421, 461));
        horizontalLayout_7 = new QHBoxLayout(groupBox_foeedInpt_Inputs);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_foeedInpt_Tol = new QLabel(groupBox_foeedInpt_Inputs);
        label_foeedInpt_Tol->setObjectName(QStringLiteral("label_foeedInpt_Tol"));

        horizontalLayout_2->addWidget(label_foeedInpt_Tol);

        lineEdit_foeedInpt_tol = new QLineEdit(groupBox_foeedInpt_Inputs);
        lineEdit_foeedInpt_tol->setObjectName(QStringLiteral("lineEdit_foeedInpt_tol"));

        horizontalLayout_2->addWidget(lineEdit_foeedInpt_tol);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_foeedInpt_timeStep = new QLabel(groupBox_foeedInpt_Inputs);
        label_foeedInpt_timeStep->setObjectName(QStringLiteral("label_foeedInpt_timeStep"));

        horizontalLayout_3->addWidget(label_foeedInpt_timeStep);

        lineEdit_foeedInpt_timeStep = new QLineEdit(groupBox_foeedInpt_Inputs);
        lineEdit_foeedInpt_timeStep->setObjectName(QStringLiteral("lineEdit_foeedInpt_timeStep"));

        horizontalLayout_3->addWidget(lineEdit_foeedInpt_timeStep);


        verticalLayout->addLayout(horizontalLayout_3);

        radioButton_foeedInpt_ExplicitScheme = new QRadioButton(groupBox_foeedInpt_Inputs);
        radioButton_foeedInpt_ExplicitScheme->setObjectName(QStringLiteral("radioButton_foeedInpt_ExplicitScheme"));

        verticalLayout->addWidget(radioButton_foeedInpt_ExplicitScheme);

        radioButton_foeedInpt_FEDScheme = new QRadioButton(groupBox_foeedInpt_Inputs);
        radioButton_foeedInpt_FEDScheme->setObjectName(QStringLiteral("radioButton_foeedInpt_FEDScheme"));

        verticalLayout->addWidget(radioButton_foeedInpt_FEDScheme);

        radioButton_foeedInpt_FSISceheme = new QRadioButton(groupBox_foeedInpt_Inputs);
        radioButton_foeedInpt_FSISceheme->setObjectName(QStringLiteral("radioButton_foeedInpt_FSISceheme"));

        verticalLayout->addWidget(radioButton_foeedInpt_FSISceheme);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        label_foeedInpt_mseTol = new QLabel(groupBox_foeedInpt_Inputs);
        label_foeedInpt_mseTol->setObjectName(QStringLiteral("label_foeedInpt_mseTol"));

        horizontalLayout_4->addWidget(label_foeedInpt_mseTol);

        lineEdit_foeedInpt_mseTol = new QLineEdit(groupBox_foeedInpt_Inputs);
        lineEdit_foeedInpt_mseTol->setObjectName(QStringLiteral("lineEdit_foeedInpt_mseTol"));

        horizontalLayout_4->addWidget(lineEdit_foeedInpt_mseTol);


        verticalLayout->addLayout(horizontalLayout_4);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_foeedInpt_FSIinnerLoopSize = new QLabel(groupBox_foeedInpt_Inputs);
        label_foeedInpt_FSIinnerLoopSize->setObjectName(QStringLiteral("label_foeedInpt_FSIinnerLoopSize"));

        horizontalLayout_5->addWidget(label_foeedInpt_FSIinnerLoopSize);

        lineEdit_foeedInpt_FSIInnerLoopSize = new QLineEdit(groupBox_foeedInpt_Inputs);
        lineEdit_foeedInpt_FSIInnerLoopSize->setObjectName(QStringLiteral("lineEdit_foeedInpt_FSIInnerLoopSize"));

        horizontalLayout_5->addWidget(lineEdit_foeedInpt_FSIInnerLoopSize);


        verticalLayout->addLayout(horizontalLayout_5);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        label_foeedInpt_FEDinnerCycle = new QLabel(groupBox_foeedInpt_Inputs);
        label_foeedInpt_FEDinnerCycle->setObjectName(QStringLiteral("label_foeedInpt_FEDinnerCycle"));

        horizontalLayout_6->addWidget(label_foeedInpt_FEDinnerCycle);

        lineEdit_foeedInpt_FEDinnerCycle = new QLineEdit(groupBox_foeedInpt_Inputs);
        lineEdit_foeedInpt_FEDinnerCycle->setObjectName(QStringLiteral("lineEdit_foeedInpt_FEDinnerCycle"));

        horizontalLayout_6->addWidget(lineEdit_foeedInpt_FEDinnerCycle);


        verticalLayout->addLayout(horizontalLayout_6);

        pushButton_foeedInpt_origDataLoad = new QPushButton(groupBox_foeedInpt_Inputs);
        pushButton_foeedInpt_origDataLoad->setObjectName(QStringLiteral("pushButton_foeedInpt_origDataLoad"));

        verticalLayout->addWidget(pushButton_foeedInpt_origDataLoad);

        pushButton_foeedInpt_loadMaskData = new QPushButton(groupBox_foeedInpt_Inputs);
        pushButton_foeedInpt_loadMaskData->setObjectName(QStringLiteral("pushButton_foeedInpt_loadMaskData"));

        verticalLayout->addWidget(pushButton_foeedInpt_loadMaskData);

        pushButton_foeedInpt_run = new QPushButton(groupBox_foeedInpt_Inputs);
        pushButton_foeedInpt_run->setObjectName(QStringLiteral("pushButton_foeedInpt_run"));

        verticalLayout->addWidget(pushButton_foeedInpt_run);

        checkBox_foeedInpt_monitor = new QCheckBox(groupBox_foeedInpt_Inputs);
        checkBox_foeedInpt_monitor->setObjectName(QStringLiteral("checkBox_foeedInpt_monitor"));

        verticalLayout->addWidget(checkBox_foeedInpt_monitor);


        horizontalLayout_7->addLayout(verticalLayout);

        layoutWidget = new QWidget(FOEEDInpaintingDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(140, 890, 441, 31));
        horizontalLayout = new QHBoxLayout(layoutWidget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_foeedInpt_saveDatName = new QLabel(layoutWidget);
        label_foeedInpt_saveDatName->setObjectName(QStringLiteral("label_foeedInpt_saveDatName"));

        horizontalLayout->addWidget(label_foeedInpt_saveDatName);

        lineEdit_foeedInpt_saveDatName = new QLineEdit(layoutWidget);
        lineEdit_foeedInpt_saveDatName->setObjectName(QStringLiteral("lineEdit_foeedInpt_saveDatName"));

        horizontalLayout->addWidget(lineEdit_foeedInpt_saveDatName);

        pushButton_foeedInpt_saveDatName = new QPushButton(layoutWidget);
        pushButton_foeedInpt_saveDatName->setObjectName(QStringLiteral("pushButton_foeedInpt_saveDatName"));

        horizontalLayout->addWidget(pushButton_foeedInpt_saveDatName);

        checkBox_foeedIntp_zeros2mask = new QCheckBox(FOEEDInpaintingDialog);
        checkBox_foeedIntp_zeros2mask->setObjectName(QStringLiteral("checkBox_foeedIntp_zeros2mask"));
        checkBox_foeedIntp_zeros2mask->setGeometry(QRect(760, 560, 121, 22));

        retranslateUi(FOEEDInpaintingDialog);

        QMetaObject::connectSlotsByName(FOEEDInpaintingDialog);
    } // setupUi

    void retranslateUi(QDialog *FOEEDInpaintingDialog)
    {
        FOEEDInpaintingDialog->setWindowTitle(QApplication::translate("FOEEDInpaintingDialog", "Inpainting by FOEED", Q_NULLPTR));
        label_foeedInpt_imgDisplay->setText(QString());
        label_foeedInpt_maskDisplay->setText(QString());
        groupBox_foeedInpt_Inputs->setTitle(QApplication::translate("FOEEDInpaintingDialog", "Inputs", Q_NULLPTR));
        label_foeedInpt_Tol->setText(QApplication::translate("FOEEDInpaintingDialog", "Tolerance", Q_NULLPTR));
        label_foeedInpt_timeStep->setText(QApplication::translate("FOEEDInpaintingDialog", "Time Step Size", Q_NULLPTR));
        radioButton_foeedInpt_ExplicitScheme->setText(QApplication::translate("FOEEDInpaintingDialog", "Explicit Scheme", Q_NULLPTR));
        radioButton_foeedInpt_FEDScheme->setText(QApplication::translate("FOEEDInpaintingDialog", "Fast Explicit Diffusion(FED) Scheme", Q_NULLPTR));
        radioButton_foeedInpt_FSISceheme->setText(QApplication::translate("FOEEDInpaintingDialog", "Fast Semi-Iterative Scheme(FSI)", Q_NULLPTR));
        label_foeedInpt_mseTol->setText(QApplication::translate("FOEEDInpaintingDialog", "MSE Tolerance", Q_NULLPTR));
        label_foeedInpt_FSIinnerLoopSize->setText(QApplication::translate("FOEEDInpaintingDialog", "Inner Loop Size", Q_NULLPTR));
        label_foeedInpt_FEDinnerCycle->setText(QApplication::translate("FOEEDInpaintingDialog", "Inner Cycle Size", Q_NULLPTR));
        pushButton_foeedInpt_origDataLoad->setText(QApplication::translate("FOEEDInpaintingDialog", "Load Original Data", Q_NULLPTR));
        pushButton_foeedInpt_loadMaskData->setText(QApplication::translate("FOEEDInpaintingDialog", "Load Mask Data", Q_NULLPTR));
        pushButton_foeedInpt_run->setText(QApplication::translate("FOEEDInpaintingDialog", "Run", Q_NULLPTR));
        checkBox_foeedInpt_monitor->setText(QApplication::translate("FOEEDInpaintingDialog", "Motoring", Q_NULLPTR));
        label_foeedInpt_saveDatName->setText(QApplication::translate("FOEEDInpaintingDialog", "Name", Q_NULLPTR));
        pushButton_foeedInpt_saveDatName->setText(QApplication::translate("FOEEDInpaintingDialog", "Save", Q_NULLPTR));
        checkBox_foeedIntp_zeros2mask->setText(QApplication::translate("FOEEDInpaintingDialog", "Zeros to Mask", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class FOEEDInpaintingDialog: public Ui_FOEEDInpaintingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_FOEEDINPAINTINGDIALOG_H
