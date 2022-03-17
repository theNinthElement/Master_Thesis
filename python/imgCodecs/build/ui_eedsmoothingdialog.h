/********************************************************************************
** Form generated from reading UI file 'eedsmoothingdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EEDSMOOTHINGDIALOG_H
#define UI_EEDSMOOTHINGDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
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

class Ui_EEDSmoothingDialog
{
public:
    QLabel *label_EEDSmoothing_ImgDisplay;
    QGroupBox *groupBox_EEDSmoothing_inputs;
    QVBoxLayout *verticalLayout_2;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label__EEDSmoothing_timeStep;
    QLineEdit *lineEdit_EEDSmoothing_timeStep;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label__EEDSmoothing_iterNum;
    QLineEdit *lineEdit_EEDSmoothing_iterNum;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label__EEDSmoothing_contPar;
    QLineEdit *lineEdit_EEDSmoothing_contPar;
    QPushButton *pushButton__EEDSmoothing_LoadImgData;
    QRadioButton *radioButton_EEDSmoothing_explScheme;
    QRadioButton *radioButton_EEDSmoothing_FED;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_EEDSmoothing_FEDInnerCycles;
    QLineEdit *lineEdit_EEDSmoothing_FEDInnerCycles;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_EEDSmoothing_totTime;
    QLineEdit *lineEdit_EEDSmoothing_TotTime;
    QPushButton *pushButton_EEDSmoothing_run;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_EEDSmoothing_saveImName;
    QLineEdit *lineEdit__EEDSmoothing_saveImgName;
    QPushButton *pushButton_EEDSmoothing_svaeImgName;

    void setupUi(QDialog *EEDSmoothingDialog)
    {
        if (EEDSmoothingDialog->objectName().isEmpty())
            EEDSmoothingDialog->setObjectName(QStringLiteral("EEDSmoothingDialog"));
        EEDSmoothingDialog->resize(1146, 1000);
        label_EEDSmoothing_ImgDisplay = new QLabel(EEDSmoothingDialog);
        label_EEDSmoothing_ImgDisplay->setObjectName(QStringLiteral("label_EEDSmoothing_ImgDisplay"));
        label_EEDSmoothing_ImgDisplay->setGeometry(QRect(40, 40, 511, 491));
        groupBox_EEDSmoothing_inputs = new QGroupBox(EEDSmoothingDialog);
        groupBox_EEDSmoothing_inputs->setObjectName(QStringLiteral("groupBox_EEDSmoothing_inputs"));
        groupBox_EEDSmoothing_inputs->setGeometry(QRect(690, 30, 371, 301));
        verticalLayout_2 = new QVBoxLayout(groupBox_EEDSmoothing_inputs);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label__EEDSmoothing_timeStep = new QLabel(groupBox_EEDSmoothing_inputs);
        label__EEDSmoothing_timeStep->setObjectName(QStringLiteral("label__EEDSmoothing_timeStep"));

        horizontalLayout_2->addWidget(label__EEDSmoothing_timeStep);

        lineEdit_EEDSmoothing_timeStep = new QLineEdit(groupBox_EEDSmoothing_inputs);
        lineEdit_EEDSmoothing_timeStep->setObjectName(QStringLiteral("lineEdit_EEDSmoothing_timeStep"));

        horizontalLayout_2->addWidget(lineEdit_EEDSmoothing_timeStep);


        verticalLayout->addLayout(horizontalLayout_2);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label__EEDSmoothing_iterNum = new QLabel(groupBox_EEDSmoothing_inputs);
        label__EEDSmoothing_iterNum->setObjectName(QStringLiteral("label__EEDSmoothing_iterNum"));

        horizontalLayout_3->addWidget(label__EEDSmoothing_iterNum);

        lineEdit_EEDSmoothing_iterNum = new QLineEdit(groupBox_EEDSmoothing_inputs);
        lineEdit_EEDSmoothing_iterNum->setObjectName(QStringLiteral("lineEdit_EEDSmoothing_iterNum"));

        horizontalLayout_3->addWidget(lineEdit_EEDSmoothing_iterNum);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        label__EEDSmoothing_contPar = new QLabel(groupBox_EEDSmoothing_inputs);
        label__EEDSmoothing_contPar->setObjectName(QStringLiteral("label__EEDSmoothing_contPar"));

        horizontalLayout_4->addWidget(label__EEDSmoothing_contPar);

        lineEdit_EEDSmoothing_contPar = new QLineEdit(groupBox_EEDSmoothing_inputs);
        lineEdit_EEDSmoothing_contPar->setObjectName(QStringLiteral("lineEdit_EEDSmoothing_contPar"));

        horizontalLayout_4->addWidget(lineEdit_EEDSmoothing_contPar);


        verticalLayout->addLayout(horizontalLayout_4);

        pushButton__EEDSmoothing_LoadImgData = new QPushButton(groupBox_EEDSmoothing_inputs);
        pushButton__EEDSmoothing_LoadImgData->setObjectName(QStringLiteral("pushButton__EEDSmoothing_LoadImgData"));

        verticalLayout->addWidget(pushButton__EEDSmoothing_LoadImgData);

        radioButton_EEDSmoothing_explScheme = new QRadioButton(groupBox_EEDSmoothing_inputs);
        radioButton_EEDSmoothing_explScheme->setObjectName(QStringLiteral("radioButton_EEDSmoothing_explScheme"));

        verticalLayout->addWidget(radioButton_EEDSmoothing_explScheme);

        radioButton_EEDSmoothing_FED = new QRadioButton(groupBox_EEDSmoothing_inputs);
        radioButton_EEDSmoothing_FED->setObjectName(QStringLiteral("radioButton_EEDSmoothing_FED"));

        verticalLayout->addWidget(radioButton_EEDSmoothing_FED);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_EEDSmoothing_FEDInnerCycles = new QLabel(groupBox_EEDSmoothing_inputs);
        label_EEDSmoothing_FEDInnerCycles->setObjectName(QStringLiteral("label_EEDSmoothing_FEDInnerCycles"));

        horizontalLayout_5->addWidget(label_EEDSmoothing_FEDInnerCycles);

        lineEdit_EEDSmoothing_FEDInnerCycles = new QLineEdit(groupBox_EEDSmoothing_inputs);
        lineEdit_EEDSmoothing_FEDInnerCycles->setObjectName(QStringLiteral("lineEdit_EEDSmoothing_FEDInnerCycles"));

        horizontalLayout_5->addWidget(lineEdit_EEDSmoothing_FEDInnerCycles);


        verticalLayout->addLayout(horizontalLayout_5);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        label_EEDSmoothing_totTime = new QLabel(groupBox_EEDSmoothing_inputs);
        label_EEDSmoothing_totTime->setObjectName(QStringLiteral("label_EEDSmoothing_totTime"));

        horizontalLayout_6->addWidget(label_EEDSmoothing_totTime);

        lineEdit_EEDSmoothing_TotTime = new QLineEdit(groupBox_EEDSmoothing_inputs);
        lineEdit_EEDSmoothing_TotTime->setObjectName(QStringLiteral("lineEdit_EEDSmoothing_TotTime"));

        horizontalLayout_6->addWidget(lineEdit_EEDSmoothing_TotTime);


        verticalLayout->addLayout(horizontalLayout_6);


        verticalLayout_2->addLayout(verticalLayout);

        pushButton_EEDSmoothing_run = new QPushButton(groupBox_EEDSmoothing_inputs);
        pushButton_EEDSmoothing_run->setObjectName(QStringLiteral("pushButton_EEDSmoothing_run"));

        verticalLayout_2->addWidget(pushButton_EEDSmoothing_run);

        widget = new QWidget(EEDSmoothingDialog);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(90, 860, 451, 31));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_EEDSmoothing_saveImName = new QLabel(widget);
        label_EEDSmoothing_saveImName->setObjectName(QStringLiteral("label_EEDSmoothing_saveImName"));

        horizontalLayout->addWidget(label_EEDSmoothing_saveImName);

        lineEdit__EEDSmoothing_saveImgName = new QLineEdit(widget);
        lineEdit__EEDSmoothing_saveImgName->setObjectName(QStringLiteral("lineEdit__EEDSmoothing_saveImgName"));

        horizontalLayout->addWidget(lineEdit__EEDSmoothing_saveImgName);

        pushButton_EEDSmoothing_svaeImgName = new QPushButton(widget);
        pushButton_EEDSmoothing_svaeImgName->setObjectName(QStringLiteral("pushButton_EEDSmoothing_svaeImgName"));

        horizontalLayout->addWidget(pushButton_EEDSmoothing_svaeImgName);


        retranslateUi(EEDSmoothingDialog);

        QMetaObject::connectSlotsByName(EEDSmoothingDialog);
    } // setupUi

    void retranslateUi(QDialog *EEDSmoothingDialog)
    {
        EEDSmoothingDialog->setWindowTitle(QApplication::translate("EEDSmoothingDialog", "EED Smoothing", Q_NULLPTR));
        label_EEDSmoothing_ImgDisplay->setText(QString());
        groupBox_EEDSmoothing_inputs->setTitle(QApplication::translate("EEDSmoothingDialog", "Inputs", Q_NULLPTR));
        label__EEDSmoothing_timeStep->setText(QApplication::translate("EEDSmoothingDialog", "Step Size", Q_NULLPTR));
        label__EEDSmoothing_iterNum->setText(QApplication::translate("EEDSmoothingDialog", "Iteration Number", Q_NULLPTR));
        label__EEDSmoothing_contPar->setText(QApplication::translate("EEDSmoothingDialog", "Contrast Parameter", Q_NULLPTR));
        pushButton__EEDSmoothing_LoadImgData->setText(QApplication::translate("EEDSmoothingDialog", "Load Data", Q_NULLPTR));
        radioButton_EEDSmoothing_explScheme->setText(QApplication::translate("EEDSmoothingDialog", "Explicit Scheme", Q_NULLPTR));
        radioButton_EEDSmoothing_FED->setText(QApplication::translate("EEDSmoothingDialog", "Fast Explicit Scheme(FED)", Q_NULLPTR));
        label_EEDSmoothing_FEDInnerCycles->setText(QApplication::translate("EEDSmoothingDialog", "Number of Inner Cycles", Q_NULLPTR));
        label_EEDSmoothing_totTime->setText(QApplication::translate("EEDSmoothingDialog", "Total Stopping Time", Q_NULLPTR));
        pushButton_EEDSmoothing_run->setText(QApplication::translate("EEDSmoothingDialog", "Run", Q_NULLPTR));
        label_EEDSmoothing_saveImName->setText(QApplication::translate("EEDSmoothingDialog", "Name", Q_NULLPTR));
        pushButton_EEDSmoothing_svaeImgName->setText(QApplication::translate("EEDSmoothingDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class EEDSmoothingDialog: public Ui_EEDSmoothingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EEDSMOOTHINGDIALOG_H
