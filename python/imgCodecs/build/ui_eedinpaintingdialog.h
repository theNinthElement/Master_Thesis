/********************************************************************************
** Form generated from reading UI file 'eedinpaintingdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_EEDINPAINTINGDIALOG_H
#define UI_EEDINPAINTINGDIALOG_H

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

class Ui_EEDInpaintingDialog
{
public:
    QGroupBox *groupBox_eedInpaint;
    QHBoxLayout *horizontalLayout_7;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout_6;
    QLabel *label_eedInpaintTol;
    QLineEdit *lineEdit_eedInpaintTol;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_eedInpaintTimeStep;
    QLineEdit *lineEdit_eedInpaintTimeStep;
    QRadioButton *radioButton_eedInpaintExplitScheme;
    QRadioButton *radioButton_eedInpaintFED;
    QRadioButton *radioButton_eedInpaintFSI;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_eedInpaintMSETol;
    QLineEdit *lineEdit_eedInpaintMSETol;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_eedInpaintInnerFEDCycSize;
    QLineEdit *lineEdit_eedInpaintInnerFEDCycSize;
    QHBoxLayout *horizontalLayout;
    QLabel *label_eedInpaintFSILoopSize;
    QLineEdit *lineEdit_EEDInpaintFSILoopSoze;
    QPushButton *pushButton_eedInpaintUplImg;
    QPushButton *pushButton_eedInpaintUplMask;
    QFrame *line;
    QPushButton *pushButton_eedInpaintRun;
    QCheckBox *checkBox_EEDInpaintMonitoring;
    QLabel *label_eedInpaintImg;
    QLabel *label_eedInpaintMaskImg;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_eedInpaintFileName;
    QLineEdit *lineEdit_eedInpaintFileName;
    QPushButton *pushButton_eedInpaintSaveImg;
    QCheckBox *checkBox_eedIntp_zeros2mask;
    QCheckBox *checkBox_load_eigenInfo;

    void setupUi(QDialog *EEDInpaintingDialog)
    {
        if (EEDInpaintingDialog->objectName().isEmpty())
            EEDInpaintingDialog->setObjectName(QStringLiteral("EEDInpaintingDialog"));
        EEDInpaintingDialog->resize(1200, 1000);
        groupBox_eedInpaint = new QGroupBox(EEDInpaintingDialog);
        groupBox_eedInpaint->setObjectName(QStringLiteral("groupBox_eedInpaint"));
        groupBox_eedInpaint->setGeometry(QRect(720, 50, 400, 380));
        horizontalLayout_7 = new QHBoxLayout(groupBox_eedInpaint);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        label_eedInpaintTol = new QLabel(groupBox_eedInpaint);
        label_eedInpaintTol->setObjectName(QStringLiteral("label_eedInpaintTol"));

        horizontalLayout_6->addWidget(label_eedInpaintTol);

        lineEdit_eedInpaintTol = new QLineEdit(groupBox_eedInpaint);
        lineEdit_eedInpaintTol->setObjectName(QStringLiteral("lineEdit_eedInpaintTol"));

        horizontalLayout_6->addWidget(lineEdit_eedInpaintTol);


        verticalLayout->addLayout(horizontalLayout_6);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_eedInpaintTimeStep = new QLabel(groupBox_eedInpaint);
        label_eedInpaintTimeStep->setObjectName(QStringLiteral("label_eedInpaintTimeStep"));

        horizontalLayout_2->addWidget(label_eedInpaintTimeStep);

        lineEdit_eedInpaintTimeStep = new QLineEdit(groupBox_eedInpaint);
        lineEdit_eedInpaintTimeStep->setObjectName(QStringLiteral("lineEdit_eedInpaintTimeStep"));

        horizontalLayout_2->addWidget(lineEdit_eedInpaintTimeStep);


        verticalLayout->addLayout(horizontalLayout_2);

        radioButton_eedInpaintExplitScheme = new QRadioButton(groupBox_eedInpaint);
        radioButton_eedInpaintExplitScheme->setObjectName(QStringLiteral("radioButton_eedInpaintExplitScheme"));

        verticalLayout->addWidget(radioButton_eedInpaintExplitScheme);

        radioButton_eedInpaintFED = new QRadioButton(groupBox_eedInpaint);
        radioButton_eedInpaintFED->setObjectName(QStringLiteral("radioButton_eedInpaintFED"));

        verticalLayout->addWidget(radioButton_eedInpaintFED);

        radioButton_eedInpaintFSI = new QRadioButton(groupBox_eedInpaint);
        radioButton_eedInpaintFSI->setObjectName(QStringLiteral("radioButton_eedInpaintFSI"));

        verticalLayout->addWidget(radioButton_eedInpaintFSI);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_eedInpaintMSETol = new QLabel(groupBox_eedInpaint);
        label_eedInpaintMSETol->setObjectName(QStringLiteral("label_eedInpaintMSETol"));

        horizontalLayout_3->addWidget(label_eedInpaintMSETol);

        lineEdit_eedInpaintMSETol = new QLineEdit(groupBox_eedInpaint);
        lineEdit_eedInpaintMSETol->setObjectName(QStringLiteral("lineEdit_eedInpaintMSETol"));

        horizontalLayout_3->addWidget(lineEdit_eedInpaintMSETol);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_eedInpaintInnerFEDCycSize = new QLabel(groupBox_eedInpaint);
        label_eedInpaintInnerFEDCycSize->setObjectName(QStringLiteral("label_eedInpaintInnerFEDCycSize"));

        horizontalLayout_5->addWidget(label_eedInpaintInnerFEDCycSize);

        lineEdit_eedInpaintInnerFEDCycSize = new QLineEdit(groupBox_eedInpaint);
        lineEdit_eedInpaintInnerFEDCycSize->setObjectName(QStringLiteral("lineEdit_eedInpaintInnerFEDCycSize"));

        horizontalLayout_5->addWidget(lineEdit_eedInpaintInnerFEDCycSize);


        verticalLayout->addLayout(horizontalLayout_5);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_eedInpaintFSILoopSize = new QLabel(groupBox_eedInpaint);
        label_eedInpaintFSILoopSize->setObjectName(QStringLiteral("label_eedInpaintFSILoopSize"));

        horizontalLayout->addWidget(label_eedInpaintFSILoopSize);

        lineEdit_EEDInpaintFSILoopSoze = new QLineEdit(groupBox_eedInpaint);
        lineEdit_EEDInpaintFSILoopSoze->setObjectName(QStringLiteral("lineEdit_EEDInpaintFSILoopSoze"));

        horizontalLayout->addWidget(lineEdit_EEDInpaintFSILoopSoze);


        verticalLayout->addLayout(horizontalLayout);

        pushButton_eedInpaintUplImg = new QPushButton(groupBox_eedInpaint);
        pushButton_eedInpaintUplImg->setObjectName(QStringLiteral("pushButton_eedInpaintUplImg"));

        verticalLayout->addWidget(pushButton_eedInpaintUplImg);

        pushButton_eedInpaintUplMask = new QPushButton(groupBox_eedInpaint);
        pushButton_eedInpaintUplMask->setObjectName(QStringLiteral("pushButton_eedInpaintUplMask"));

        verticalLayout->addWidget(pushButton_eedInpaintUplMask);

        line = new QFrame(groupBox_eedInpaint);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_eedInpaintRun = new QPushButton(groupBox_eedInpaint);
        pushButton_eedInpaintRun->setObjectName(QStringLiteral("pushButton_eedInpaintRun"));

        verticalLayout->addWidget(pushButton_eedInpaintRun);

        checkBox_EEDInpaintMonitoring = new QCheckBox(groupBox_eedInpaint);
        checkBox_EEDInpaintMonitoring->setObjectName(QStringLiteral("checkBox_EEDInpaintMonitoring"));

        verticalLayout->addWidget(checkBox_EEDInpaintMonitoring);


        horizontalLayout_7->addLayout(verticalLayout);

        label_eedInpaintImg = new QLabel(EEDInpaintingDialog);
        label_eedInpaintImg->setObjectName(QStringLiteral("label_eedInpaintImg"));
        label_eedInpaintImg->setGeometry(QRect(80, 60, 500, 500));
        label_eedInpaintMaskImg = new QLabel(EEDInpaintingDialog);
        label_eedInpaintMaskImg->setObjectName(QStringLiteral("label_eedInpaintMaskImg"));
        label_eedInpaintMaskImg->setGeometry(QRect(680, 490, 450, 450));
        layoutWidget = new QWidget(EEDInpaintingDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(90, 860, 471, 31));
        horizontalLayout_4 = new QHBoxLayout(layoutWidget);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, 0, 0);
        label_eedInpaintFileName = new QLabel(layoutWidget);
        label_eedInpaintFileName->setObjectName(QStringLiteral("label_eedInpaintFileName"));

        horizontalLayout_4->addWidget(label_eedInpaintFileName);

        lineEdit_eedInpaintFileName = new QLineEdit(layoutWidget);
        lineEdit_eedInpaintFileName->setObjectName(QStringLiteral("lineEdit_eedInpaintFileName"));

        horizontalLayout_4->addWidget(lineEdit_eedInpaintFileName);

        pushButton_eedInpaintSaveImg = new QPushButton(layoutWidget);
        pushButton_eedInpaintSaveImg->setObjectName(QStringLiteral("pushButton_eedInpaintSaveImg"));

        horizontalLayout_4->addWidget(pushButton_eedInpaintSaveImg);

        checkBox_eedIntp_zeros2mask = new QCheckBox(EEDInpaintingDialog);
        checkBox_eedIntp_zeros2mask->setObjectName(QStringLiteral("checkBox_eedIntp_zeros2mask"));
        checkBox_eedIntp_zeros2mask->setGeometry(QRect(740, 430, 121, 22));
        checkBox_load_eigenInfo = new QCheckBox(EEDInpaintingDialog);
        checkBox_load_eigenInfo->setObjectName(QStringLiteral("checkBox_load_eigenInfo"));
        checkBox_load_eigenInfo->setGeometry(QRect(740, 470, 171, 22));

        retranslateUi(EEDInpaintingDialog);

        QMetaObject::connectSlotsByName(EEDInpaintingDialog);
    } // setupUi

    void retranslateUi(QDialog *EEDInpaintingDialog)
    {
        EEDInpaintingDialog->setWindowTitle(QApplication::translate("EEDInpaintingDialog", "Inpainting by Edge Enhancing Anisotropy Diffusion", Q_NULLPTR));
        groupBox_eedInpaint->setTitle(QApplication::translate("EEDInpaintingDialog", "Inputs", Q_NULLPTR));
        label_eedInpaintTol->setText(QApplication::translate("EEDInpaintingDialog", "Tolerance", Q_NULLPTR));
        label_eedInpaintTimeStep->setText(QApplication::translate("EEDInpaintingDialog", "Time Step Size", Q_NULLPTR));
        radioButton_eedInpaintExplitScheme->setText(QApplication::translate("EEDInpaintingDialog", "Explicit Scheme", Q_NULLPTR));
        radioButton_eedInpaintFED->setText(QApplication::translate("EEDInpaintingDialog", "Fast Explicit Diffusion (FED) Scheme", Q_NULLPTR));
        radioButton_eedInpaintFSI->setText(QApplication::translate("EEDInpaintingDialog", "Fast Semi_iterative Scheme (FSI)", Q_NULLPTR));
        label_eedInpaintMSETol->setText(QApplication::translate("EEDInpaintingDialog", "MSE Tolerance", Q_NULLPTR));
        label_eedInpaintInnerFEDCycSize->setText(QApplication::translate("EEDInpaintingDialog", "Inner Cycle Size", Q_NULLPTR));
        label_eedInpaintFSILoopSize->setText(QApplication::translate("EEDInpaintingDialog", "Inner Loop Size", Q_NULLPTR));
        pushButton_eedInpaintUplImg->setText(QApplication::translate("EEDInpaintingDialog", "Load Image Data", Q_NULLPTR));
        pushButton_eedInpaintUplMask->setText(QApplication::translate("EEDInpaintingDialog", "Load Mask Data", Q_NULLPTR));
        pushButton_eedInpaintRun->setText(QApplication::translate("EEDInpaintingDialog", "Run", Q_NULLPTR));
        checkBox_EEDInpaintMonitoring->setText(QApplication::translate("EEDInpaintingDialog", "Monitoring", Q_NULLPTR));
        label_eedInpaintImg->setText(QString());
        label_eedInpaintMaskImg->setText(QString());
        label_eedInpaintFileName->setText(QApplication::translate("EEDInpaintingDialog", "File Name", Q_NULLPTR));
        pushButton_eedInpaintSaveImg->setText(QApplication::translate("EEDInpaintingDialog", "Save", Q_NULLPTR));
        checkBox_eedIntp_zeros2mask->setText(QApplication::translate("EEDInpaintingDialog", "Zeros to Mask", Q_NULLPTR));
        checkBox_load_eigenInfo->setText(QApplication::translate("EEDInpaintingDialog", "Load eigen information", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class EEDInpaintingDialog: public Ui_EEDInpaintingDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_EEDINPAINTINGDIALOG_H
