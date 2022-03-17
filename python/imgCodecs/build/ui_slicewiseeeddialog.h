/********************************************************************************
** Form generated from reading UI file 'slicewiseeeddialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_SLICEWISEEEDDIALOG_H
#define UI_SLICEWISEEEDDIALOG_H

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
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_SlicewiseEEDDialog
{
public:
    QGroupBox *groupBox_2dEED_inputs;
    QWidget *widget;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label_2deed_tol;
    QLineEdit *lineEdit_2deed_tol;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2deed_timeStep;
    QLineEdit *lineEdit_2deed_timeStep;
    QRadioButton *radioButton_2deed_ExScheme;
    QRadioButton *radioButton_FEDScheme;
    QRadioButton *radioButton_FSIScheme;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_mse_tol;
    QLineEdit *lineEdit_mse_tol;
    QHBoxLayout *horizontalLayout_4;
    QLabel *label_innerCycleSize;
    QLineEdit *lineEdit_innerCycleSize;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_innerLoopSIze;
    QLineEdit *lineEdit_innerLoopSIze;
    QPushButton *pushButton_loadImg;
    QPushButton *pushButton_mask;
    QFrame *line_2dEED;
    QPushButton *pushButton_run2dEED;
    QHBoxLayout *horizontalLayout_6;
    QCheckBox *checkBox_monitoring;
    QSpinBox *spinBox_2dEED;
    QLabel *label_imgShow;
    QLabel *label_maskShow;
    QWidget *widget1;
    QHBoxLayout *horizontalLayout_7;
    QLabel *label_saveImgName;
    QLineEdit *lineEdit_saveImgName;
    QPushButton *pushButton_saveImgName;

    void setupUi(QDialog *SlicewiseEEDDialog)
    {
        if (SlicewiseEEDDialog->objectName().isEmpty())
            SlicewiseEEDDialog->setObjectName(QStringLiteral("SlicewiseEEDDialog"));
        SlicewiseEEDDialog->resize(1281, 1142);
        groupBox_2dEED_inputs = new QGroupBox(SlicewiseEEDDialog);
        groupBox_2dEED_inputs->setObjectName(QStringLiteral("groupBox_2dEED_inputs"));
        groupBox_2dEED_inputs->setGeometry(QRect(820, 30, 381, 501));
        widget = new QWidget(groupBox_2dEED_inputs);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(20, 40, 341, 441));
        verticalLayout = new QVBoxLayout(widget);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_2deed_tol = new QLabel(widget);
        label_2deed_tol->setObjectName(QStringLiteral("label_2deed_tol"));

        horizontalLayout->addWidget(label_2deed_tol);

        lineEdit_2deed_tol = new QLineEdit(widget);
        lineEdit_2deed_tol->setObjectName(QStringLiteral("lineEdit_2deed_tol"));

        horizontalLayout->addWidget(lineEdit_2deed_tol);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_2deed_timeStep = new QLabel(widget);
        label_2deed_timeStep->setObjectName(QStringLiteral("label_2deed_timeStep"));

        horizontalLayout_2->addWidget(label_2deed_timeStep);

        lineEdit_2deed_timeStep = new QLineEdit(widget);
        lineEdit_2deed_timeStep->setObjectName(QStringLiteral("lineEdit_2deed_timeStep"));

        horizontalLayout_2->addWidget(lineEdit_2deed_timeStep);


        verticalLayout->addLayout(horizontalLayout_2);

        radioButton_2deed_ExScheme = new QRadioButton(widget);
        radioButton_2deed_ExScheme->setObjectName(QStringLiteral("radioButton_2deed_ExScheme"));

        verticalLayout->addWidget(radioButton_2deed_ExScheme);

        radioButton_FEDScheme = new QRadioButton(widget);
        radioButton_FEDScheme->setObjectName(QStringLiteral("radioButton_FEDScheme"));

        verticalLayout->addWidget(radioButton_FEDScheme);

        radioButton_FSIScheme = new QRadioButton(widget);
        radioButton_FSIScheme->setObjectName(QStringLiteral("radioButton_FSIScheme"));

        verticalLayout->addWidget(radioButton_FSIScheme);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        label_mse_tol = new QLabel(widget);
        label_mse_tol->setObjectName(QStringLiteral("label_mse_tol"));

        horizontalLayout_3->addWidget(label_mse_tol);

        lineEdit_mse_tol = new QLineEdit(widget);
        lineEdit_mse_tol->setObjectName(QStringLiteral("lineEdit_mse_tol"));

        horizontalLayout_3->addWidget(lineEdit_mse_tol);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        label_innerCycleSize = new QLabel(widget);
        label_innerCycleSize->setObjectName(QStringLiteral("label_innerCycleSize"));

        horizontalLayout_4->addWidget(label_innerCycleSize);

        lineEdit_innerCycleSize = new QLineEdit(widget);
        lineEdit_innerCycleSize->setObjectName(QStringLiteral("lineEdit_innerCycleSize"));

        horizontalLayout_4->addWidget(lineEdit_innerCycleSize);


        verticalLayout->addLayout(horizontalLayout_4);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QStringLiteral("horizontalLayout_5"));
        label_innerLoopSIze = new QLabel(widget);
        label_innerLoopSIze->setObjectName(QStringLiteral("label_innerLoopSIze"));

        horizontalLayout_5->addWidget(label_innerLoopSIze);

        lineEdit_innerLoopSIze = new QLineEdit(widget);
        lineEdit_innerLoopSIze->setObjectName(QStringLiteral("lineEdit_innerLoopSIze"));

        horizontalLayout_5->addWidget(lineEdit_innerLoopSIze);


        verticalLayout->addLayout(horizontalLayout_5);

        pushButton_loadImg = new QPushButton(widget);
        pushButton_loadImg->setObjectName(QStringLiteral("pushButton_loadImg"));

        verticalLayout->addWidget(pushButton_loadImg);

        pushButton_mask = new QPushButton(widget);
        pushButton_mask->setObjectName(QStringLiteral("pushButton_mask"));

        verticalLayout->addWidget(pushButton_mask);

        line_2dEED = new QFrame(widget);
        line_2dEED->setObjectName(QStringLiteral("line_2dEED"));
        line_2dEED->setFrameShape(QFrame::HLine);
        line_2dEED->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line_2dEED);

        pushButton_run2dEED = new QPushButton(widget);
        pushButton_run2dEED->setObjectName(QStringLiteral("pushButton_run2dEED"));

        verticalLayout->addWidget(pushButton_run2dEED);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QStringLiteral("horizontalLayout_6"));
        checkBox_monitoring = new QCheckBox(widget);
        checkBox_monitoring->setObjectName(QStringLiteral("checkBox_monitoring"));

        horizontalLayout_6->addWidget(checkBox_monitoring);

        spinBox_2dEED = new QSpinBox(widget);
        spinBox_2dEED->setObjectName(QStringLiteral("spinBox_2dEED"));

        horizontalLayout_6->addWidget(spinBox_2dEED);


        verticalLayout->addLayout(horizontalLayout_6);

        label_imgShow = new QLabel(SlicewiseEEDDialog);
        label_imgShow->setObjectName(QStringLiteral("label_imgShow"));
        label_imgShow->setGeometry(QRect(100, 60, 541, 531));
        label_maskShow = new QLabel(SlicewiseEEDDialog);
        label_maskShow->setObjectName(QStringLiteral("label_maskShow"));
        label_maskShow->setGeometry(QRect(700, 600, 521, 461));
        widget1 = new QWidget(SlicewiseEEDDialog);
        widget1->setObjectName(QStringLiteral("widget1"));
        widget1->setGeometry(QRect(150, 1020, 421, 31));
        horizontalLayout_7 = new QHBoxLayout(widget1);
        horizontalLayout_7->setObjectName(QStringLiteral("horizontalLayout_7"));
        horizontalLayout_7->setContentsMargins(0, 0, 0, 0);
        label_saveImgName = new QLabel(widget1);
        label_saveImgName->setObjectName(QStringLiteral("label_saveImgName"));

        horizontalLayout_7->addWidget(label_saveImgName);

        lineEdit_saveImgName = new QLineEdit(widget1);
        lineEdit_saveImgName->setObjectName(QStringLiteral("lineEdit_saveImgName"));

        horizontalLayout_7->addWidget(lineEdit_saveImgName);

        pushButton_saveImgName = new QPushButton(widget1);
        pushButton_saveImgName->setObjectName(QStringLiteral("pushButton_saveImgName"));

        horizontalLayout_7->addWidget(pushButton_saveImgName);


        retranslateUi(SlicewiseEEDDialog);

        QMetaObject::connectSlotsByName(SlicewiseEEDDialog);
    } // setupUi

    void retranslateUi(QDialog *SlicewiseEEDDialog)
    {
        SlicewiseEEDDialog->setWindowTitle(QApplication::translate("SlicewiseEEDDialog", "Slice-wise EED Inpainting", Q_NULLPTR));
        groupBox_2dEED_inputs->setTitle(QApplication::translate("SlicewiseEEDDialog", "Inputs", Q_NULLPTR));
        label_2deed_tol->setText(QApplication::translate("SlicewiseEEDDialog", "Tolerance", Q_NULLPTR));
        label_2deed_timeStep->setText(QApplication::translate("SlicewiseEEDDialog", "Time Step Size", Q_NULLPTR));
        radioButton_2deed_ExScheme->setText(QApplication::translate("SlicewiseEEDDialog", "Explicit Scheme", Q_NULLPTR));
        radioButton_FEDScheme->setText(QApplication::translate("SlicewiseEEDDialog", "Fast Explicit Scheme (FED)", Q_NULLPTR));
        radioButton_FSIScheme->setText(QApplication::translate("SlicewiseEEDDialog", "Fast Semi Iterative Scheme (FSI)", Q_NULLPTR));
        label_mse_tol->setText(QApplication::translate("SlicewiseEEDDialog", "MSE Tolerance", Q_NULLPTR));
        label_innerCycleSize->setText(QApplication::translate("SlicewiseEEDDialog", "Inner Cycle Size", Q_NULLPTR));
        label_innerLoopSIze->setText(QApplication::translate("SlicewiseEEDDialog", "Inner Loop Size", Q_NULLPTR));
        pushButton_loadImg->setText(QApplication::translate("SlicewiseEEDDialog", "Image", Q_NULLPTR));
        pushButton_mask->setText(QApplication::translate("SlicewiseEEDDialog", "Mask", Q_NULLPTR));
        pushButton_run2dEED->setText(QApplication::translate("SlicewiseEEDDialog", "Run", Q_NULLPTR));
        checkBox_monitoring->setText(QApplication::translate("SlicewiseEEDDialog", "Monitoring", Q_NULLPTR));
        label_imgShow->setText(QString());
        label_maskShow->setText(QString());
        label_saveImgName->setText(QApplication::translate("SlicewiseEEDDialog", "Name", Q_NULLPTR));
        pushButton_saveImgName->setText(QApplication::translate("SlicewiseEEDDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class SlicewiseEEDDialog: public Ui_SlicewiseEEDDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_SLICEWISEEEDDIALOG_H
