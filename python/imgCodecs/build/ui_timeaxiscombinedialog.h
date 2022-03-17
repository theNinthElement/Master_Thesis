/********************************************************************************
** Form generated from reading UI file 'timeaxiscombinedialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TIMEAXISCOMBINEDIALOG_H
#define UI_TIMEAXISCOMBINEDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_TimeAxisCombineDialog
{
public:
    QLabel *label_timeComb_imgData1;
    QLabel *label_timeComb_imgData2;
    QGroupBox *groupBox_timeComb_inputs;
    QHBoxLayout *horizontalLayout;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_timeComb_UplData1;
    QPushButton *pushButton_timeComb_UplData2;
    QPushButton *pushButton_outData;
    QFrame *line;
    QPushButton *pushButton_timeComb_combine;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_timeComb_svaeImgName;
    QLineEdit *lineEdit_timeComb_saveImgName;
    QPushButton *pushButton_timeComb_saveImg;

    void setupUi(QDialog *TimeAxisCombineDialog)
    {
        if (TimeAxisCombineDialog->objectName().isEmpty())
            TimeAxisCombineDialog->setObjectName(QStringLiteral("TimeAxisCombineDialog"));
        TimeAxisCombineDialog->resize(1000, 1000);
        label_timeComb_imgData1 = new QLabel(TimeAxisCombineDialog);
        label_timeComb_imgData1->setObjectName(QStringLiteral("label_timeComb_imgData1"));
        label_timeComb_imgData1->setGeometry(QRect(80, 70, 451, 431));
        label_timeComb_imgData2 = new QLabel(TimeAxisCombineDialog);
        label_timeComb_imgData2->setObjectName(QStringLiteral("label_timeComb_imgData2"));
        label_timeComb_imgData2->setGeometry(QRect(510, 530, 421, 411));
        groupBox_timeComb_inputs = new QGroupBox(TimeAxisCombineDialog);
        groupBox_timeComb_inputs->setObjectName(QStringLiteral("groupBox_timeComb_inputs"));
        groupBox_timeComb_inputs->setGeometry(QRect(600, 70, 311, 191));
        horizontalLayout = new QHBoxLayout(groupBox_timeComb_inputs);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        pushButton_timeComb_UplData1 = new QPushButton(groupBox_timeComb_inputs);
        pushButton_timeComb_UplData1->setObjectName(QStringLiteral("pushButton_timeComb_UplData1"));

        verticalLayout->addWidget(pushButton_timeComb_UplData1);

        pushButton_timeComb_UplData2 = new QPushButton(groupBox_timeComb_inputs);
        pushButton_timeComb_UplData2->setObjectName(QStringLiteral("pushButton_timeComb_UplData2"));

        verticalLayout->addWidget(pushButton_timeComb_UplData2);

        pushButton_outData = new QPushButton(groupBox_timeComb_inputs);
        pushButton_outData->setObjectName(QStringLiteral("pushButton_outData"));

        verticalLayout->addWidget(pushButton_outData);

        line = new QFrame(groupBox_timeComb_inputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_timeComb_combine = new QPushButton(groupBox_timeComb_inputs);
        pushButton_timeComb_combine->setObjectName(QStringLiteral("pushButton_timeComb_combine"));

        verticalLayout->addWidget(pushButton_timeComb_combine);


        horizontalLayout->addLayout(verticalLayout);

        widget = new QWidget(TimeAxisCombineDialog);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(70, 880, 401, 31));
        horizontalLayout_2 = new QHBoxLayout(widget);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_timeComb_svaeImgName = new QLabel(widget);
        label_timeComb_svaeImgName->setObjectName(QStringLiteral("label_timeComb_svaeImgName"));

        horizontalLayout_2->addWidget(label_timeComb_svaeImgName);

        lineEdit_timeComb_saveImgName = new QLineEdit(widget);
        lineEdit_timeComb_saveImgName->setObjectName(QStringLiteral("lineEdit_timeComb_saveImgName"));

        horizontalLayout_2->addWidget(lineEdit_timeComb_saveImgName);

        pushButton_timeComb_saveImg = new QPushButton(widget);
        pushButton_timeComb_saveImg->setObjectName(QStringLiteral("pushButton_timeComb_saveImg"));

        horizontalLayout_2->addWidget(pushButton_timeComb_saveImg);


        retranslateUi(TimeAxisCombineDialog);

        QMetaObject::connectSlotsByName(TimeAxisCombineDialog);
    } // setupUi

    void retranslateUi(QDialog *TimeAxisCombineDialog)
    {
        TimeAxisCombineDialog->setWindowTitle(QApplication::translate("TimeAxisCombineDialog", "Time-axis Combine Two Nifti Files", Q_NULLPTR));
        label_timeComb_imgData1->setText(QString());
        label_timeComb_imgData2->setText(QString());
        groupBox_timeComb_inputs->setTitle(QApplication::translate("TimeAxisCombineDialog", "Inputs", Q_NULLPTR));
        pushButton_timeComb_UplData1->setText(QApplication::translate("TimeAxisCombineDialog", "Load Data 1", Q_NULLPTR));
        pushButton_timeComb_UplData2->setText(QApplication::translate("TimeAxisCombineDialog", "Load Data2", Q_NULLPTR));
        pushButton_outData->setText(QApplication::translate("TimeAxisCombineDialog", "Load Output Data", Q_NULLPTR));
        pushButton_timeComb_combine->setText(QApplication::translate("TimeAxisCombineDialog", "Combine", Q_NULLPTR));
        label_timeComb_svaeImgName->setText(QApplication::translate("TimeAxisCombineDialog", "Name", Q_NULLPTR));
        pushButton_timeComb_saveImg->setText(QApplication::translate("TimeAxisCombineDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class TimeAxisCombineDialog: public Ui_TimeAxisCombineDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TIMEAXISCOMBINEDIALOG_H
