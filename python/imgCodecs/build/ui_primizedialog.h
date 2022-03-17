/********************************************************************************
** Form generated from reading UI file 'primizedialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PRIMIZEDIALOG_H
#define UI_PRIMIZEDIALOG_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PrimizeDialog
{
public:
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_primImg_ImgName;
    QLineEdit *lineEdit_primImgName;
    QPushButton *pushButton_primImgSave;
    QLabel *label_primImg_origImg;
    QGroupBox *groupBox_primImgInputs;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout;
    QComboBox *comboBox_primImg_ImgDim;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButton_primImg_imgUpl;
    QSpacerItem *verticalSpacer_2;
    QFrame *line;
    QPushButton *pushButton_primImg_run;

    void setupUi(QDialog *PrimizeDialog)
    {
        if (PrimizeDialog->objectName().isEmpty())
            PrimizeDialog->setObjectName(QStringLiteral("PrimizeDialog"));
        PrimizeDialog->resize(1239, 1012);
        layoutWidget = new QWidget(PrimizeDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(100, 840, 401, 31));
        horizontalLayout_3 = new QHBoxLayout(layoutWidget);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_primImg_ImgName = new QLabel(layoutWidget);
        label_primImg_ImgName->setObjectName(QStringLiteral("label_primImg_ImgName"));

        horizontalLayout_3->addWidget(label_primImg_ImgName);

        lineEdit_primImgName = new QLineEdit(layoutWidget);
        lineEdit_primImgName->setObjectName(QStringLiteral("lineEdit_primImgName"));

        horizontalLayout_3->addWidget(lineEdit_primImgName);

        pushButton_primImgSave = new QPushButton(layoutWidget);
        pushButton_primImgSave->setObjectName(QStringLiteral("pushButton_primImgSave"));

        horizontalLayout_3->addWidget(pushButton_primImgSave);

        label_primImg_origImg = new QLabel(PrimizeDialog);
        label_primImg_origImg->setObjectName(QStringLiteral("label_primImg_origImg"));
        label_primImg_origImg->setGeometry(QRect(90, 70, 420, 420));
        groupBox_primImgInputs = new QGroupBox(PrimizeDialog);
        groupBox_primImgInputs->setObjectName(QStringLiteral("groupBox_primImgInputs"));
        groupBox_primImgInputs->setGeometry(QRect(920, 80, 220, 220));
        horizontalLayout_2 = new QHBoxLayout(groupBox_primImgInputs);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        comboBox_primImg_ImgDim = new QComboBox(groupBox_primImgInputs);
        comboBox_primImg_ImgDim->setObjectName(QStringLiteral("comboBox_primImg_ImgDim"));

        verticalLayout->addWidget(comboBox_primImg_ImgDim);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));

        verticalLayout->addLayout(horizontalLayout);

        pushButton_primImg_imgUpl = new QPushButton(groupBox_primImgInputs);
        pushButton_primImg_imgUpl->setObjectName(QStringLiteral("pushButton_primImg_imgUpl"));

        verticalLayout->addWidget(pushButton_primImg_imgUpl);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_2);

        line = new QFrame(groupBox_primImgInputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_primImg_run = new QPushButton(groupBox_primImgInputs);
        pushButton_primImg_run->setObjectName(QStringLiteral("pushButton_primImg_run"));

        verticalLayout->addWidget(pushButton_primImg_run);


        horizontalLayout_2->addLayout(verticalLayout);


        retranslateUi(PrimizeDialog);

        QMetaObject::connectSlotsByName(PrimizeDialog);
    } // setupUi

    void retranslateUi(QDialog *PrimizeDialog)
    {
        PrimizeDialog->setWindowTitle(QApplication::translate("PrimizeDialog", "Image Primize Dialog", Q_NULLPTR));
        label_primImg_ImgName->setText(QApplication::translate("PrimizeDialog", "File Name", Q_NULLPTR));
        pushButton_primImgSave->setText(QApplication::translate("PrimizeDialog", "Save", Q_NULLPTR));
        label_primImg_origImg->setText(QString());
        groupBox_primImgInputs->setTitle(QApplication::translate("PrimizeDialog", "Inputs", Q_NULLPTR));
        comboBox_primImg_ImgDim->clear();
        comboBox_primImg_ImgDim->insertItems(0, QStringList()
         << QApplication::translate("PrimizeDialog", "2D Data", Q_NULLPTR)
         << QApplication::translate("PrimizeDialog", "3D Data", Q_NULLPTR)
         << QApplication::translate("PrimizeDialog", "4D Data", Q_NULLPTR)
        );
        pushButton_primImg_imgUpl->setText(QApplication::translate("PrimizeDialog", "Open Image Data", Q_NULLPTR));
        pushButton_primImg_run->setText(QApplication::translate("PrimizeDialog", "Run", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class PrimizeDialog: public Ui_PrimizeDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PRIMIZEDIALOG_H
