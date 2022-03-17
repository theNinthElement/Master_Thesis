/********************************************************************************
** Form generated from reading UI file 'xorresdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_XORRESDIALOG_H
#define UI_XORRESDIALOG_H

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

class Ui_XORResDialog
{
public:
    QLabel *label_xorRes_recImg;
    QLabel *label_xorRes_OrigImg;
    QGroupBox *groupBox;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_xorRes_origImgLoad;
    QPushButton *pushButton_xorRes_recImgLoad;
    QFrame *line_xorRes;
    QPushButton *pushButton_xorRes_run;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_xorRes_savedImgName;
    QLineEdit *lineEdit_xorRes_ImgName;
    QPushButton *pushButton_xorres_saveImg;

    void setupUi(QDialog *XORResDialog)
    {
        if (XORResDialog->objectName().isEmpty())
            XORResDialog->setObjectName(QStringLiteral("XORResDialog"));
        XORResDialog->resize(1000, 1000);
        label_xorRes_recImg = new QLabel(XORResDialog);
        label_xorRes_recImg->setObjectName(QStringLiteral("label_xorRes_recImg"));
        label_xorRes_recImg->setGeometry(QRect(100, 90, 420, 420));
        label_xorRes_OrigImg = new QLabel(XORResDialog);
        label_xorRes_OrigImg->setObjectName(QStringLiteral("label_xorRes_OrigImg"));
        label_xorRes_OrigImg->setGeometry(QRect(560, 460, 390, 390));
        groupBox = new QGroupBox(XORResDialog);
        groupBox->setObjectName(QStringLiteral("groupBox"));
        groupBox->setGeometry(QRect(680, 90, 200, 200));
        horizontalLayout_2 = new QHBoxLayout(groupBox);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        pushButton_xorRes_origImgLoad = new QPushButton(groupBox);
        pushButton_xorRes_origImgLoad->setObjectName(QStringLiteral("pushButton_xorRes_origImgLoad"));

        verticalLayout->addWidget(pushButton_xorRes_origImgLoad);

        pushButton_xorRes_recImgLoad = new QPushButton(groupBox);
        pushButton_xorRes_recImgLoad->setObjectName(QStringLiteral("pushButton_xorRes_recImgLoad"));

        verticalLayout->addWidget(pushButton_xorRes_recImgLoad);

        line_xorRes = new QFrame(groupBox);
        line_xorRes->setObjectName(QStringLiteral("line_xorRes"));
        line_xorRes->setFrameShape(QFrame::HLine);
        line_xorRes->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line_xorRes);

        pushButton_xorRes_run = new QPushButton(groupBox);
        pushButton_xorRes_run->setObjectName(QStringLiteral("pushButton_xorRes_run"));

        verticalLayout->addWidget(pushButton_xorRes_run);


        horizontalLayout_2->addLayout(verticalLayout);

        widget = new QWidget(XORResDialog);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(80, 860, 411, 31));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_xorRes_savedImgName = new QLabel(widget);
        label_xorRes_savedImgName->setObjectName(QStringLiteral("label_xorRes_savedImgName"));

        horizontalLayout->addWidget(label_xorRes_savedImgName);

        lineEdit_xorRes_ImgName = new QLineEdit(widget);
        lineEdit_xorRes_ImgName->setObjectName(QStringLiteral("lineEdit_xorRes_ImgName"));

        horizontalLayout->addWidget(lineEdit_xorRes_ImgName);

        pushButton_xorres_saveImg = new QPushButton(widget);
        pushButton_xorres_saveImg->setObjectName(QStringLiteral("pushButton_xorres_saveImg"));

        horizontalLayout->addWidget(pushButton_xorres_saveImg);


        retranslateUi(XORResDialog);

        QMetaObject::connectSlotsByName(XORResDialog);
    } // setupUi

    void retranslateUi(QDialog *XORResDialog)
    {
        XORResDialog->setWindowTitle(QApplication::translate("XORResDialog", "XOR type residual", Q_NULLPTR));
        label_xorRes_recImg->setText(QString());
        label_xorRes_OrigImg->setText(QString());
        groupBox->setTitle(QApplication::translate("XORResDialog", "Inputs", Q_NULLPTR));
        pushButton_xorRes_origImgLoad->setText(QApplication::translate("XORResDialog", "Original Data", Q_NULLPTR));
        pushButton_xorRes_recImgLoad->setText(QApplication::translate("XORResDialog", "Reconstructed Data", Q_NULLPTR));
        pushButton_xorRes_run->setText(QApplication::translate("XORResDialog", "Run", Q_NULLPTR));
        label_xorRes_savedImgName->setText(QApplication::translate("XORResDialog", "File Name", Q_NULLPTR));
        pushButton_xorres_saveImg->setText(QApplication::translate("XORResDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class XORResDialog: public Ui_XORResDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_XORRESDIALOG_H
