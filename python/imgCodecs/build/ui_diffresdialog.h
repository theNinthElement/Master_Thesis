/********************************************************************************
** Form generated from reading UI file 'diffresdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIFFRESDIALOG_H
#define UI_DIFFRESDIALOG_H

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

class Ui_DiffResDialog
{
public:
    QLabel *label_diffRes_recImg;
    QLabel *label_diffRes_OrigImg;
    QGroupBox *groupBox_diffRes_inputs;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout;
    QPushButton *pushButton_diffRes_origImgLoad;
    QPushButton *pushButton_diffRes_recImgLoad;
    QFrame *line;
    QPushButton *pushButton_diffRes_run;
    QWidget *widget;
    QHBoxLayout *horizontalLayout;
    QLabel *label_diffRes_saveFileName;
    QLineEdit *lineEdit_diffRes_ImgName;
    QPushButton *pushButton_diffRes_saveImg;

    void setupUi(QDialog *DiffResDialog)
    {
        if (DiffResDialog->objectName().isEmpty())
            DiffResDialog->setObjectName(QStringLiteral("DiffResDialog"));
        DiffResDialog->resize(1000, 1000);
        label_diffRes_recImg = new QLabel(DiffResDialog);
        label_diffRes_recImg->setObjectName(QStringLiteral("label_diffRes_recImg"));
        label_diffRes_recImg->setGeometry(QRect(70, 80, 420, 420));
        label_diffRes_OrigImg = new QLabel(DiffResDialog);
        label_diffRes_OrigImg->setObjectName(QStringLiteral("label_diffRes_OrigImg"));
        label_diffRes_OrigImg->setGeometry(QRect(560, 540, 390, 390));
        groupBox_diffRes_inputs = new QGroupBox(DiffResDialog);
        groupBox_diffRes_inputs->setObjectName(QStringLiteral("groupBox_diffRes_inputs"));
        groupBox_diffRes_inputs->setGeometry(QRect(690, 70, 200, 200));
        horizontalLayout_2 = new QHBoxLayout(groupBox_diffRes_inputs);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        pushButton_diffRes_origImgLoad = new QPushButton(groupBox_diffRes_inputs);
        pushButton_diffRes_origImgLoad->setObjectName(QStringLiteral("pushButton_diffRes_origImgLoad"));

        verticalLayout->addWidget(pushButton_diffRes_origImgLoad);

        pushButton_diffRes_recImgLoad = new QPushButton(groupBox_diffRes_inputs);
        pushButton_diffRes_recImgLoad->setObjectName(QStringLiteral("pushButton_diffRes_recImgLoad"));

        verticalLayout->addWidget(pushButton_diffRes_recImgLoad);

        line = new QFrame(groupBox_diffRes_inputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_diffRes_run = new QPushButton(groupBox_diffRes_inputs);
        pushButton_diffRes_run->setObjectName(QStringLiteral("pushButton_diffRes_run"));

        verticalLayout->addWidget(pushButton_diffRes_run);


        horizontalLayout_2->addLayout(verticalLayout);

        widget = new QWidget(DiffResDialog);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(90, 870, 401, 31));
        horizontalLayout = new QHBoxLayout(widget);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        label_diffRes_saveFileName = new QLabel(widget);
        label_diffRes_saveFileName->setObjectName(QStringLiteral("label_diffRes_saveFileName"));

        horizontalLayout->addWidget(label_diffRes_saveFileName);

        lineEdit_diffRes_ImgName = new QLineEdit(widget);
        lineEdit_diffRes_ImgName->setObjectName(QStringLiteral("lineEdit_diffRes_ImgName"));

        horizontalLayout->addWidget(lineEdit_diffRes_ImgName);

        pushButton_diffRes_saveImg = new QPushButton(widget);
        pushButton_diffRes_saveImg->setObjectName(QStringLiteral("pushButton_diffRes_saveImg"));

        horizontalLayout->addWidget(pushButton_diffRes_saveImg);


        retranslateUi(DiffResDialog);

        QMetaObject::connectSlotsByName(DiffResDialog);
    } // setupUi

    void retranslateUi(QDialog *DiffResDialog)
    {
        DiffResDialog->setWindowTitle(QApplication::translate("DiffResDialog", "Difference type residual", Q_NULLPTR));
        label_diffRes_recImg->setText(QString());
        label_diffRes_OrigImg->setText(QString());
        groupBox_diffRes_inputs->setTitle(QApplication::translate("DiffResDialog", "Input", Q_NULLPTR));
        pushButton_diffRes_origImgLoad->setText(QApplication::translate("DiffResDialog", "Original Data", Q_NULLPTR));
        pushButton_diffRes_recImgLoad->setText(QApplication::translate("DiffResDialog", "Reconstructed Data", Q_NULLPTR));
        pushButton_diffRes_run->setText(QApplication::translate("DiffResDialog", "Run", Q_NULLPTR));
        label_diffRes_saveFileName->setText(QApplication::translate("DiffResDialog", "File Name", Q_NULLPTR));
        pushButton_diffRes_saveImg->setText(QApplication::translate("DiffResDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class DiffResDialog: public Ui_DiffResDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIFFRESDIALOG_H
