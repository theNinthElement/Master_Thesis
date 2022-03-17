/********************************************************************************
** Form generated from reading UI file 'reggridmaskdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_REGGRIDMASKDIALOG_H
#define UI_REGGRIDMASKDIALOG_H

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
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_RegGridMaskDialog
{
public:
    QLabel *label_regGridMaskProcImg;
    QGroupBox *groupBox_regMaskInputs;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout;
    QComboBox *comboBox_regMask_ImgDim;
    QHBoxLayout *horizontalLayout;
    QLabel *label_regMask_percentage;
    QLineEdit *lineEdit_regMask_percentage;
    QPushButton *pushButton_regMask_imgUpl;
    QSpacerItem *verticalSpacer_2;
    QFrame *line;
    QPushButton *pushButton_regMask_run;
    QLabel *label_regMaskOrigImg;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_regmask_ImgName;
    QLineEdit *lineEdit_regMaskImgName;
    QPushButton *pushButton_regMaskImgSave;
    QSlider *horizontalSlider_regGridMaskProcImg;
    QSlider *horizontalSlider_regGridOrigImg;

    void setupUi(QDialog *RegGridMaskDialog)
    {
        if (RegGridMaskDialog->objectName().isEmpty())
            RegGridMaskDialog->setObjectName(QStringLiteral("RegGridMaskDialog"));
        RegGridMaskDialog->resize(1000, 1000);
        label_regGridMaskProcImg = new QLabel(RegGridMaskDialog);
        label_regGridMaskProcImg->setObjectName(QStringLiteral("label_regGridMaskProcImg"));
        label_regGridMaskProcImg->setGeometry(QRect(80, 70, 420, 420));
        groupBox_regMaskInputs = new QGroupBox(RegGridMaskDialog);
        groupBox_regMaskInputs->setObjectName(QStringLiteral("groupBox_regMaskInputs"));
        groupBox_regMaskInputs->setGeometry(QRect(670, 70, 220, 220));
        horizontalLayout_2 = new QHBoxLayout(groupBox_regMaskInputs);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        comboBox_regMask_ImgDim = new QComboBox(groupBox_regMaskInputs);
        comboBox_regMask_ImgDim->setObjectName(QStringLiteral("comboBox_regMask_ImgDim"));

        verticalLayout->addWidget(comboBox_regMask_ImgDim);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_regMask_percentage = new QLabel(groupBox_regMaskInputs);
        label_regMask_percentage->setObjectName(QStringLiteral("label_regMask_percentage"));

        horizontalLayout->addWidget(label_regMask_percentage);

        lineEdit_regMask_percentage = new QLineEdit(groupBox_regMaskInputs);
        lineEdit_regMask_percentage->setObjectName(QStringLiteral("lineEdit_regMask_percentage"));

        horizontalLayout->addWidget(lineEdit_regMask_percentage);


        verticalLayout->addLayout(horizontalLayout);

        pushButton_regMask_imgUpl = new QPushButton(groupBox_regMaskInputs);
        pushButton_regMask_imgUpl->setObjectName(QStringLiteral("pushButton_regMask_imgUpl"));

        verticalLayout->addWidget(pushButton_regMask_imgUpl);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_2);

        line = new QFrame(groupBox_regMaskInputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_regMask_run = new QPushButton(groupBox_regMaskInputs);
        pushButton_regMask_run->setObjectName(QStringLiteral("pushButton_regMask_run"));

        verticalLayout->addWidget(pushButton_regMask_run);


        horizontalLayout_2->addLayout(verticalLayout);

        label_regMaskOrigImg = new QLabel(RegGridMaskDialog);
        label_regMaskOrigImg->setObjectName(QStringLiteral("label_regMaskOrigImg"));
        label_regMaskOrigImg->setGeometry(QRect(540, 480, 390, 390));
        layoutWidget = new QWidget(RegGridMaskDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(90, 860, 401, 31));
        horizontalLayout_3 = new QHBoxLayout(layoutWidget);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_regmask_ImgName = new QLabel(layoutWidget);
        label_regmask_ImgName->setObjectName(QStringLiteral("label_regmask_ImgName"));

        horizontalLayout_3->addWidget(label_regmask_ImgName);

        lineEdit_regMaskImgName = new QLineEdit(layoutWidget);
        lineEdit_regMaskImgName->setObjectName(QStringLiteral("lineEdit_regMaskImgName"));

        horizontalLayout_3->addWidget(lineEdit_regMaskImgName);

        pushButton_regMaskImgSave = new QPushButton(layoutWidget);
        pushButton_regMaskImgSave->setObjectName(QStringLiteral("pushButton_regMaskImgSave"));

        horizontalLayout_3->addWidget(pushButton_regMaskImgSave);

        horizontalSlider_regGridMaskProcImg = new QSlider(RegGridMaskDialog);
        horizontalSlider_regGridMaskProcImg->setObjectName(QStringLiteral("horizontalSlider_regGridMaskProcImg"));
        horizontalSlider_regGridMaskProcImg->setGeometry(QRect(100, 530, 400, 17));
        horizontalSlider_regGridMaskProcImg->setOrientation(Qt::Horizontal);
        horizontalSlider_regGridOrigImg = new QSlider(RegGridMaskDialog);
        horizontalSlider_regGridOrigImg->setObjectName(QStringLiteral("horizontalSlider_regGridOrigImg"));
        horizontalSlider_regGridOrigImg->setGeometry(QRect(550, 890, 370, 17));
        horizontalSlider_regGridOrigImg->setOrientation(Qt::Horizontal);

        retranslateUi(RegGridMaskDialog);

        QMetaObject::connectSlotsByName(RegGridMaskDialog);
    } // setupUi

    void retranslateUi(QDialog *RegGridMaskDialog)
    {
        RegGridMaskDialog->setWindowTitle(QApplication::translate("RegGridMaskDialog", "Regular Grid Mask Selection", Q_NULLPTR));
        label_regGridMaskProcImg->setText(QString());
        groupBox_regMaskInputs->setTitle(QApplication::translate("RegGridMaskDialog", "Inputs", Q_NULLPTR));
        comboBox_regMask_ImgDim->clear();
        comboBox_regMask_ImgDim->insertItems(0, QStringList()
         << QApplication::translate("RegGridMaskDialog", "2D Data", Q_NULLPTR)
         << QApplication::translate("RegGridMaskDialog", "3D Data", Q_NULLPTR)
         << QApplication::translate("RegGridMaskDialog", "4D Data", Q_NULLPTR)
        );
        label_regMask_percentage->setText(QApplication::translate("RegGridMaskDialog", "Ratio", Q_NULLPTR));
        pushButton_regMask_imgUpl->setText(QApplication::translate("RegGridMaskDialog", "Open Image Data", Q_NULLPTR));
        pushButton_regMask_run->setText(QApplication::translate("RegGridMaskDialog", "Run", Q_NULLPTR));
        label_regMaskOrigImg->setText(QString());
        label_regmask_ImgName->setText(QApplication::translate("RegGridMaskDialog", "File Name", Q_NULLPTR));
        pushButton_regMaskImgSave->setText(QApplication::translate("RegGridMaskDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class RegGridMaskDialog: public Ui_RegGridMaskDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_REGGRIDMASKDIALOG_H
