/********************************************************************************
** Form generated from reading UI file 'directneighbordilationdialog.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_DIRECTNEIGHBORDILATIONDIALOG_H
#define UI_DIRECTNEIGHBORDILATIONDIALOG_H

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

class Ui_DirectNeighborDilationDialog
{
public:
    QLabel *label_directDilationOrigImg;
    QLabel *label_directDilationMaskImg;
    QLabel *label_directDilationRecImg;
    QGroupBox *groupBox_regMaskInputs;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout;
    QComboBox *comboBox_directDilation_ImgDim;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButton_directDilation_imgUpl;
    QPushButton *pushButton_directDilation_mskUpl;
    QPushButton *pushButton_directDilation_recUpl;
    QSpacerItem *verticalSpacer_2;
    QFrame *line;
    QPushButton *pushButton_directDilation_run;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_directDilation_ImgName;
    QLineEdit *lineEdit_directDilationImgName;
    QPushButton *pushButton_directDilationImgSave;

    void setupUi(QDialog *DirectNeighborDilationDialog)
    {
        if (DirectNeighborDilationDialog->objectName().isEmpty())
            DirectNeighborDilationDialog->setObjectName(QStringLiteral("DirectNeighborDilationDialog"));
        DirectNeighborDilationDialog->resize(1567, 997);
        label_directDilationOrigImg = new QLabel(DirectNeighborDilationDialog);
        label_directDilationOrigImg->setObjectName(QStringLiteral("label_directDilationOrigImg"));
        label_directDilationOrigImg->setGeometry(QRect(20, 30, 531, 431));
        label_directDilationMaskImg = new QLabel(DirectNeighborDilationDialog);
        label_directDilationMaskImg->setObjectName(QStringLiteral("label_directDilationMaskImg"));
        label_directDilationMaskImg->setGeometry(QRect(1070, 550, 421, 351));
        label_directDilationRecImg = new QLabel(DirectNeighborDilationDialog);
        label_directDilationRecImg->setObjectName(QStringLiteral("label_directDilationRecImg"));
        label_directDilationRecImg->setGeometry(QRect(630, 30, 411, 381));
        groupBox_regMaskInputs = new QGroupBox(DirectNeighborDilationDialog);
        groupBox_regMaskInputs->setObjectName(QStringLiteral("groupBox_regMaskInputs"));
        groupBox_regMaskInputs->setGeometry(QRect(1210, 40, 220, 220));
        horizontalLayout_2 = new QHBoxLayout(groupBox_regMaskInputs);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        comboBox_directDilation_ImgDim = new QComboBox(groupBox_regMaskInputs);
        comboBox_directDilation_ImgDim->setObjectName(QStringLiteral("comboBox_directDilation_ImgDim"));

        verticalLayout->addWidget(comboBox_directDilation_ImgDim);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));

        verticalLayout->addLayout(horizontalLayout);

        pushButton_directDilation_imgUpl = new QPushButton(groupBox_regMaskInputs);
        pushButton_directDilation_imgUpl->setObjectName(QStringLiteral("pushButton_directDilation_imgUpl"));

        verticalLayout->addWidget(pushButton_directDilation_imgUpl);

        pushButton_directDilation_mskUpl = new QPushButton(groupBox_regMaskInputs);
        pushButton_directDilation_mskUpl->setObjectName(QStringLiteral("pushButton_directDilation_mskUpl"));

        verticalLayout->addWidget(pushButton_directDilation_mskUpl);

        pushButton_directDilation_recUpl = new QPushButton(groupBox_regMaskInputs);
        pushButton_directDilation_recUpl->setObjectName(QStringLiteral("pushButton_directDilation_recUpl"));

        verticalLayout->addWidget(pushButton_directDilation_recUpl);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer_2);

        line = new QFrame(groupBox_regMaskInputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_directDilation_run = new QPushButton(groupBox_regMaskInputs);
        pushButton_directDilation_run->setObjectName(QStringLiteral("pushButton_directDilation_run"));

        verticalLayout->addWidget(pushButton_directDilation_run);


        horizontalLayout_2->addLayout(verticalLayout);

        layoutWidget = new QWidget(DirectNeighborDilationDialog);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(140, 840, 401, 31));
        horizontalLayout_3 = new QHBoxLayout(layoutWidget);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_directDilation_ImgName = new QLabel(layoutWidget);
        label_directDilation_ImgName->setObjectName(QStringLiteral("label_directDilation_ImgName"));

        horizontalLayout_3->addWidget(label_directDilation_ImgName);

        lineEdit_directDilationImgName = new QLineEdit(layoutWidget);
        lineEdit_directDilationImgName->setObjectName(QStringLiteral("lineEdit_directDilationImgName"));

        horizontalLayout_3->addWidget(lineEdit_directDilationImgName);

        pushButton_directDilationImgSave = new QPushButton(layoutWidget);
        pushButton_directDilationImgSave->setObjectName(QStringLiteral("pushButton_directDilationImgSave"));

        horizontalLayout_3->addWidget(pushButton_directDilationImgSave);


        retranslateUi(DirectNeighborDilationDialog);

        QMetaObject::connectSlotsByName(DirectNeighborDilationDialog);
    } // setupUi

    void retranslateUi(QDialog *DirectNeighborDilationDialog)
    {
        DirectNeighborDilationDialog->setWindowTitle(QApplication::translate("DirectNeighborDilationDialog", "Direct Neighbor Dilation Dialog", Q_NULLPTR));
        label_directDilationOrigImg->setText(QString());
        label_directDilationMaskImg->setText(QString());
        label_directDilationRecImg->setText(QString());
        groupBox_regMaskInputs->setTitle(QApplication::translate("DirectNeighborDilationDialog", "Inputs", Q_NULLPTR));
        comboBox_directDilation_ImgDim->clear();
        comboBox_directDilation_ImgDim->insertItems(0, QStringList()
         << QApplication::translate("DirectNeighborDilationDialog", "2D Data", Q_NULLPTR)
         << QApplication::translate("DirectNeighborDilationDialog", "3D Data", Q_NULLPTR)
         << QApplication::translate("DirectNeighborDilationDialog", "4D Data", Q_NULLPTR)
        );
        pushButton_directDilation_imgUpl->setText(QApplication::translate("DirectNeighborDilationDialog", "Open Image Data", Q_NULLPTR));
        pushButton_directDilation_mskUpl->setText(QApplication::translate("DirectNeighborDilationDialog", "Open Prev. Mask Data", Q_NULLPTR));
        pushButton_directDilation_recUpl->setText(QApplication::translate("DirectNeighborDilationDialog", "Open Reconst. Data", Q_NULLPTR));
        pushButton_directDilation_run->setText(QApplication::translate("DirectNeighborDilationDialog", "Run", Q_NULLPTR));
        label_directDilation_ImgName->setText(QApplication::translate("DirectNeighborDilationDialog", "File Name", Q_NULLPTR));
        pushButton_directDilationImgSave->setText(QApplication::translate("DirectNeighborDilationDialog", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class DirectNeighborDilationDialog: public Ui_DirectNeighborDilationDialog {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_DIRECTNEIGHBORDILATIONDIALOG_H
