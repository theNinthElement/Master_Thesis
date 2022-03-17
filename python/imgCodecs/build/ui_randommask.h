/********************************************************************************
** Form generated from reading UI file 'randommask.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_RANDOMMASK_H
#define UI_RANDOMMASK_H

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

class Ui_RandomMask
{
public:
    QGroupBox *groupBox_randMask_inputs;
    QVBoxLayout *verticalLayout_2;
    QVBoxLayout *verticalLayout;
    QComboBox *comboBox_randMask;
    QHBoxLayout *horizontalLayout;
    QLabel *label_randMask_percentage;
    QLineEdit *lineEdit_randMask_percentage;
    QPushButton *pushButton_randMask_imgUpl;
    QSpacerItem *verticalSpacer;
    QFrame *line;
    QPushButton *pushButton_randMask_run;
    QLabel *label_randMaks_ProcImg;
    QLabel *label_randMask_origImg;
    QWidget *layoutWidget;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_randMask_fileName;
    QLineEdit *lineEdit_randMask_fileName;
    QPushButton *pushButton_randMask_fileNameSave;
    QSlider *horizontalSlider_radnMaskProcImg;
    QSlider *horizontalSlider_randMaskOrigImg;

    void setupUi(QDialog *RandomMask)
    {
        if (RandomMask->objectName().isEmpty())
            RandomMask->setObjectName(QStringLiteral("RandomMask"));
        RandomMask->resize(1000, 1000);
        groupBox_randMask_inputs = new QGroupBox(RandomMask);
        groupBox_randMask_inputs->setObjectName(QStringLiteral("groupBox_randMask_inputs"));
        groupBox_randMask_inputs->setGeometry(QRect(670, 90, 220, 220));
        verticalLayout_2 = new QVBoxLayout(groupBox_randMask_inputs);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        comboBox_randMask = new QComboBox(groupBox_randMask_inputs);
        comboBox_randMask->setObjectName(QStringLiteral("comboBox_randMask"));

        verticalLayout->addWidget(comboBox_randMask);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_randMask_percentage = new QLabel(groupBox_randMask_inputs);
        label_randMask_percentage->setObjectName(QStringLiteral("label_randMask_percentage"));

        horizontalLayout->addWidget(label_randMask_percentage);

        lineEdit_randMask_percentage = new QLineEdit(groupBox_randMask_inputs);
        lineEdit_randMask_percentage->setObjectName(QStringLiteral("lineEdit_randMask_percentage"));

        horizontalLayout->addWidget(lineEdit_randMask_percentage);


        verticalLayout->addLayout(horizontalLayout);

        pushButton_randMask_imgUpl = new QPushButton(groupBox_randMask_inputs);
        pushButton_randMask_imgUpl->setObjectName(QStringLiteral("pushButton_randMask_imgUpl"));

        verticalLayout->addWidget(pushButton_randMask_imgUpl);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        line = new QFrame(groupBox_randMask_inputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_randMask_run = new QPushButton(groupBox_randMask_inputs);
        pushButton_randMask_run->setObjectName(QStringLiteral("pushButton_randMask_run"));

        verticalLayout->addWidget(pushButton_randMask_run);


        verticalLayout_2->addLayout(verticalLayout);

        label_randMaks_ProcImg = new QLabel(RandomMask);
        label_randMaks_ProcImg->setObjectName(QStringLiteral("label_randMaks_ProcImg"));
        label_randMaks_ProcImg->setGeometry(QRect(100, 100, 420, 420));
        label_randMask_origImg = new QLabel(RandomMask);
        label_randMask_origImg->setObjectName(QStringLiteral("label_randMask_origImg"));
        label_randMask_origImg->setGeometry(QRect(550, 540, 390, 390));
        layoutWidget = new QWidget(RandomMask);
        layoutWidget->setObjectName(QStringLiteral("layoutWidget"));
        layoutWidget->setGeometry(QRect(90, 850, 401, 31));
        horizontalLayout_2 = new QHBoxLayout(layoutWidget);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        label_randMask_fileName = new QLabel(layoutWidget);
        label_randMask_fileName->setObjectName(QStringLiteral("label_randMask_fileName"));

        horizontalLayout_2->addWidget(label_randMask_fileName);

        lineEdit_randMask_fileName = new QLineEdit(layoutWidget);
        lineEdit_randMask_fileName->setObjectName(QStringLiteral("lineEdit_randMask_fileName"));

        horizontalLayout_2->addWidget(lineEdit_randMask_fileName);

        pushButton_randMask_fileNameSave = new QPushButton(layoutWidget);
        pushButton_randMask_fileNameSave->setObjectName(QStringLiteral("pushButton_randMask_fileNameSave"));

        horizontalLayout_2->addWidget(pushButton_randMask_fileNameSave);

        horizontalSlider_radnMaskProcImg = new QSlider(RandomMask);
        horizontalSlider_radnMaskProcImg->setObjectName(QStringLiteral("horizontalSlider_radnMaskProcImg"));
        horizontalSlider_radnMaskProcImg->setGeometry(QRect(110, 550, 400, 17));
        horizontalSlider_radnMaskProcImg->setOrientation(Qt::Horizontal);
        horizontalSlider_randMaskOrigImg = new QSlider(RandomMask);
        horizontalSlider_randMaskOrigImg->setObjectName(QStringLiteral("horizontalSlider_randMaskOrigImg"));
        horizontalSlider_randMaskOrigImg->setGeometry(QRect(560, 950, 370, 17));
        horizontalSlider_randMaskOrigImg->setOrientation(Qt::Horizontal);

        retranslateUi(RandomMask);

        QMetaObject::connectSlotsByName(RandomMask);
    } // setupUi

    void retranslateUi(QDialog *RandomMask)
    {
        RandomMask->setWindowTitle(QApplication::translate("RandomMask", "Randomly Data Selection", Q_NULLPTR));
        groupBox_randMask_inputs->setTitle(QApplication::translate("RandomMask", "Inputs", Q_NULLPTR));
        comboBox_randMask->clear();
        comboBox_randMask->insertItems(0, QStringList()
         << QApplication::translate("RandomMask", "2D Image", Q_NULLPTR)
         << QApplication::translate("RandomMask", "3D Image", Q_NULLPTR)
         << QApplication::translate("RandomMask", "4D Image", Q_NULLPTR)
        );
        label_randMask_percentage->setText(QApplication::translate("RandomMask", "Percentage", Q_NULLPTR));
        pushButton_randMask_imgUpl->setText(QApplication::translate("RandomMask", "Open Image Data", Q_NULLPTR));
        pushButton_randMask_run->setText(QApplication::translate("RandomMask", "Run", Q_NULLPTR));
        label_randMaks_ProcImg->setText(QString());
        label_randMask_origImg->setText(QString());
        label_randMask_fileName->setText(QApplication::translate("RandomMask", "File Name", Q_NULLPTR));
        pushButton_randMask_fileNameSave->setText(QApplication::translate("RandomMask", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class RandomMask: public Ui_RandomMask {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RANDOMMASK_H
