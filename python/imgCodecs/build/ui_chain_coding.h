/********************************************************************************
** Form generated from reading UI file 'chain_coding.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_CHAIN_CODING_H
#define UI_CHAIN_CODING_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_Chain_coding
{
public:
    QGroupBox *chainCode_inputs;
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label_chainCode_cannyLowerThres;
    QLineEdit *lineEdit_chainCode_cannyLowerThres;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_chainCode_upperThres;
    QLineEdit *lineEdit_chainCode_upperThres;
    QPushButton *pushButton_chainCode_loadData;
    QFrame *line;
    QPushButton *pushButton_chainCode_run;
    QLabel *label_chainCode_data;
    QWidget *widget;
    QHBoxLayout *horizontalLayout_3;
    QLabel *label_chainCode_Output;
    QLineEdit *lineEdit_chainCode_Output;
    QPushButton *pushButton_chainCode_Output;

    void setupUi(QDialog *Chain_coding)
    {
        if (Chain_coding->objectName().isEmpty())
            Chain_coding->setObjectName(QStringLiteral("Chain_coding"));
        Chain_coding->resize(1162, 969);
        chainCode_inputs = new QGroupBox(Chain_coding);
        chainCode_inputs->setObjectName(QStringLiteral("chainCode_inputs"));
        chainCode_inputs->setGeometry(QRect(770, 60, 321, 221));
        gridLayout = new QGridLayout(chainCode_inputs);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        label_chainCode_cannyLowerThres = new QLabel(chainCode_inputs);
        label_chainCode_cannyLowerThres->setObjectName(QStringLiteral("label_chainCode_cannyLowerThres"));

        horizontalLayout->addWidget(label_chainCode_cannyLowerThres);

        lineEdit_chainCode_cannyLowerThres = new QLineEdit(chainCode_inputs);
        lineEdit_chainCode_cannyLowerThres->setObjectName(QStringLiteral("lineEdit_chainCode_cannyLowerThres"));

        horizontalLayout->addWidget(lineEdit_chainCode_cannyLowerThres);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        label_chainCode_upperThres = new QLabel(chainCode_inputs);
        label_chainCode_upperThres->setObjectName(QStringLiteral("label_chainCode_upperThres"));

        horizontalLayout_2->addWidget(label_chainCode_upperThres);

        lineEdit_chainCode_upperThres = new QLineEdit(chainCode_inputs);
        lineEdit_chainCode_upperThres->setObjectName(QStringLiteral("lineEdit_chainCode_upperThres"));

        horizontalLayout_2->addWidget(lineEdit_chainCode_upperThres);


        verticalLayout->addLayout(horizontalLayout_2);

        pushButton_chainCode_loadData = new QPushButton(chainCode_inputs);
        pushButton_chainCode_loadData->setObjectName(QStringLiteral("pushButton_chainCode_loadData"));

        verticalLayout->addWidget(pushButton_chainCode_loadData);

        line = new QFrame(chainCode_inputs);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);

        pushButton_chainCode_run = new QPushButton(chainCode_inputs);
        pushButton_chainCode_run->setObjectName(QStringLiteral("pushButton_chainCode_run"));

        verticalLayout->addWidget(pushButton_chainCode_run);


        gridLayout->addLayout(verticalLayout, 0, 0, 1, 1);

        label_chainCode_data = new QLabel(Chain_coding);
        label_chainCode_data->setObjectName(QStringLiteral("label_chainCode_data"));
        label_chainCode_data->setGeometry(QRect(60, 50, 611, 541));
        widget = new QWidget(Chain_coding);
        widget->setObjectName(QStringLiteral("widget"));
        widget->setGeometry(QRect(110, 780, 461, 31));
        horizontalLayout_3 = new QHBoxLayout(widget);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        label_chainCode_Output = new QLabel(widget);
        label_chainCode_Output->setObjectName(QStringLiteral("label_chainCode_Output"));

        horizontalLayout_3->addWidget(label_chainCode_Output);

        lineEdit_chainCode_Output = new QLineEdit(widget);
        lineEdit_chainCode_Output->setObjectName(QStringLiteral("lineEdit_chainCode_Output"));

        horizontalLayout_3->addWidget(lineEdit_chainCode_Output);

        pushButton_chainCode_Output = new QPushButton(widget);
        pushButton_chainCode_Output->setObjectName(QStringLiteral("pushButton_chainCode_Output"));

        horizontalLayout_3->addWidget(pushButton_chainCode_Output);


        retranslateUi(Chain_coding);

        QMetaObject::connectSlotsByName(Chain_coding);
    } // setupUi

    void retranslateUi(QDialog *Chain_coding)
    {
        Chain_coding->setWindowTitle(QApplication::translate("Chain_coding", "Chain Coding", Q_NULLPTR));
        chainCode_inputs->setTitle(QApplication::translate("Chain_coding", "Inputs", Q_NULLPTR));
        label_chainCode_cannyLowerThres->setText(QApplication::translate("Chain_coding", "Lower Threshold", Q_NULLPTR));
        label_chainCode_upperThres->setText(QApplication::translate("Chain_coding", "Upper Threshold", Q_NULLPTR));
        pushButton_chainCode_loadData->setText(QApplication::translate("Chain_coding", "Load Data", Q_NULLPTR));
        pushButton_chainCode_run->setText(QApplication::translate("Chain_coding", "Run", Q_NULLPTR));
        label_chainCode_data->setText(QString());
        label_chainCode_Output->setText(QApplication::translate("Chain_coding", "Name", Q_NULLPTR));
        pushButton_chainCode_Output->setText(QApplication::translate("Chain_coding", "Save", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class Chain_coding: public Ui_Chain_coding {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_CHAIN_CODING_H
