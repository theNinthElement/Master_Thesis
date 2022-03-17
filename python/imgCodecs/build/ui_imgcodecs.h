/********************************************************************************
** Form generated from reading UI file 'imgcodecs.ui'
**
** Created by: Qt User Interface Compiler version 5.9.7
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_IMGCODECS_H
#define UI_IMGCODECS_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QFrame>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTabWidget>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ImgCodecs
{
public:
    QWidget *centralWidget;
    QHBoxLayout *horizontalLayout;
    QVBoxLayout *verticalLayout;
    QTabWidget *tabWidget_codec;
    QWidget *tab;
    QHBoxLayout *horizontalLayout_2;
    QVBoxLayout *verticalLayout_2;
    QSpacerItem *verticalSpacer;
    QPushButton *pushButton_RandMask;
    QSpacerItem *verticalSpacer_3;
    QPushButton *pushButton_RegGridMask;
    QSpacerItem *verticalSpacer_4;
    QPushButton *pushButton_RectSubMask;
    QSpacerItem *verticalSpacer_2;
    QPushButton *pushButton_chainCode;
    QWidget *tab_2;
    QHBoxLayout *horizontalLayout_3;
    QVBoxLayout *verticalLayout_3;
    QSpacerItem *verticalSpacer_5;
    QPushButton *pushButton;
    QPushButton *pushButton_eed_4_4D;
    QSpacerItem *verticalSpacer_6;
    QPushButton *pushButton_2dEED_inpt;
    QSpacerItem *verticalSpacer_11;
    QPushButton *pushButton_2;
    QSpacerItem *verticalSpacer_7;
    QWidget *tab_3;
    QHBoxLayout *horizontalLayout_4;
    QVBoxLayout *verticalLayout_4;
    QSpacerItem *verticalSpacer_8;
    QPushButton *pushButton_diffRes;
    QSpacerItem *verticalSpacer_10;
    QPushButton *pushButton_xorRes;
    QSpacerItem *verticalSpacer_9;
    QWidget *tab_4;
    QVBoxLayout *verticalLayout_6;
    QVBoxLayout *verticalLayout_5;
    QPushButton *pushButton__EEDSmoothing;
    QPushButton *pushButton_3;
    QPushButton *pushButton_ZaxisCombine;
    QWidget *tab_5;
    QPushButton *pushButton_4;
    QPushButton *pushButton_5;
    QWidget *tab_6;
    QPushButton *pushButton_codecPDE_LinHom;
    QFrame *line;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ImgCodecs)
    {
        if (ImgCodecs->objectName().isEmpty())
            ImgCodecs->setObjectName(QStringLiteral("ImgCodecs"));
        ImgCodecs->resize(482, 426);
        centralWidget = new QWidget(ImgCodecs);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        horizontalLayout = new QHBoxLayout(centralWidget);
        horizontalLayout->setSpacing(6);
        horizontalLayout->setContentsMargins(11, 11, 11, 11);
        horizontalLayout->setObjectName(QStringLiteral("horizontalLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QStringLiteral("verticalLayout"));
        tabWidget_codec = new QTabWidget(centralWidget);
        tabWidget_codec->setObjectName(QStringLiteral("tabWidget_codec"));
        tab = new QWidget();
        tab->setObjectName(QStringLiteral("tab"));
        horizontalLayout_2 = new QHBoxLayout(tab);
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_2->setObjectName(QStringLiteral("horizontalLayout_2"));
        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QStringLiteral("verticalLayout_2"));
        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer);

        pushButton_RandMask = new QPushButton(tab);
        pushButton_RandMask->setObjectName(QStringLiteral("pushButton_RandMask"));

        verticalLayout_2->addWidget(pushButton_RandMask);

        verticalSpacer_3 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer_3);

        pushButton_RegGridMask = new QPushButton(tab);
        pushButton_RegGridMask->setObjectName(QStringLiteral("pushButton_RegGridMask"));

        verticalLayout_2->addWidget(pushButton_RegGridMask);

        verticalSpacer_4 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer_4);

        pushButton_RectSubMask = new QPushButton(tab);
        pushButton_RectSubMask->setObjectName(QStringLiteral("pushButton_RectSubMask"));

        verticalLayout_2->addWidget(pushButton_RectSubMask);

        verticalSpacer_2 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer_2);

        pushButton_chainCode = new QPushButton(tab);
        pushButton_chainCode->setObjectName(QStringLiteral("pushButton_chainCode"));

        verticalLayout_2->addWidget(pushButton_chainCode);


        horizontalLayout_2->addLayout(verticalLayout_2);

        tabWidget_codec->addTab(tab, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QStringLiteral("tab_2"));
        horizontalLayout_3 = new QHBoxLayout(tab_2);
        horizontalLayout_3->setSpacing(6);
        horizontalLayout_3->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_3->setObjectName(QStringLiteral("horizontalLayout_3"));
        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QStringLiteral("verticalLayout_3"));
        verticalSpacer_5 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer_5);

        pushButton = new QPushButton(tab_2);
        pushButton->setObjectName(QStringLiteral("pushButton"));

        verticalLayout_3->addWidget(pushButton);

        pushButton_eed_4_4D = new QPushButton(tab_2);
        pushButton_eed_4_4D->setObjectName(QStringLiteral("pushButton_eed_4_4D"));

        verticalLayout_3->addWidget(pushButton_eed_4_4D);

        verticalSpacer_6 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer_6);

        pushButton_2dEED_inpt = new QPushButton(tab_2);
        pushButton_2dEED_inpt->setObjectName(QStringLiteral("pushButton_2dEED_inpt"));

        verticalLayout_3->addWidget(pushButton_2dEED_inpt);

        verticalSpacer_11 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer_11);

        pushButton_2 = new QPushButton(tab_2);
        pushButton_2->setObjectName(QStringLiteral("pushButton_2"));

        verticalLayout_3->addWidget(pushButton_2);

        verticalSpacer_7 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer_7);


        horizontalLayout_3->addLayout(verticalLayout_3);

        tabWidget_codec->addTab(tab_2, QString());
        tab_3 = new QWidget();
        tab_3->setObjectName(QStringLiteral("tab_3"));
        horizontalLayout_4 = new QHBoxLayout(tab_3);
        horizontalLayout_4->setSpacing(6);
        horizontalLayout_4->setContentsMargins(11, 11, 11, 11);
        horizontalLayout_4->setObjectName(QStringLiteral("horizontalLayout_4"));
        verticalLayout_4 = new QVBoxLayout();
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setObjectName(QStringLiteral("verticalLayout_4"));
        verticalSpacer_8 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_4->addItem(verticalSpacer_8);

        pushButton_diffRes = new QPushButton(tab_3);
        pushButton_diffRes->setObjectName(QStringLiteral("pushButton_diffRes"));

        verticalLayout_4->addWidget(pushButton_diffRes);

        verticalSpacer_10 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_4->addItem(verticalSpacer_10);

        pushButton_xorRes = new QPushButton(tab_3);
        pushButton_xorRes->setObjectName(QStringLiteral("pushButton_xorRes"));

        verticalLayout_4->addWidget(pushButton_xorRes);

        verticalSpacer_9 = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_4->addItem(verticalSpacer_9);


        horizontalLayout_4->addLayout(verticalLayout_4);

        tabWidget_codec->addTab(tab_3, QString());
        tab_4 = new QWidget();
        tab_4->setObjectName(QStringLiteral("tab_4"));
        verticalLayout_6 = new QVBoxLayout(tab_4);
        verticalLayout_6->setSpacing(6);
        verticalLayout_6->setContentsMargins(11, 11, 11, 11);
        verticalLayout_6->setObjectName(QStringLiteral("verticalLayout_6"));
        verticalLayout_5 = new QVBoxLayout();
        verticalLayout_5->setSpacing(6);
        verticalLayout_5->setObjectName(QStringLiteral("verticalLayout_5"));
        pushButton__EEDSmoothing = new QPushButton(tab_4);
        pushButton__EEDSmoothing->setObjectName(QStringLiteral("pushButton__EEDSmoothing"));

        verticalLayout_5->addWidget(pushButton__EEDSmoothing);

        pushButton_3 = new QPushButton(tab_4);
        pushButton_3->setObjectName(QStringLiteral("pushButton_3"));

        verticalLayout_5->addWidget(pushButton_3);

        pushButton_ZaxisCombine = new QPushButton(tab_4);
        pushButton_ZaxisCombine->setObjectName(QStringLiteral("pushButton_ZaxisCombine"));

        verticalLayout_5->addWidget(pushButton_ZaxisCombine);


        verticalLayout_6->addLayout(verticalLayout_5);

        tabWidget_codec->addTab(tab_4, QString());
        tab_5 = new QWidget();
        tab_5->setObjectName(QStringLiteral("tab_5"));
        pushButton_4 = new QPushButton(tab_5);
        pushButton_4->setObjectName(QStringLiteral("pushButton_4"));
        pushButton_4->setGeometry(QRect(140, 80, 87, 29));
        pushButton_5 = new QPushButton(tab_5);
        pushButton_5->setObjectName(QStringLiteral("pushButton_5"));
        pushButton_5->setGeometry(QRect(140, 130, 87, 29));
        tabWidget_codec->addTab(tab_5, QString());
        tab_6 = new QWidget();
        tab_6->setObjectName(QStringLiteral("tab_6"));
        pushButton_codecPDE_LinHom = new QPushButton(tab_6);
        pushButton_codecPDE_LinHom->setObjectName(QStringLiteral("pushButton_codecPDE_LinHom"));
        pushButton_codecPDE_LinHom->setGeometry(QRect(170, 50, 91, 29));
        tabWidget_codec->addTab(tab_6, QString());

        verticalLayout->addWidget(tabWidget_codec);

        line = new QFrame(centralWidget);
        line->setObjectName(QStringLiteral("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayout->addWidget(line);


        horizontalLayout->addLayout(verticalLayout);

        ImgCodecs->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ImgCodecs);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 482, 23));
        ImgCodecs->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ImgCodecs);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        ImgCodecs->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ImgCodecs);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        ImgCodecs->setStatusBar(statusBar);

        retranslateUi(ImgCodecs);

        tabWidget_codec->setCurrentIndex(5);


        QMetaObject::connectSlotsByName(ImgCodecs);
    } // setupUi

    void retranslateUi(QMainWindow *ImgCodecs)
    {
        ImgCodecs->setWindowTitle(QApplication::translate("ImgCodecs", "Medical Image Codecs", Q_NULLPTR));
        pushButton_RandMask->setText(QApplication::translate("ImgCodecs", "Random", Q_NULLPTR));
        pushButton_RegGridMask->setText(QApplication::translate("ImgCodecs", "Regular Grid", Q_NULLPTR));
        pushButton_RectSubMask->setText(QApplication::translate("ImgCodecs", "Rectangular Subdivision", Q_NULLPTR));
        pushButton_chainCode->setText(QApplication::translate("ImgCodecs", "Chain Code", Q_NULLPTR));
        tabWidget_codec->setTabText(tabWidget_codec->indexOf(tab), QApplication::translate("ImgCodecs", "Encoding", Q_NULLPTR));
        pushButton->setText(QApplication::translate("ImgCodecs", "Edge Enhancing Diffusion", Q_NULLPTR));
        pushButton_eed_4_4D->setText(QApplication::translate("ImgCodecs", "EED for 4D", Q_NULLPTR));
        pushButton_2dEED_inpt->setText(QApplication::translate("ImgCodecs", "Slice-wise EED", Q_NULLPTR));
        pushButton_2->setText(QApplication::translate("ImgCodecs", "Fourth-Order Edge Enhancing Diffusion", Q_NULLPTR));
        tabWidget_codec->setTabText(tabWidget_codec->indexOf(tab_2), QApplication::translate("ImgCodecs", "Decoding", Q_NULLPTR));
        pushButton_diffRes->setText(QApplication::translate("ImgCodecs", "Difference", Q_NULLPTR));
        pushButton_xorRes->setText(QApplication::translate("ImgCodecs", "XOR", Q_NULLPTR));
        tabWidget_codec->setTabText(tabWidget_codec->indexOf(tab_3), QApplication::translate("ImgCodecs", "Residuals", Q_NULLPTR));
        pushButton__EEDSmoothing->setText(QApplication::translate("ImgCodecs", "EED Smoothing", Q_NULLPTR));
        pushButton_3->setText(QApplication::translate("ImgCodecs", "FOEED Smoothing", Q_NULLPTR));
        pushButton_ZaxisCombine->setText(QApplication::translate("ImgCodecs", "Time-axis Combine", Q_NULLPTR));
        tabWidget_codec->setTabText(tabWidget_codec->indexOf(tab_4), QApplication::translate("ImgCodecs", "Pre-processing", Q_NULLPTR));
        pushButton_4->setText(QApplication::translate("ImgCodecs", "Dilation", Q_NULLPTR));
        pushButton_5->setText(QApplication::translate("ImgCodecs", "Primize", Q_NULLPTR));
        tabWidget_codec->setTabText(tabWidget_codec->indexOf(tab_5), QApplication::translate("ImgCodecs", "Tools", Q_NULLPTR));
        pushButton_codecPDE_LinHom->setText(QApplication::translate("ImgCodecs", "R-ILH-0", Q_NULLPTR));
        tabWidget_codec->setTabText(tabWidget_codec->indexOf(tab_6), QApplication::translate("ImgCodecs", "Decompression", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class ImgCodecs: public Ui_ImgCodecs {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_IMGCODECS_H
