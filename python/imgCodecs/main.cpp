#include "imgcodecs.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    ImgCodecs w;
    w.show();

    return a.exec();
}
