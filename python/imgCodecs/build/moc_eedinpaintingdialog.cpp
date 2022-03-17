/****************************************************************************
** Meta object code from reading C++ file 'eedinpaintingdialog.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../eedinpaintingdialog.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'eedinpaintingdialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_EEDInpaintingDialog_t {
    QByteArrayData data[9];
    char stringdata0[296];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_EEDInpaintingDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_EEDInpaintingDialog_t qt_meta_stringdata_EEDInpaintingDialog = {
    {
QT_MOC_LITERAL(0, 0, 19), // "EEDInpaintingDialog"
QT_MOC_LITERAL(1, 20, 45), // "on_radioButton_eedInpaintExpl..."
QT_MOC_LITERAL(2, 66, 0), // ""
QT_MOC_LITERAL(3, 67, 36), // "on_radioButton_eedInpaintFED_..."
QT_MOC_LITERAL(4, 104, 36), // "on_radioButton_eedInpaintFSI_..."
QT_MOC_LITERAL(5, 141, 38), // "on_pushButton_eedInpaintUplIm..."
QT_MOC_LITERAL(6, 180, 39), // "on_pushButton_eedInpaintUplMa..."
QT_MOC_LITERAL(7, 220, 35), // "on_pushButton_eedInpaintRun_c..."
QT_MOC_LITERAL(8, 256, 39) // "on_pushButton_eedInpaintSaveI..."

    },
    "EEDInpaintingDialog\0"
    "on_radioButton_eedInpaintExplitScheme_clicked\0"
    "\0on_radioButton_eedInpaintFED_clicked\0"
    "on_radioButton_eedInpaintFSI_clicked\0"
    "on_pushButton_eedInpaintUplImg_clicked\0"
    "on_pushButton_eedInpaintUplMask_clicked\0"
    "on_pushButton_eedInpaintRun_clicked\0"
    "on_pushButton_eedInpaintSaveImg_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_EEDInpaintingDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x08 /* Private */,
       3,    0,   50,    2, 0x08 /* Private */,
       4,    0,   51,    2, 0x08 /* Private */,
       5,    0,   52,    2, 0x08 /* Private */,
       6,    0,   53,    2, 0x08 /* Private */,
       7,    0,   54,    2, 0x08 /* Private */,
       8,    0,   55,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void EEDInpaintingDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        EEDInpaintingDialog *_t = static_cast<EEDInpaintingDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->on_radioButton_eedInpaintExplitScheme_clicked(); break;
        case 1: _t->on_radioButton_eedInpaintFED_clicked(); break;
        case 2: _t->on_radioButton_eedInpaintFSI_clicked(); break;
        case 3: _t->on_pushButton_eedInpaintUplImg_clicked(); break;
        case 4: _t->on_pushButton_eedInpaintUplMask_clicked(); break;
        case 5: _t->on_pushButton_eedInpaintRun_clicked(); break;
        case 6: _t->on_pushButton_eedInpaintSaveImg_clicked(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject EEDInpaintingDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_EEDInpaintingDialog.data,
      qt_meta_data_EEDInpaintingDialog,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *EEDInpaintingDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *EEDInpaintingDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_EEDInpaintingDialog.stringdata0))
        return static_cast<void*>(this);
    return QDialog::qt_metacast(_clname);
}

int EEDInpaintingDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
