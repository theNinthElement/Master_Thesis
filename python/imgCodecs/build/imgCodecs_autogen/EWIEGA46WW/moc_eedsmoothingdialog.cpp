/****************************************************************************
** Meta object code from reading C++ file 'eedsmoothingdialog.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.7)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../eedsmoothingdialog.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'eedsmoothingdialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.7. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_EEDSmoothingDialog_t {
    QByteArrayData data[7];
    char stringdata0[241];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_EEDSmoothingDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_EEDSmoothingDialog_t qt_meta_stringdata_EEDSmoothingDialog = {
    {
QT_MOC_LITERAL(0, 0, 18), // "EEDSmoothingDialog"
QT_MOC_LITERAL(1, 19, 46), // "on_radioButton_EEDSmoothing_e..."
QT_MOC_LITERAL(2, 66, 0), // ""
QT_MOC_LITERAL(3, 67, 39), // "on_radioButton_EEDSmoothing_F..."
QT_MOC_LITERAL(4, 107, 46), // "on_pushButton_EEDSmoothing_sv..."
QT_MOC_LITERAL(5, 154, 47), // "on_pushButton__EEDSmoothing_L..."
QT_MOC_LITERAL(6, 202, 38) // "on_pushButton_EEDSmoothing_ru..."

    },
    "EEDSmoothingDialog\0"
    "on_radioButton_EEDSmoothing_explScheme_clicked\0"
    "\0on_radioButton_EEDSmoothing_FED_clicked\0"
    "on_pushButton_EEDSmoothing_svaeImgName_clicked\0"
    "on_pushButton__EEDSmoothing_LoadImgData_clicked\0"
    "on_pushButton_EEDSmoothing_run_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_EEDSmoothingDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   39,    2, 0x08 /* Private */,
       3,    0,   40,    2, 0x08 /* Private */,
       4,    0,   41,    2, 0x08 /* Private */,
       5,    0,   42,    2, 0x08 /* Private */,
       6,    0,   43,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void EEDSmoothingDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        EEDSmoothingDialog *_t = static_cast<EEDSmoothingDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->on_radioButton_EEDSmoothing_explScheme_clicked(); break;
        case 1: _t->on_radioButton_EEDSmoothing_FED_clicked(); break;
        case 2: _t->on_pushButton_EEDSmoothing_svaeImgName_clicked(); break;
        case 3: _t->on_pushButton__EEDSmoothing_LoadImgData_clicked(); break;
        case 4: _t->on_pushButton_EEDSmoothing_run_clicked(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject EEDSmoothingDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_EEDSmoothingDialog.data,
      qt_meta_data_EEDSmoothingDialog,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *EEDSmoothingDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *EEDSmoothingDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_EEDSmoothingDialog.stringdata0))
        return static_cast<void*>(this);
    return QDialog::qt_metacast(_clname);
}

int EEDSmoothingDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 5)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 5;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 5)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 5;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
