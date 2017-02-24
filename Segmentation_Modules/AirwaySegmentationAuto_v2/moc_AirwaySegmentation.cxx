/****************************************************************************
** AirwaySegmentation meta object code from reading C++ file 'AirwaySegmentation.h'
**
** Created: Thu Jun 26 10:57:06 2008
**      by: The Qt MOC ($Id: qt/moc_yacc.cpp   3.3.3   edited Aug 5 16:40 $)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#undef QT_NO_COMPAT
#include "../../../mitk/QFunctionalities/AirwaySegmentation/AirwaySegmentation.h"
#include <qmetaobject.h>
#include <qapplication.h>

#include <private/qucomextra_p.h>
#if !defined(Q_MOC_OUTPUT_REVISION) || (Q_MOC_OUTPUT_REVISION != 26)
#error "This file was generated using the moc from 3.3.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

const char *AirwaySegmentation::className() const
{
    return "AirwaySegmentation";
}

QMetaObject *AirwaySegmentation::metaObj = 0;
static QMetaObjectCleanUp cleanUp_AirwaySegmentation( "AirwaySegmentation", &AirwaySegmentation::staticMetaObject );

#ifndef QT_NO_TRANSLATION
QString AirwaySegmentation::tr( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "AirwaySegmentation", s, c, QApplication::DefaultCodec );
    else
	return QString::fromLatin1( s );
}
#ifndef QT_NO_TRANSLATION_UTF8
QString AirwaySegmentation::trUtf8( const char *s, const char *c )
{
    if ( qApp )
	return qApp->translate( "AirwaySegmentation", s, c, QApplication::UnicodeUTF8 );
    else
	return QString::fromUtf8( s );
}
#endif // QT_NO_TRANSLATION_UTF8

#endif // QT_NO_TRANSLATION

QMetaObject* AirwaySegmentation::staticMetaObject()
{
    if ( metaObj )
	return metaObj;
    QMetaObject* parentObject = QmitkFunctionality::staticMetaObject();
    static const QUMethod slot_0 = {"TreeChanged", 0, 0 };
    static const QUParameter param_slot_1[] = {
	{ "imageIt", &static_QUType_ptr, "mitk::DataTreeIteratorClone", QUParameter::In }
    };
    static const QUMethod slot_1 = {"ImageSelected", 1, param_slot_1 };
    static const QUMethod slot_2 = {"StartButtonClicked", 0, 0 };
    static const QMetaData slot_tbl[] = {
	{ "TreeChanged()", &slot_0, QMetaData::Protected },
	{ "ImageSelected(mitk::DataTreeIteratorClone)", &slot_1, QMetaData::Protected },
	{ "StartButtonClicked()", &slot_2, QMetaData::Protected }
    };
    metaObj = QMetaObject::new_metaobject(
	"AirwaySegmentation", parentObject,
	slot_tbl, 3,
	0, 0,
#ifndef QT_NO_PROPERTIES
	0, 0,
	0, 0,
#endif // QT_NO_PROPERTIES
	0, 0 );
    cleanUp_AirwaySegmentation.setMetaObject( metaObj );
    return metaObj;
}

void* AirwaySegmentation::qt_cast( const char* clname )
{
    if ( !qstrcmp( clname, "AirwaySegmentation" ) )
	return this;
    return QmitkFunctionality::qt_cast( clname );
}

bool AirwaySegmentation::qt_invoke( int _id, QUObject* _o )
{
    switch ( _id - staticMetaObject()->slotOffset() ) {
    case 0: TreeChanged(); break;
    case 1: ImageSelected((mitk::DataTreeIteratorClone)(*((mitk::DataTreeIteratorClone*)static_QUType_ptr.get(_o+1)))); break;
    case 2: StartButtonClicked(); break;
    default:
	return QmitkFunctionality::qt_invoke( _id, _o );
    }
    return TRUE;
}

bool AirwaySegmentation::qt_emit( int _id, QUObject* _o )
{
    return QmitkFunctionality::qt_emit(_id,_o);
}
#ifndef QT_NO_PROPERTIES

bool AirwaySegmentation::qt_property( int id, int f, QVariant* v)
{
    return QmitkFunctionality::qt_property( id, f, v);
}

bool AirwaySegmentation::qt_static_property( QObject* , int , int , QVariant* ){ return FALSE; }
#endif // QT_NO_PROPERTIES
