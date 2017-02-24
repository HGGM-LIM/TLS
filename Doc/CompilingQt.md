# Compiling Qt {#CompilingQt}
@brief Instructions for building Qt from source.

## Ubuntu

~~~{.sh}
QTDIR="/path/to/Qt"
QT_SOURCE_DIR="/path/to/source"

cd ${QT_SOURCE_DIR}
./configure -prefix ${QTDIR} -confirm-license -fast -no-qt3support -nomake demos -nomake docs -nomake examples -opensource -release -webkit
make
make install
~~~

## OS X

~~~{.sh}
QTDIR="/path/to/Qt"
QT_SOURCE_DIR="/path/to/source"

cd ${QT_SOURCE_DIR}
sdk=`xcodebuild -version \`xcodebuild -showsdks | awk '/^$/{p=0};p; /OS X SDKs:/{p=1}' | tail -1 | cut -f3\` Path`
./configure -prefix ${QTDIR} -arch x86_64 -sdk ${sdk} -confirm-license -fast -no-qt3support -nomake demos -nomake docs -nomake examples -opensource -release -webkit
make
make install
~~~

## Windows

If using Visual Studio, it is recommended to use [jom](http://qt-project.org/wiki/jom) instead of nmake.

From Visual Studio x64 Command Prompt:

~~~{.bat}
set QTDIR="C:\path\to\Qt"

cd %QTDIR%
configure -confirm-license -fast -no-qt3support -nomake demos -nomake docs -nomake examples -opensource -release -webkit
nmake
nmake clean
~~~
