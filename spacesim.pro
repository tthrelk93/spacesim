QT       += core gui
QT += core gui widgets
QT += core gui widgets openglwidgets
INCLUDEPATH += /usr/local/include
INCLUDEPATH += /opt/homebrew/Cellar/sdl2/2.0.22/include/SDL2/
INCLUDEPATH += /opt/homebrew/Cellar/glew/2.2.0_1/include/
INCLUDEPATH += /opt/homebrew/Cellar/glm/0.9.9.8/include/
LIBS += -L/opt/homebrew/lib -lSDL2


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    simulation.cpp

HEADERS += \
    mainwindow.h \
    simulation.h

FORMS += \
    mainwindow.ui

TRANSLATIONS += \
    spacesim_en_US.ts
CONFIG += lrelease
CONFIG += embed_translations

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
