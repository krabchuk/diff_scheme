TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += -std=c++11

SOURCES += \
    src/main.cpp \
    src/diff_scheme_solver.cpp

HEADERS += \
    src/cpp_headers.h \
    src/diff_scheme_solver.h

INCLUDEPATH += src
