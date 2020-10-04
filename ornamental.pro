TEMPLATE = app
CONFIG += console c++1z
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXX = /usr/bin/g++-7

USER_CXXFLAGS = -std=c++17
QMAKE_CXXFLAGS_DEBUG += $$USER_CXXFLAGS
QMAKE_CXXFLAGS_RELEASE += $$USER_CXXFLAGS -O3 -g

INCLUDEPATH += ./include
DEPENDPATH += ./include
LIBS += -lgsl -lgslcblas -lcurses
SOURCES += main.cpp
HEADERS += \
    include/container/vector.h \
    include/defs.h \
    include/gnu/colors.h \
    include/gnu/gnuplot.h \
    include/gsl/minimizer.h \
    include/math/WLSM.h \
    include/math/curve.h \
    include/math/distribution.h \
    include/math/math.h \
    include/math/maxmin.h \
    include/math/point.h \
    include/math/random.h \
    include/math/statistics.h \
    include/orn/command.h \
    include/orn/defs.h \
    include/orn/flower.h \
    include/orn/flowerdome.h \
    include/orn/flowerline.h \
    include/orn/gnuplot.h \
    include/orn/orn.defs.h \
    include/orn/phase.h \
    include/orn/recovery.h \
    include/orn/scheme.h \
    include/orn/sprout.h \
    include/orn/sproutmodel.h \
    include/orn/svg.h \
    include/orn/time.h \
    include/os/command.h \
    include/os/os.h \
    include/container/string.h \
    include/svg/defs.h \
    include/svg/primitive.h \
    include/test/logger.h
