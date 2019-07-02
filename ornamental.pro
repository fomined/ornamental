TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
TARGET = ornamental

USER_CXXFLAGS = -std=c++11 -Wall

QMAKE_CXXFLAGS_DEBUG += $$USER_CXXFLAGS -D_DEBUG -D_IGNORE
QMAKE_CXXFLAGS_RELEASE += $$USER_CXXFLAGS -O3 -g 

INCLUDEPATH += ./include ./include/prg
DEPENDPATH += ./include

LIBS += -lrt -lgsl -lgslcblas
LIBS += \
	-lboost_system \
	-lboost_filesystem \
	-lboost_date_time \

SOURCES += main.cpp \
    src/flower.structure.cpp \
	src/test.input.file.cpp \
	src/data.recovery.cpp \
    src/math.WLSM.cpp \
    src/flower.ensemble.cpp \
    src/flower.dome.cpp \
    src/plot.flowering.cpp \
    src/flower.generator.cpp

HEADERS += \
    include/flower.structure.h \
    include/prg/math/prg.distribution.h \
    include/prg/prg.array.h \
    include/prg/prg.curve.h \
    include/prg/prg.defs.h \
    include/prg/prg.iterator.h \
    include/prg/prg.logger.h \
    include/prg/prg.math.h \
    include/prg/prg.maxmin.h \
    include/prg/prg.os.h \
    include/prg/prg.point.h \
    include/prg/prg.random.h \
    include/prg/prg.string.h \
	include/readme.h \
	include/defs.h \
	include/math.defs.h \
	include/math.linalg.h \
	include/math.statistics.h \
	include/record.flower.h \
	include/record.phase.h \
	include/record.time.h \
	include/data.recovery.h \
	include/test.defs.h \
	include/test.input.file.h \
    include/math.gnuplot.h \
    include/math.WLSM.h \
    include/flower.ensemble.h \
    include/flower.dome.h \
    include/readme.experiments.h \
    include/plot.flowering.h \
    include/flower.generator.h \
    include/math.gnuplot.h
