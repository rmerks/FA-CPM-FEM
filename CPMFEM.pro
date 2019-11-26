######################################################################
# Automatically generated by qmake (2.01a) ma mei 6 17:10:58 2013
######################################################################

TEMPLATE = app
TARGET = CPMFEM
DEPENDPATH += .
INCLUDEPATH += .
CONFIG -= app_bundle
QT += widgets

# Input
HEADERS += def.h \
           functions.h \
           graph.h \
           output.h \
           myparameters.h \
           parse.h \
           qtgraph.h \
           structures.h \
           warning.h \
    plotcpm.h
SOURCES += cellforces.cpp \
           cellmoves.cpp \
           CPM_dH.cpp \
           FE_assembly.cpp \
           FE_local.cpp \
           FE_nodes2dofs.cpp \
           FE_solver.cpp \
           PDE_solver.cpp \
           init.cpp \
           mt.cpp \
           mylib.cpp \
           output.cpp \
           myparameters.cpp \
           parse.cpp \
           qtgraph.cpp \
           read.cpp \
           warning.cpp \
           write.cpp \
    sandbox.cpp \
    plotcpm.cpp
