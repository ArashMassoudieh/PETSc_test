TEMPLATE = app
CONFIG += c++17 core gui
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    Matrix.cpp \
    Matrix_arma.cpp \
    Vector.cpp \
    Vector_arma.cpp \
    grid.cpp \
    petscmatrix.cpp \
    petscsolver.cpp \
    petscvector.cpp

# Use the MPI wrapper for both compile and link steps
QMAKE_CXX = mpicxx
QMAKE_LINK = mpicxx
QMAKE_LINK_SHLIB = mpicxx

# PETSc paths (from your environment)
PETSC_DIR = $$(PETSC_DIR)
PETSC_ARCH = $$(PETSC_ARCH)

INCLUDEPATH += $$PETSC_DIR/include $$PETSC_DIR/$$PETSC_ARCH/include
LIBS += -L$$PETSC_DIR/$$PETSC_ARCH/lib -lpetsc

# Add MPI explicitly in case the wrapper isn't used at link time
LIBS += -lmpi

# Ensure the runtime can find PETSc (and optionally CUDA) without LD_LIBRARY_PATH
QMAKE_LFLAGS += -Wl,-rpath,$$PETSC_DIR/$$PETSC_ARCH/lib
# If needed for CUDA:
# QMAKE_LFLAGS += -Wl,-rpath,/usr/local/cuda-12.8/lib64

# (Optional) warnings
QMAKE_CXXFLAGS += -Wall -Wextra -Wpedantic

HEADERS += \
    Matrix.h \
    Matrix_arma.h \
    TimeSeries.h \
    TimeSeries.hpp \
    TimeSeriesSet.h \
    TimeSeriesSet.hpp \
    Vector.h \
    Vector_arma.h \
    grid.h \
    petsc_init.h \
    petscmatrix.h \
    petscsolver.h \
    petscvector.h

DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
LIBS += -larmadillo -llapack -lblas -lgsl
