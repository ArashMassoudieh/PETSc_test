TEMPLATE = app
CONFIG += c++17 core gui
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

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
