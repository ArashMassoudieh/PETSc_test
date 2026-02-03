TEMPLATE = app
CONFIG  += c++17 console
QT      -= core gui widgets

DEFINES += GRID_USE_VTK
CONFIG  += GRID_USE_VTK
DEFINES += GSL

# ============================================================
# Host-config (canonical): enable one
# ============================================================
#CONFIG += Behzad
#DEFINES += Behzad

#CONFIG += PowerEdge
#DEFINES += PowerEdge

#CONFIG += Arash
#DEFINES += Arash

#CONFIG += SligoCreek
#DEFINES += SligoCreek

CONFIG  += Jason
DEFINES += Jason

#CONFIG += WSL
#DEFINES += WSL


# ============================================================
# VTK paths
# ============================================================
contains(DEFINES, Behzad) {
    VTKBUILDPATH   = /home/behzad/Projects/VTK-9.3.1/VTK-build
    VTKHEADERPATH  = /home/behzad/Projects/VTK-9.3.1
    VTK_V          = -9.3
}

contains(DEFINES, PowerEdge) {
    VTKBUILDPATH   = ../VTK-build
    VTKHEADERPATH  = ../VTK
    VTK_V          = -9.0
}

contains(DEFINES, Arash) {
    VTKBUILDPATH   = /home/arash/Projects/VTK-9.4.1/build
    VTKHEADERPATH  = /home/arash/Projects/VTK-9.4.1
    VTK_V          = -9.4
}

contains(DEFINES, SligoCreek) {
    VTKBUILDPATH   = /media/arash/E/Projects/VTK-9.1.0/VTK-build
    VTKHEADERPATH  = /media/arash/E/Projects/VTK-9.1.0
    VTK_V          = -9.1
}

# ---- Jason: VTK INSTALL prefix (stable) ----
contains(DEFINES, Jason) {
    VTKINSTALLPATH = /home/arash/Projects/VTK-install
    VTK_INCNAME    = vtk-9.3
    VTK_V          = -9.3

    VTKHEADERPATH  = $$VTKINSTALLPATH/include/$$VTK_INCNAME
    VTKLIBPATH     = $$VTKINSTALLPATH/lib

    !exists($$VTKHEADERPATH/vtkSmartPointer.h) {
        error("VTK not found: $$VTKHEADERPATH/vtkSmartPointer.h")
    }
    !exists($$VTKLIBPATH) {
        error("VTK lib path not found: $$VTKLIBPATH")
    }
}

contains(DEFINES, WSL) {
    VTKINSTALLPATH = /home/behzad/Projects/VTK-install
    VTK_INCNAME    = vtk-9.4
    VTK_V          = -9.4
    VTKHEADERPATH  = $$VTKINSTALLPATH/include/$$VTK_INCNAME
    VTKLIBPATH     = $$VTKINSTALLPATH/lib
}


# ============================================================
# PETSc (universal, qmake-friendly)
#   Order: pkg-config -> petscvariables -> manual -lpetsc
#   Also: auto-detect PETSC_ARCH if wrong/empty.
# ============================================================
PETSC_DIR  = $$(PETSC_DIR)
PETSC_ARCH = $$(PETSC_ARCH)

# Host defaults only if env empty
contains(DEFINES, PowerEdge) {
    isEmpty(PETSC_DIR)  { PETSC_DIR  = /mnt/3rd900/Projects/petsc }
    isEmpty(PETSC_ARCH) { PETSC_ARCH = arch-linux-c-opt }
} else:contains(DEFINES, Jason)|contains(DEFINES, Arash) {
    isEmpty(PETSC_DIR)  { PETSC_DIR  = /home/arash/Projects/petsc }
    isEmpty(PETSC_ARCH) { PETSC_ARCH = arch-linux-c-opt }
}

message("PETSC_DIR  = $$PETSC_DIR")
message("PETSC_ARCH = $$PETSC_ARCH")

# If PETSC_ARCH doesn't exist, auto-pick the first arch-* folder
!exists($$PETSC_DIR/$$PETSC_ARCH) {
    message("NOTE: PETSC_ARCH not found at $$PETSC_DIR/$$PETSC_ARCH. Auto-detecting...")
    PETSC_ARCH = $$system(ls -1 $$PETSC_DIR 2>/dev/null | grep '^arch-' | head -n 1)
    isEmpty(PETSC_ARCH) {
        error("Could not auto-detect PETSC_ARCH under PETSC_DIR=$$PETSC_DIR")
    }
    message("Auto-detected PETSC_ARCH = $$PETSC_ARCH")
}

PETSC_LIBDIR = $$PETSC_DIR/$$PETSC_ARCH/lib
PETSC_PKGCFG = $$PETSC_LIBDIR/pkgconfig
PETSC_VARIABLES = $$PETSC_LIBDIR/petsc/conf/petscvariables

# ---- 1) pkg-config ----
PETSC_CFLAGS = $$system(PKG_CONFIG_PATH=$$PETSC_PKGCFG pkg-config --cflags petsc 2>/dev/null)
PETSC_LIBS   = $$system(PKG_CONFIG_PATH=$$PETSC_PKGCFG pkg-config --libs   petsc 2>/dev/null)

!isEmpty(PETSC_CFLAGS) {
    message("Using PETSc via pkg-config from: $$PETSC_PKGCFG")
    QMAKE_CXXFLAGS += $$PETSC_CFLAGS
    LIBS          += $$PETSC_LIBS
} else:exists($$PETSC_VARIABLES) {
    # ---- 2) petscvariables ----
    message("Using PETSc via petscvariables: $$PETSC_VARIABLES")
    include($$PETSC_VARIABLES)

    # Compile includes from PETSc (fixes petscsys.h)
    QMAKE_CXXFLAGS += $$PETSC_CC_INCLUDES

    # Important: some PETSc installs don't populate PETSC_LIB well for qmake.
    # So still force -lpetsc explicitly + PETSc externals if provided.
    LIBS += $$PETSC_EXTERNAL_LIB_BASIC $$PETSC_EXTERNAL_LIB
    LIBS += -L$$PETSC_LIBDIR -lpetsc
} else {
    # ---- 3) manual fallback ----
    message("NOTE: PETSc pkg-config and petscvariables not found; using manual include/lib flags")
    INCLUDEPATH += $$PETSC_DIR/include
    INCLUDEPATH += $$PETSC_DIR/$$PETSC_ARCH/include
    LIBS        += -L$$PETSC_LIBDIR -lpetsc
}

# runtime path
QMAKE_LFLAGS += -Wl,-rpath,$$PETSC_LIBDIR


# ============================================================
# MPI wrapper for compile & link (universal, robust)
# ============================================================
MPI_CXX = $$(MPICXX)

# Optional per-host hard overrides (only if you want):
# contains(DEFINES, PowerEdge) { isEmpty(MPI_CXX) { MPI_CXX = /mnt/3rd900/Projects/petsc-install/bin/mpicxx } }
# contains(DEFINES, Jason)     { isEmpty(MPI_CXX) { MPI_CXX = /usr/bin/mpicxx } }

# Auto-detect from PATH
isEmpty(MPI_CXX) {
    MPI_CXX = $$system(which mpicxx 2>/dev/null)
}

isEmpty(MPI_CXX) {
    error("mpicxx not found. Install mpich/openmpi or set MPICXX=/full/path/to/mpicxx")
}

!exists($$MPI_CXX) {
    error("MPI compiler not found at: $$MPI_CXX (check MPICXX or PATH)")
}

message("Using MPI_CXX = $$MPI_CXX")

QMAKE_CXX        = $$MPI_CXX
QMAKE_LINK       = $$MPI_CXX
QMAKE_LINK_SHLIB = $$MPI_CXX


# ============================================================
# Project includes
# ============================================================
INCLUDEPATH += Utilities/


# ============================================================
# Warnings / optimizations
# ============================================================
QMAKE_CXXFLAGS += -Wall -Wextra -Wpedantic -O3


# ============================================================
# Sources / Headers
# ============================================================
SOURCES += \
    Particle.cpp \
    Pathway.cpp \
    PathwaySet.cpp \
    Utilities/QuickSort.cpp \
    Utilities/Utilities.cpp \
    main.cpp \
    Utilities/Matrix.cpp \
    Utilities/Matrix_arma.cpp \
    Utilities/Vector.cpp \
    Utilities/Vector_arma.cpp \
    grid.cpp \
    petscmatrix.cpp \
    petscsolver.cpp \
    petscvector.cpp \
    plotter.cpp \
    sim_helpers.cpp \
    sim_runner.cpp

HEADERS += \
    Particle.h \
    Pathway.h \
    PathwaySet.h \
    Utilities/Matrix.h \
    Utilities/Matrix_arma.h \
    Utilities/QuickSort.h \
    Utilities/TimeSeries.h \
    Utilities/TimeSeries.hpp \
    Utilities/TimeSeriesSet.h \
    Utilities/TimeSeriesSet.hpp \
    Utilities/Utilities.h \
    Utilities/Vector.h \
    Utilities/Vector_arma.h \
    grid.h \
    petsc_init.h \
    petscmatrix.h \
    petscsolver.h \
    petscvector.h \
    plotter.h \
    sim_helpers.h \
    sim_runner.h


# ============================================================
# Armadillo / GSL / FFTW
# ============================================================
DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
LIBS    += -larmadillo -lgsl -lgslcblas -lfftw3


# ============================================================
# VTK (PowerEdge build-tree OR Jason install-prefix)
# ============================================================
GRID_USE_VTK {
    # ---- Lib path + rpath ----
    contains(DEFINES, Jason)|contains(DEFINES, WSL) {
        LIBS += -L$$VTKLIBPATH
        QMAKE_LFLAGS += -Wl,-rpath,$$VTKLIBPATH
    } else {
        LIBS += -L$$VTKBUILDPATH/lib
        QMAKE_LFLAGS += -Wl,-rpath,$$VTKBUILDPATH/lib
    }

    # ---- VTK libs (unchanged; depends on VTK_V) ----
    LIBS += -lvtkChartsCore$$VTK_V \
            -lvtkCommonColor$$VTK_V \
            -lvtkCommonComputationalGeometry$$VTK_V \
            -lvtkCommonCore$$VTK_V \
            -lvtkCommonDataModel$$VTK_V \
            -lvtkCommonExecutionModel$$VTK_V \
            -lvtkCommonMath$$VTK_V \
            -lvtkCommonMisc$$VTK_V \
            -lvtkCommonSystem$$VTK_V \
            -lvtkCommonTransforms$$VTK_V \
            -lvtkDICOMParser$$VTK_V \
            -lvtkexpat$$VTK_V \
            -lvtkFiltersAMR$$VTK_V \
            -lvtkFiltersCore$$VTK_V \
            -lvtkFiltersExtraction$$VTK_V \
            -lvtkFiltersFlowPaths$$VTK_V \
            -lvtkFiltersGeneral$$VTK_V \
            -lvtkFiltersGeneric$$VTK_V \
            -lvtkFiltersGeometry$$VTK_V \
            -lvtkFiltersHybrid$$VTK_V \
            -lvtkFiltersHyperTree$$VTK_V \
            -lvtkFiltersImaging$$VTK_V \
            -lvtkFiltersModeling$$VTK_V \
            -lvtkFiltersParallel$$VTK_V \
            -lvtkFiltersParallelImaging$$VTK_V \
            -lvtkFiltersPoints$$VTK_V \
            -lvtkFiltersProgrammable$$VTK_V \
            -lvtkFiltersSelection$$VTK_V \
            -lvtkFiltersSMP$$VTK_V \
            -lvtkFiltersSources$$VTK_V \
            -lvtkFiltersStatistics$$VTK_V \
            -lvtkFiltersTexture$$VTK_V \
            -lvtkFiltersTopology$$VTK_V \
            -lvtkFiltersVerdict$$VTK_V \
            -lvtkfreetype$$VTK_V \
            -lvtkGeovisCore$$VTK_V \
            -lvtkgl2ps$$VTK_V \
            -lvtkhdf5$$VTK_V \
            -lvtkhdf5_hl$$VTK_V \
            -lvtkImagingColor$$VTK_V \
            -lvtkImagingCore$$VTK_V \
            -lvtkImagingFourier$$VTK_V \
            -lvtkImagingGeneral$$VTK_V \
            -lvtkImagingHybrid$$VTK_V \
            -lvtkImagingMath$$VTK_V \
            -lvtkImagingMorphological$$VTK_V \
            -lvtkImagingSources$$VTK_V \
            -lvtkImagingStatistics$$VTK_V \
            -lvtkImagingStencil$$VTK_V \
            -lvtkInfovisCore$$VTK_V \
            -lvtkInfovisLayout$$VTK_V \
            -lvtkInteractionImage$$VTK_V \
            -lvtkInteractionStyle$$VTK_V \
            -lvtkInteractionWidgets$$VTK_V \
            -lvtkIOAMR$$VTK_V \
            -lvtkIOCore$$VTK_V \
            -lvtkIOEnSight$$VTK_V \
            -lvtkIOExodus$$VTK_V \
            -lvtkIOGeometry$$VTK_V \
            -lvtkIOImage$$VTK_V \
            -lvtkIOImport$$VTK_V \
            -lvtkIOInfovis$$VTK_V \
            -lvtkIOLegacy$$VTK_V \
            -lvtkIOLSDyna$$VTK_V \
            -lvtkIOMINC$$VTK_V \
            -lvtkIOMovie$$VTK_V \
            -lvtkIONetCDF$$VTK_V \
            -lvtkIOParallel$$VTK_V \
            -lvtkIOParallelXML$$VTK_V \
            -lvtkIOPLY$$VTK_V \
            -lvtkIOSQL$$VTK_V \
            -lvtkIOTecplotTable$$VTK_V \
            -lvtkIOVideo$$VTK_V \
            -lvtkIOXML$$VTK_V \
            -lvtkIOXMLParser$$VTK_V \
            -lvtkjpeg$$VTK_V \
            -lvtkjsoncpp$$VTK_V \
            -lvtklibharu$$VTK_V \
            -lvtklibxml2$$VTK_V \
            -lvtklz4$$VTK_V \
            -lvtkmetaio$$VTK_V \
            -lvtkParallelCore$$VTK_V \
            -lvtkpng$$VTK_V \
            -lvtkRenderingAnnotation$$VTK_V \
            -lvtkRenderingContext2D$$VTK_V \
            -lvtkRenderingCore$$VTK_V \
            -lvtkRenderingFreeType$$VTK_V \
            -lvtkRenderingGL2PSOpenGL2$$VTK_V \
            -lvtkRenderingImage$$VTK_V \
            -lvtkRenderingLabel$$VTK_V \
            -lvtkRenderingLOD$$VTK_V \
            -lvtkRenderingOpenGL2$$VTK_V \
            -lvtkRenderingVolume$$VTK_V \
            -lvtkRenderingVolumeOpenGL2$$VTK_V \
            -lvtksqlite$$VTK_V \
            -lvtksys$$VTK_V \
            -lvtktiff$$VTK_V \
            -lvtkverdict$$VTK_V \
            -lvtkViewsContext2D$$VTK_V \
            -lvtkViewsCore$$VTK_V \
            -lvtkViewsInfovis$$VTK_V \
            -lvtkzlib$$VTK_V

    # Base include is enough for Jason install-prefix
    INCLUDEPATH += $$VTKHEADERPATH

    # PowerEdge build-tree needs extra include dirs; keep only for PowerEdge
    contains(DEFINES, PowerEdge) {
        INCLUDEPATH += $$VTKHEADERPATH \
                       $$VTKHEADERPATH/Common/Core \
                       $$VTKBUILDPATH/Common/Core \
                       $$VTKHEADERPATH/Common/Color \
                       $$VTKHEADERPATH/Common/Transforms \
                       $$VTKBUILDPATH/Common/Transforms \
                       $$VTKBUILDPATH/Common/Color \
                       $$VTKBUILDPATH/Common/DataModel \
                       $$VTKBUILDPATH/Utilities/KWIML \
                       $$VTKHEADERPATH/Utilities/KWIML \
                       $$VTKHEADERPATH/Rendering/Core \
                       $$VTKBUILDPATH/Rendering/Core \
                       $$VTKBUILDPATH/Filters/Core \
                       $$VTKHEADERPATH/Charts/Core \
                       $$VTKBUILDPATH/Charts/Core \
                       $$VTKBUILDPATH/Filters/General \
                       $$VTKBUILDPATH/Rendering/Context2D \
                       $$VTKHEADERPATH/Rendering/Context2D \
                       $$VTKHEADERPATH/Common/DataModel \
                       $$VTKHEADERPATH/Common/Math \
                       $$VTKHEADERPATH/Views/Context2D \
                       $$VTKBUILDPATH/Views/Context2D \
                       $$VTKBUILDPATH/Views/Core \
                       $$VTKBUILDPATH/Interaction/Widgets \
                       $$VTKHEADERPATH/Views/Core \
                       $$VTKHEADERPATH/Interaction/Style \
                       $$VTKBUILDPATH/Interaction/Style \
                       $$VTKHEADERPATH/Filters/Modeling \
                       $$VTKBUILDPATH/Filters/Modeling \
                       $$VTKHEADERPATH/Common/ExecutionModel \
                       $$VTKBUILDPATH/Common/ExecutionModel \
                       $$VTKHEADERPATH/Interaction/Widgets \
                       $$VTKHEADERPATH/Filters/Core \
                       $$VTKHEADERPATH/Common/Misc \
                       $$VTKBUILDPATH/Common/Misc \
                       $$VTKHEADERPATH/IO/XML \
                       $$VTKBUILDPATH/IO/XML \
                       $$VTKHEADERPATH/Filters/Sources \
                       $$VTKBUILDPATH/Filters/Sources \
                       $$VTKHEADERPATH/Filters/General \
                       $$VTKHEADERPATH/IO/Image \
                       $$VTKBUILDPATH/IO/Image \
                       $$VTKHEADERPATH/Imaging/Core \
                       $$VTKBUILDPATH/Imaging/Core \
                       $$VTKBUILDPATH/Utilities/KWSys \
                       $$VTKBUILDPATH/ThirdParty/nlohmannjson \
                       $$VTKHEADERPATH/ThirdParty/nlohmannjson \
                       $$VTKBUILDPATH/Common/Math
    }
}

message("Final INCLUDEPATH = $$INCLUDEPATH")
message("Final LIBS        = $$LIBS")
message("Final QMAKE_LFLAGS= $$QMAKE_LFLAGS")
