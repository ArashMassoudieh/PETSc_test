TEMPLATE = app
CONFIG  += c++17 console
QT      -= core gui widgets

DEFINES += GRID_USE_VTK
CONFIG  += GRID_USE_VTK
DEFINES += GSL

# ============================================================
# Host-config (PowerEdge active)
# ============================================================
#CONFIG += Behzad
#DEFINES += Behzad

CONFIG += PowerEdge
DEFINES += PowerEdge

#CONFIG += Arash
#DEFINES += Arash

#CONFIG += SligoCreek
#DEFINES += SligoCreek

#CONFIG  += Jason
#DEFINES += Jason

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

# Keep Jason/WSL blocks here if you want (inactive now)
contains(DEFINES, Jason) {
    VTKBUILDPATH   = /home/arash/Projects/VTK-9.3.0/build
    VTKHEADERPATH  = /home/arash/Projects/VTK-9.3.0
    VTK_V          = -9.3
}

contains(DEFINES, WSL) {
    VTKINSTALLPATH = /home/behzad/Projects/VTK-install
    VTK_V          = -9.4
    VTK_INCNAME    = vtk-9.4

    !exists($$VTKINSTALLPATH/include/$$VTK_INCNAME/vtkSmartPointer.h) {
        error("VTK not found: $$VTKINSTALLPATH/include/$$VTK_INCNAME/vtkSmartPointer.h")
    }
}


# ============================================================
# PETSc (PowerEdge)  -- prefer petscvariables if present
# ============================================================
PETSC_DIR  = $$(PETSC_DIR)
PETSC_ARCH = $$(PETSC_ARCH)

message("PETSC_DIR from env  = $$PETSC_DIR")
message("PETSC_ARCH from env = $$PETSC_ARCH")

contains(DEFINES, PowerEdge) {
    # ### EDIT ME (PowerEdge PETSc root)
    isEmpty(PETSC_DIR)  { PETSC_DIR  = /mnt/3rd900/Projects/petsc }
    # ### EDIT ME (PowerEdge PETSc arch dir name)
    isEmpty(PETSC_ARCH) { PETSC_ARCH = arch-linux-c-opt }
} else {
    # fallback (kept for safety)
    isEmpty(PETSC_DIR)  { PETSC_DIR  = /home/arash/petsc }
    isEmpty(PETSC_ARCH) { PETSC_ARCH = arch-linux-c-opt }
}

message("PETSC_DIR final  = $$PETSC_DIR")
message("PETSC_ARCH final = $$PETSC_ARCH")

PETSC_VARIABLES = $$PETSC_DIR/$$PETSC_ARCH/lib/petsc/conf/petscvariables

exists($$PETSC_VARIABLES) {
    message("Including PETSc variables from: $$PETSC_VARIABLES")
    include($$PETSC_VARIABLES)

    # Compile flags & include paths from PETSc
    QMAKE_CXXFLAGS += $$PETSC_CC_INCLUDES

    # Link PETSc + its external deps
    LIBS += $$PETSC_LIB $$PETSC_EXTERNAL_LIB_BASIC $$PETSC_EXTERNAL_LIB
} else {
    # Fallback: manual include/lib (your current working style)
    message("NOTE: PETSc petscvariables not found; using manual include/lib flags")
    INCLUDEPATH += $$PETSC_DIR/include
    INCLUDEPATH += $$PETSC_DIR/$$PETSC_ARCH/include
    LIBS += -L$$PETSC_DIR/$$PETSC_ARCH/lib -lpetsc
}

# Ensure runtime finds PETSc shared libs
QMAKE_LFLAGS += -Wl,-rpath,$$PETSC_DIR/$$PETSC_ARCH/lib

# Project includes
INCLUDEPATH += Utilities/


# ============================================================
# MPI wrapper for compile & link (PowerEdge)
# ============================================================
# Prefer environment MPICXX if set; else host defaults.
MPI_CXX = $$(MPICXX)

contains(DEFINES, PowerEdge) {
    # ### EDIT ME if you truly want PETSc-bundled mpicxx on PowerEdge:
    # MPI_CXX = /mnt/3rd900/Projects/petsc-install/bin/mpicxx
}

isEmpty(MPI_CXX) {
    MPI_CXX = mpicxx
}

message("Using MPI_CXX = $$MPI_CXX")

QMAKE_CXX        = $$MPI_CXX
QMAKE_LINK       = $$MPI_CXX
QMAKE_LINK_SHLIB = $$MPI_CXX


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
# VTK (PowerEdge build-tree)
# ============================================================
GRID_USE_VTK {
    LIBS += -L$$VTKBUILDPATH/lib
    QMAKE_LFLAGS += -Wl,-rpath,$$VTKBUILDPATH/lib

    # VTK libs (as you had; depends on VTK_V)
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

    # Keep your include spam (since it’s proven working on build-tree VTK)
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

message("Final INCLUDEPATH = $$INCLUDEPATH")
message("Final LIBS        = $$LIBS")
message("Final QMAKE_LFLAGS= $$QMAKE_LFLAGS")
