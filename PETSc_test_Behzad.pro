TEMPLATE = app
CONFIG  += c++17 console
QT      -= core gui widgets

DEFINES += GRID_USE_VTK
CONFIG  += GRID_USE_VTK
DEFINES += GSL

# ==== VTK paths ====
VTKBUILDPATH  = /home/arash/Projects/VTK-build
VTKHEADERPATH = /home/arash/Projects/VTK
VTK_V         = -9.0

# ==== PETSc (build-tree layout) ====
PETSC_DIR  = $$(PETSC_DIR)
PETSC_ARCH = $$(PETSC_ARCH)

message("PETSC_DIR from env = $$PETSC_DIR")
message("PETSC_ARCH from env = $$PETSC_ARCH")

isEmpty(PETSC_DIR) {
    PETSC_DIR = /home/arash/petsc
}
isEmpty(PETSC_ARCH) {
    PETSC_ARCH = arch-linux-c-opt
}

INCLUDEPATH += $$PETSC_DIR/include
INCLUDEPATH += $$PETSC_DIR/$$PETSC_ARCH/include

LIBS += -L$$PETSC_DIR/$$PETSC_ARCH/lib -lpetsc
QMAKE_LFLAGS += -Wl,-rpath,$$PETSC_DIR/$$PETSC_ARCH/lib


# ==== MPI wrapper for compile & link ====
# ==== MPI wrapper (use PETSc's MPICH) ====
QMAKE_CXX        = /home/arash/petsc-install/bin/mpicxx
QMAKE_LINK       = /home/arash/petsc-install/bin/mpicxx
QMAKE_LINK_SHLIB = /home/arash/petsc-install/bin/mpicxx
# do NOT add -lmpi manually — mpicxx already adds MPI libs

# ==== Warnings / optimizations ====
QMAKE_CXXFLAGS += -Wall -Wextra -Wpedantic -O3

# ==== Sources / Headers ====
SOURCES += \
    main.cpp \
    Matrix.cpp \
    Matrix_arma.cpp \
    Vector.cpp \
    Vector_arma.cpp \
    grid.cpp \
    petscmatrix.cpp \
    petscsolver.cpp \
    petscvector.cpp

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

# ==== Armadillo / GSL ====
DEFINES += ARMA_USE_LAPACK ARMA_USE_BLAS
LIBS    += -larmadillo -lgsl -lgslcblas

# ==== VTK ====
GRID_USE_VTK {
    LIBS += -L$$VTKBUILDPATH/lib
    QMAKE_LFLAGS += -Wl,-rpath,$$VTKBUILDPATH/lib

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

            INCLUDEPATH +=$${VTKHEADERPATH}
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Common/Core
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/Color
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/Transforms
            INCLUDEPATH +=$${VTKBUILDPATH}/Common/Transforms
            INCLUDEPATH +=$${VTKBUILDPATH}/Common/Color
            INCLUDEPATH +=$${VTKBUILDPATH}/Common/DataModel
            INCLUDEPATH +=$${VTKBUILDPATH}/Utilities/KWIML
            INCLUDEPATH +=$${VTKHEADERPATH}/Utilities/KWIML
            INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Core
            INCLUDEPATH +=$${VTKHEADERPATH}/Charts/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Charts/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Filters/General
            INCLUDEPATH +=$${VTKBUILDPATH}/Rendering/Context2D
            INCLUDEPATH +=$${VTKHEADERPATH}/Rendering/Context2D
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/DataModel
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/Math
            INCLUDEPATH +=$${VTKHEADERPATH}/Views/Context2D
            INCLUDEPATH +=$${VTKBUILDPATH}/Views/Context2D
            INCLUDEPATH +=$${VTKBUILDPATH}/Views/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Widgets
            INCLUDEPATH +=$${VTKHEADERPATH}/Views/Core
            INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Style
            INCLUDEPATH +=$${VTKBUILDPATH}/Interaction/Style
            INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Modeling
            INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Modeling
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/ExecutionModel
            INCLUDEPATH +=$${VTKBUILDPATH}/Common/ExecutionModel
            INCLUDEPATH +=$${VTKHEADERPATH}/Interaction/Widgets/
            INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Core/
            INCLUDEPATH +=$${VTKHEADERPATH}/Common/Misc/
            INCLUDEPATH +=$${VTKBUILDPATH}/Common/Misc
            INCLUDEPATH +=$${VTKHEADERPATH}/IO/XML/
            INCLUDEPATH +=$${VTKBUILDPATH}/IO/XML
            INCLUDEPATH +=$${VTKHEADERPATH}/Filters/Sources
            INCLUDEPATH +=$${VTKBUILDPATH}/Filters/Sources
            INCLUDEPATH +=$${VTKHEADERPATH}/Filters/General
            INCLUDEPATH +=$${VTKHEADERPATH}/IO/Image
            INCLUDEPATH +=$${VTKBUILDPATH}/IO/Image
            INCLUDEPATH +=$${VTKHEADERPATH}/Imaging/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Imaging/Core
            INCLUDEPATH +=$${VTKBUILDPATH}/Utilities/KWSys
            INCLUDEPATH += $${VTKBUILDPATH}/ThirdParty/nlohmannjson
            INCLUDEPATH += $${VTKHEADERPATH}/ThirdParty/nlohmannjson
            INCLUDEPATH += $${VTKBUILDPATH}/Common/Math
}

message("Final INCLUDEPATH = $$INCLUDEPATH")
message("Final LIBS        = $$LIBS")
message("Final QMAKE_LFLAGS= $$QMAKE_LFLAGS")

