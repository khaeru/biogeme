ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/linearAlgebra


noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += libtrustRegion.la
endif

if BUILD_GUI
noinst_LIBRARIES = libtrustRegionfltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += libtrustRegionpython.la
endif

if BUILD_STATIC
libtrustRegion_la_LDFLAGS = -static
libtrustRegionpython_la_LDFLAGS = -static
endif

libtrustRegionfltk_a_SOURCES = $(libtrustRegion_la_SOURCES)
libtrustRegionfltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

libtrustRegionpython_la_SOURCES = $(libtrustRegion_la_SOURCES)
libtrustRegionpython_la_CXXFLAGS = $(PYTHON_CXXFLAGS)


libtrustRegion_la_SOURCES =  \
patNonLinearProblem.cc \
patSimBasedOptimizationProblem.cc \
trBFGS.cc \
trBounds.cc \
trBoxGCP.cc \
trFunction.cc \
trHessian.cc \
trHybridMatrix.cc \
trLineSearchAlgo.cc \
trNonLinearAlgo.cc \
trPrecond.cc \
trPrecondCG.cc \
trQuadratic.cc \
trSR1.cc \
trSchnabelEskow.cc \
trSecantUpdate.cc \
trSimBasedSimpleBoundsAlgo.cc \
trSimpleBoundsAlgo.cc \
trTheoryGCP.cc \
trTointSteihaug.cc \
trParameters.cc \
patNonLinearProblem.h \
patSimBasedOptimizationProblem.h \
trBFGS.h \
trBounds.h \
trBoxGCP.h \
trFunction.h \
trHessian.h \
trHybridMatrix.h \
trLineSearchAlgo.h \
trMatrixVector.h \
trParameters.h \
trPrecond.h \
trPrecondCG.h \
trQuadratic.h \
trSR1.h \
trSchnabelEskow.h \
trSecantUpdate.h \
trSimBasedSimpleBoundsAlgo.h \
trSimpleBoundsAlgo.h \
trTheoryGCP.h \
trTointSteihaug.h \
trVector.h 
