ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils


noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += liblinearAlgebra.la
endif

if BUILD_GUI
noinst_LIBRARIES = liblinearAlgebrafltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += liblinearAlgebrapython.la
endif

if BUILD_STATIC
liblinearAlgebra_la_LDFLAGS = -static
liblinearAlgebrapython_la_LDFLAGS = -static
endif

liblinearAlgebrafltk_a_SOURCES = $(liblinearAlgebra_la_SOURCES)
liblinearAlgebrafltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

liblinearAlgebrapython_la_SOURCES = $(liblinearAlgebra_la_SOURCES)
liblinearAlgebrapython_la_CXXFLAGS = $(PYTHON_CXXFLAGS)

liblinearAlgebra_la_SOURCES =  \
patEigenVectors.cc \
patGeneralizedInverseIteration.cc \
patHouseholder.cc \
patHybridMatrix.cc \
patInverse.cc \
patInverseIteration.cc \
patLu.cc \
patMatrix.cc \
patMyMatrix.cc \
patPythag.cc \
patQR.cc \
patSvd.cc \
patKalman.cc \
patEigenVectors.h \
patGeneralizedInverseIteration.h \
patHouseholder.h \
patHybridMatrix.h \
patInverse.h \
patInverseIteration.h \
patLu.h \
patMatrix.h \
patMyMatrix.h \
patNrSgn.h \
patPythag.h \
patQR.h \
patSvd.h \
patKalman.h 