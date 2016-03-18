ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/linearAlgebra

noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += librandomNumbers.la
endif

if BUILD_GUI
noinst_LIBRARIES = librandomNumbersfltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += librandomNumberspython.la
endif

if BUILD_STATIC
librandomNumbers_la_LDFLAGS = -static
librandomNumberspython_la_LDFLAGS = -static
endif

librandomNumbersfltk_a_SOURCES = $(librandomNumbers_la_SOURCES)
librandomNumbersfltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

librandomNumberspython_la_SOURCES = $(librandomNumbers_la_SOURCES)
librandomNumberspython_la_CXXFLAGS = $(PYTHON_CXXFLAGS)

librandomNumbers_la_SOURCES =  \
patRandomNumberGenerator.cc \
patBinomial.cc \
patCenteredUniform.cc \
patDiscreteDistribution.cc \
patDrawsFromFile.cc \
patGenerateNormalDraws.cc \
patHalton.cc \
patHessTrain.cc \
patNormal.cc \
patNormalWichura.cc \
patRandomInteger.cc \
patUnixUniform.cc \
patUniform.cc \
patZeroOneUniform.cc \
patMultivariateNormal.cc \
patBinomial.h \
patCenteredUniform.h \
patDiscreteDistribution.h \
patDrawsFromFile.h \
patGenerateNormalDraws.h \
patHalton.h \
patHessTrain.h \
patNormal.h \
patNormalWichura.h \
patRandomInteger.h \
patRandomNumberGenerator.h \
patUniform.h \
patUnixUniform.h \
patZeroOneUniform.h \
patMultivariateNormal.h
