ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/linearAlgebra  -I$(top_srcdir)/libraries/trustRegion 

noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += libsolvopt.la
endif

if BUILD_GUI
noinst_LIBRARIES = libsolvoptfltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += libsolvoptpython.la
endif

if BUILD_STATIC
libsolvopt_la_LDFLAGS = -static
libsolvoptpython_la_LDFLAGS = -static
endif

libsolvoptfltk_a_SOURCES = $(libsolvopt_la_SOURCES)
libsolvoptfltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

libsolvoptpython_la_SOURCES = $(libsolvopt_la_SOURCES)
libsolvoptpython_la_CXXFLAGS = $(PYTHON_CXXFLAGS)

libsolvopt_la_SOURCES = \
myapprgrdn.cc \
mysolvopt.cc \
patSolvOpt.cc \
patSolvOpt.h \
solvopt.h \
solvoptParameters.h 
