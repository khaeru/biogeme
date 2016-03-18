ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/linearAlgebra  -I$(top_srcdir)/libraries/trustRegion 

noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += libcfsqp.la
endif

if BUILD_GUI
noinst_LIBRARIES = libcfsqpfltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += libcfsqppython.la
endif

if BUILD_STATIC
libcfsqp_la_LDFLAGS = -static
libcfsqppython_la_LDFLAGS = -static
endif

libcfsqpfltk_a_SOURCES = $(libcfsqp_la_SOURCES)
libcfsqpfltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

libcfsqppython_la_SOURCES = $(libcfsqp_la_SOURCES)
libcfsqppython_la_CXXFLAGS = $(PYTHON_CXXFLAGS)

libcfsqp_la_SOURCES = \
mycfsqp.cc \
myqld.cc \
patCfsqp.cc \
cfsqpusr.h \
patCfsqp.h 