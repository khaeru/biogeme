ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/linearAlgebra  -I$(top_srcdir)/libraries/trustRegion 

noinst_LTLIBRARIES =

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES += libdonlp2.la
endif

if BUILD_GUI
noinst_LIBRARIES = libdonlp2fltk.a
endif

if BUILD_PYTHON
noinst_LTLIBRARIES += libdonlp2python.la
endif

if BUILD_STATIC
libdonlp2_la_LDFLAGS = -static
libdonlp2python_la_LDFLAGS = -static
endif

libdonlp2fltk_a_SOURCES = $(libdonlp2_la_SOURCES)
libdonlp2fltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)

libdonlp2python_la_SOURCES = $(libdonlp2_la_SOURCES)
libdonlp2python_la_CXXFLAGS = $(PYTHON_CXXFLAGS)

libdonlp2_la_SOURCES = \
mydonlp2.cc \
myuser_eval.cc \
patDonlp2.cc \
donlp2Parameters.h \
o8comm.h \
o8cons.h \
o8fint.h \
o8fuco.h \
o8gene.h \
o8para.h \
patDonlp2.h

