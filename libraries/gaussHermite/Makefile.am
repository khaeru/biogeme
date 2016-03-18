ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils 

noinst_LTLIBRARIES =

if BUILD_PYTHON
noinst_LTLIBRARIES += libgaussHermite.la
endif

if BUILD_STATIC
libgaussHermite_la_LDFLAGS = -static
endif

libgaussHermite_la_CXXFLAGS = $(PYTHON_CXXFLAGS)

libgaussHermite_la_SOURCES = \
bioGaussHermite.cc \
bioGhFunction.cc \
bioGaussHermite.h \
bioGhFunction.h
