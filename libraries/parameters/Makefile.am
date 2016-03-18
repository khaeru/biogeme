ACLOCAL_AMFLAGS = -I m4
AM_YFLAGS = -d
AM_LFLAGS = -+ -olex.yy.c

AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils 

if BUILD_BISONBIOGEME
noinst_LTLIBRARIES = libparameters.la
endif

libparameters_la_SOURCES = patParserParam.yy patScannerParam.ll patParameters.cc  patParameters.h

if BUILD_GUI 
noinst_LIBRARIES = libparametersfltk.a
endif

if BUILD_STATIC
libparameters_la_LDFLAGS = -static
endif

libparametersfltk_a_SOURCES = $(libparameters_la_SOURCES)
libparametersfltk_a_CXXFLAGS = $(FLTK_CXXFLAGS)


