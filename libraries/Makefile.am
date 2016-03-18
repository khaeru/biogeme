ACLOCAL_AMFLAGS = -I m4
pkglibrdir = ${pkglibdir}$(INSTALL_VERSION_SUFFIX)

pkglibr_LTLIBRARIES = 

if BUILD_BISONBIOGEME
pkglibr_LTLIBRARIES += libbiogeme.la
endif

if BUILD_PYTHON
pkglibr_LTLIBRARIES += libpythonbiogeme.la
endif

# add here the libs that are needed for everything that is compiled natively
libbiogeme_la_LIBADD = \
           $(shell pwd)/bisonbiogeme/libbisonbiogeme.la \
           $(shell pwd)/randomNumbers/librandomNumbers.la \
           $(shell pwd)/cfsqp/libcfsqp.la \
           $(shell pwd)/donlp2/libdonlp2.la \
           $(shell pwd)/solvopt/libsolvopt.la \
           $(shell pwd)/trustRegion/libtrustRegion.la \
           $(shell pwd)/linearAlgebra/liblinearAlgebra.la \
           $(shell pwd)/parameters/libparameters.la \
           $(shell pwd)/utils/libutils.la

libbiogeme_la_SOURCES =

libbiogeme_la_LDFLAGS = $(INSTALL_RELEASE)

libbiogeme_la_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/randomNumbers -I$(top_srcdir)/libraries/parameters -I$(top_srcdir)/libraries/cfsqp -I$(top_srcdir)/libraries/donlp2 -I$(top_srcdir)/libraries/solvopt -I$(top_srcdir)/libraries/linearAlgebra -I$(top_srcdir)/libraries/trustRegion -I$(top_srcdir)/libraries/bisonbiogeme 

libpythonbiogeme_la_LIBADD = \
           $(shell pwd)/pythonbiogeme/libpythonbiogeme.la \
           $(shell pwd)/randomNumbers/librandomNumberspython.la \
           $(shell pwd)/cfsqp/libcfsqppython.la \
           $(shell pwd)/donlp2/libdonlp2python.la \
           $(shell pwd)/solvopt/libsolvoptpython.la \
           $(shell pwd)/trustRegion/libtrustRegionpython.la \
           $(shell pwd)/gaussHermite/libgaussHermite.la \
           $(shell pwd)/linearAlgebra/liblinearAlgebrapython.la \
           $(shell pwd)/utils/libutilspython.la

libpythonbiogeme_la_SOURCES =

libpythonbiogeme_la_LDFLAGS = $(INSTALL_RELEASE) $(PYTHON_LDFLAGS)

libpythonbiogeme_la_CPPFLAGS = -I$(top_srcdir)/libraries/utils -I$(top_srcdir)/libraries/randomNumbers -I$(top_srcdir)/libraries/cfsqp -I$(top_srcdir)/libraries/donlp2 -I$(top_srcdir)/libraries/solvopt -I$(top_srcdir)/libraries/gaussHermite -I$(top_srcdir)/libraries/linearAlgebra -I$(top_srcdir)/libraries/trustRegion -I$(top_srcdir)/libraries/pythonbiogeme

