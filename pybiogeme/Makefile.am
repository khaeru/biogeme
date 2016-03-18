ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils  -I$(top_srcdir)/libraries/randomNumbers -I$(top_srcdir)/libraries/cfsqp -I$(top_srcdir)/libraries/donlp2 -I$(top_srcdir)/libraries/solvopt -I$(top_srcdir)/libraries/linearAlgebra -I$(top_srcdir)/libraries/trustRegion -I$(top_srcdir)/libraries/pythonbiogeme $(PYTHON_INCLUDES)

pybiogeme_CXXFLAGS = $(PYTHON_CXXFLAGS)
pybiogeme_LDADD = $(shell pwd)/../libraries/libpythonbiogeme.la
pybiogeme_LDFLAGS = $(PYTHON_LDFLAGS)
if APPLE
pybiogeme_LDFLAGS += -Wl,-stack_size -Wl,1000000
endif
LIBS = $(PYTHON_LIBS)

pybiogeme_static_LDFLAGS = -static $(PYTHON_LDFLAGS)
if APPLE
pybiogeme_static_LDFLAGS += -Wl,-stack_size -Wl,1000000
endif
pybiogeme_static_SOURCES = $(pybiogeme_SOURCES)
pybiogeme_static_LDADD = $(shell pwd)/../libraries/pythonbiogeme/libpythonbiogeme.la   $(shell pwd)/../libraries/gaussHermite/libgaussHermite.la  $(shell pwd)/../libraries/randomNumbers/librandomNumberspython.la $(shell pwd)/../libraries/cfsqp/libcfsqppython.la $(shell pwd)/../libraries/solvopt/libsolvoptpython.la   $(shell pwd)/../libraries/trustRegion/libtrustRegionpython.la  $(shell pwd)/../libraries/linearAlgebra/liblinearAlgebrapython.la  $(shell pwd)/../libraries/utils/libutilspython.la $(shell pwd)/../libraries/donlp2/libdonlp2python.la

#$(shell pwd)/../libraries/donlp2/libdonlp2python.la 

if BUILD_PYTHON
bin_PROGRAMS = pybiogeme
if BUILD_STATIC
bin_PROGRAMS += pybiogeme_static
endif
if INSTALL_VERSIONED
install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	for p in $(bin_PROGRAMS); do \
	mv -f $${p}$(EXEEXT) $${p}$(INSTALL_VERSION_SUFFIX)$(EXEEXT) && \
	$(LN_S) $${p}$(INSTALL_VERSION_SUFFIX)$(EXEEXT) $${p}$(EXEEXT); \
	done
endif
endif

pybiogeme_SOURCES = \
patDisplay.cc \
pybiogeme.cc

