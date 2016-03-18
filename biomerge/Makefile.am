ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils  -I$(top_srcdir)/libraries/randomNumbers -I$(top_srcdir)/libraries/parameters -I$(top_srcdir)/libraries/cfsqp -I$(top_srcdir)/libraries/donlp2 -I$(top_srcdir)/libraries/solvopt -I$(top_srcdir)/libraries/linearAlgebra -I$(top_srcdir)/libraries/trustRegion -I$(top_srcdir)/libraries/bisonbiogeme

biomerge_LDADD =  $(shell pwd)/../libraries/libbiogeme.la

biomerge_static_LDFLAGS = -static  -static-libstdc++ -static-libgcc
biomerge_static_SOURCES = $(biomerge_SOURCES)


biomerge_static_LDADD =  $(shell pwd)/../libraries/bisonbiogeme/libbisonbiogeme.la $(shell pwd)/../libraries/randomNumbers/librandomNumbers.la $(shell pwd)/../libraries/cfsqp/libcfsqp.la $(shell pwd)/../libraries/donlp2/libdonlp2.la $(shell pwd)/../libraries/solvopt/libsolvopt.la $(shell pwd)/../libraries/trustRegion/libtrustRegion.la  $(shell pwd)/../libraries/linearAlgebra/liblinearAlgebra.la $(shell pwd)/../libraries/parameters/libparameters.la $(shell pwd)/../libraries/utils/libutils.la 


if BUILD_BISONBIOGEME
bin_PROGRAMS = biomerge 
if BUILD_STATIC
bin_PROGRAMS += biomerge_static
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

biomerge_SOURCES = biomerge.cc patDisplay.cc

