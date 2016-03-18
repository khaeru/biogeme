ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -I$(top_srcdir)/libraries/utils  -I$(top_srcdir)/libraries/randomNumbers -I$(top_srcdir)/libraries/parameters -I$(top_srcdir)/libraries/cfsqp -I$(top_srcdir)/libraries/donlp2 -I$(top_srcdir)/libraries/solvopt -I$(top_srcdir)/libraries/linearAlgebra -I$(top_srcdir)/libraries/trustRegion -I$(top_srcdir)/libraries/bisonbiogeme

guibiogeme_CXXFLAGS = $(FLTK_CXXFLAGS)

guibiogeme_LDADD = $(FLTK_LDFLAGS) $(shell pwd)/../libraries/libbiogeme.la

guibiogeme_static_CXXFLAGS = $(guibiogeme_CXXFLAGS)
guibiogeme_static_LDFLAGS = -static -static-libgcc -static-libstdc++
guibiogeme_static_SOURCES = $(guibiogeme_SOURCES)
guibiogeme_static_LDADD = $(FLTK_LDSTATICFLAGS) $(shell pwd)/../libraries/bisonbiogeme/libbisonbiogemefltk.a  $(shell pwd)/../libraries/randomNumbers/librandomNumbersfltk.a $(shell pwd)/../libraries/cfsqp/libcfsqpfltk.a $(shell pwd)/../libraries/donlp2/libdonlp2fltk.a $(shell pwd)/../libraries/solvopt/libsolvoptfltk.a $(shell pwd)/../libraries/trustRegion/libtrustRegionfltk.a $(shell pwd)/../libraries/linearAlgebra/liblinearAlgebrafltk.a   $(shell pwd)/../libraries/parameters/libparametersfltk.a  $(shell pwd)/../libraries/utils/libutilsfltk.a

if BUILD_GUI
if BUILD_BISONBIOGEME
bin_PROGRAMS = guibiogeme
if BUILD_STATIC
bin_PROGRAMS += guibiogeme_static
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
endif
guibiogeme_SOURCES = \
patWinBiogeme.cc \
patDisplay.cc \
winbiogeme.cc \
patWinBiogeme.h
