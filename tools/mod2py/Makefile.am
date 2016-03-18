ACLOCAL_AMFLAGS = -I m4
pkgdatdir = ${pkgdatadir}$(INSTALL_VERSION_SUFFIX)
bin_SCRIPTS = mod2py

if INSTALL_VERSIONED
install-exec-hook:
	cd $(DESTDIR)$(bindir) && \
	for p in $(bin_SCRIPTS); do \
	mv -f $${p}$(EXEEXT) $${p}$(INSTALL_VERSION_SUFFIX)$(EXEEXT) && \
	$(LN_S) $${p}$(INSTALL_VERSION_SUFFIX)$(EXEEXT) $${p}$(EXEEXT); \
	done
endif

CLEANFILES = $(bin_SCRIPTS)


