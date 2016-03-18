import sys
import sysconfig
option=sys.argv[1].lstrip('-')
if option=='includes':
 print('-I' + sysconfig.get_path('include') + ' -I' + sysconfig.get_path('platinclude'))
if option=='cflags':
 flags=sysconfig.get_config_var('CFLAGS').split()
 for unflag in ['-Wno-unused-result','-Wstrict-prototypes']:
  if flags.count(unflag) > 0:
   flags.remove(unflag)
 print(' '.join(flags))
if option=='libs':
 print('-lpython' + sysconfig.get_config_var('VERSION') + sys.abiflags + ' ' + sysconfig.get_config_var('LIBS') + ' ' + sysconfig.get_config_var('SYSLIBS'))
if option=='ldflags':
 flags=sysconfig.get_config_var('LDFLAGS')
 if not sysconfig.get_config_var('Py_ENABLE_SHARED'):
  flags+=' -L'+sysconfig.get_config_var('LIBPL')
 flags+=' ' + sysconfig.get_config_var('LINKFORSHARED')
 print(flags)
if option=='headerdir':
 print(sysconfig.get_path('include'))
if option=='version':
 print(sysconfig.get_config_var('VERSION'))
