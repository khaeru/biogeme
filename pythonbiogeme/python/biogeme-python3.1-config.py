import sys
from distutils import sysconfig
option=sys.argv[1]
if option=='includes':
 print('-I' + sysconfig.get_python_inc() + ' -I' + sysconfig.get_python_inc(plat_specific=True))
if option=='cxxflags':
 flags=sysconfig.get_config_var('CFLAGS').split()
 for unflag in ['-Wno-unused-result','-Wstrict-prototypes']:
  if flags.count(unflag) > 0:
   flags.remove(unflag)
 print(' '.join(flags))
if option=='libs':
 print('-lpython' + sysconfig.get_config_var('VERSION') + ' ' + sysconfig.get_config_var('LIBS') + ' ' + sysconfig.get_config_var('SYSLIBS'))
if option=='ldflags':
 print(sysconfig.get_config_var('LDFLAGS') + ' -L' + sysconfig.get_config_var('LIBPL'))
if option=='headerdir':
 print(sysconfig.get_python_inc())
if option=='version':
 print(sysconfig.get_config_var('VERSION'))
