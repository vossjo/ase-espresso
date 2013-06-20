from sys import argv,executable,exit,stderr

if len(argv)<2 or len(argv)>3 or argv[1]!='install' or len(argv)==3 and argv[2][:9]!='--prefix=':
    print >>stderr, 'usage: '+executable+' '+argv[0]+' install [--prefix=installation-directory]'
    exit(1)

if len(argv)==3:
    instdir = argv[2][9:]
else:
    from distutils.sysconfig import get_python_lib
    instdir = get_python_lib()

from os import system
exit(system('make doinstall instdir='+instdir+' py='+executable))
