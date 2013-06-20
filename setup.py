from sys import argv,executable,exit,stderr

if len(argv)<2 or len(argv)>3 or argv[1]!='install':
    print >>stderr, 'usage: '+executable+' '+argv[0]+' install [installation-directory]'
    exit(1)

if len(argv)==3:
    instdir = argv[2]
else:
    from distutils.sysconfig import get_python_lib
    instdir = get_python_lib()

from os import system
exit(system('make doinstall instdir='+instdir+' py='+executable))
